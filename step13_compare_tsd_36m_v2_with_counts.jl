using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using CSV, DataFrames, Statistics, Printf, LinearAlgebra, Random
using StochasticDominance
using JuMP, Ipopt
import MathOptInterface as MOI

# =========================================================
# 0. 固定随机种子
# =========================================================
Random.seed!(1234)

# =========================================================
# 1. 参数设置
# =========================================================

sourcefile = joinpath(@__DIR__, "ff12_selected5_all_months.csv")
windowfile = joinpath(@__DIR__, "ff12_selected5_last36.csv")

n_keep = 36

# high_order 方法：先 dummy 跑 2，再正式跑 3
dummy_SDorder = 2
target_SDorder = 3

# high_order 数值参数
eps_tol = 1e-8
max_iter_try = 800
n_particles_try = 600

# 如果你之后想按论文 Table 3 的口径手工填 3，
# 就把 nothing 改成 3
high_order_constraints_manual = 3

# Post-Kopa baseline 参数
threshold_mode = :equal25
n_equal_grid = 25

# 求解器参数
ipopt_tol = 1e-8
ipopt_max_iter = 5000
ipopt_print_level = 0

# 后验 TSD 检查
validate_grid_n = 400
nontrivial_eps = 1e-6

# active TSD risk levels 数值容差
active_level_atol = 1e-8

# =========================================================
# 2. 读取全样本并截取最近 36 个月
# =========================================================

df_all = CSV.read(sourcefile, DataFrame)

if nrow(df_all) < n_keep
    error("样本不足 $(n_keep) 个时期，当前只有 $(nrow(df_all)) 行。")
end

df = last(df_all, n_keep)
CSV.write(windowfile, df)

colnames = String.(names(df))
if length(colnames) < 3
    error("列数不足，至少需要 Date + 若干资产列 + Benchmark。")
end

date_col = colnames[1]
asset_cols = colnames[2:end-1]
benchmark_col = colnames[end]

println("All columns: ", colnames)
println("Selected industries: ", asset_cols)
println("Date range: ", df[1, date_col], " to ", df[end, date_col])
println("Window length: ", nrow(df), " months")
println("Saved 36-month window to: ", windowfile)
println()

# R: scenarios × assets
R = Matrix{Float64}(df[:, asset_cols])

# ξ: assets × scenarios
ξ = permutedims(R)

# benchmark return series
ξ_0 = Vector{Float64}(df[:, benchmark_col])

T = size(R, 1)
K = size(R, 2)

# 等概率
p = fill(1.0 / T, T)

println("size(R) = ", size(R), "   (scenarios × assets)")
println("size(ξ) = ", size(ξ), "   (assets × scenarios)")
println("length(ξ_0) = ", length(ξ_0))
println()

# =========================================================
# 3. 工具函数
# =========================================================

function pretty_weights(x::AbstractVector, asset_cols)
    parts = String[]
    for (name, w) in zip(asset_cols, x)
        push!(parts, @sprintf("%s=%.4f", name, w))
    end
    return join(parts, ", ")
end

portfolio_returns(R::AbstractMatrix, w::AbstractVector) = vec(R * w)

mean_return(series::AbstractVector, p::AbstractVector) = dot(p, series)

function expected_shortfall_lpm1(series::AbstractVector, thr::Float64, p::AbstractVector)
    return dot(p, max.(thr .- series, 0.0))
end

function semivariance_lpm2(series::AbstractVector, thr::Float64, p::AbstractVector)
    return dot(p, max.(thr .- series, 0.0).^2)
end

function build_thresholds(y_sorted::AbstractVector; mode::Symbol=:equal25, n_equal_grid::Int=25)
    if mode == :benchmark
        return collect(y_sorted)
    elseif mode == :equal25
        lo, hi = minimum(y_sorted), maximum(y_sorted)
        return collect(range(lo, hi; length=n_equal_grid))
    else
        error("unknown threshold_mode = $mode")
    end
end

function compute_postkopa_tolerances(y_sorted::AbstractVector, thresholds::AbstractVector, p::AbstractVector)
    m = length(thresholds)
    E_b = [expected_shortfall_lpm1(y_sorted, thresholds[s], p) for s in 1:m]
    S2_b = [semivariance_lpm2(y_sorted, thresholds[s], p) for s in 1:m]

    delta = zeros(m)
    if m >= 1
        delta[1] = 0.0
    end
    if m >= 2
        delta[2] = 0.0
    end
    for s in 3:m
        denom = S2_b[s-1] + 2.0 * E_b[s-1] * (thresholds[s] - thresholds[s-1])
        if denom <= 1e-14
            delta[s] = 0.0
        else
            delta[s] = max(S2_b[s] / denom - 1.0, 0.0)
        end
    end
    return E_b, S2_b, delta
end

function dense_validation_grid(bench::AbstractVector, port1::AbstractVector, port2::AbstractVector; n::Int=400)
    lo = minimum(vcat(bench, port1, port2))
    hi = maximum(vcat(bench, port1, port2))
    if isapprox(lo, hi; atol=1e-12)
        return [lo]
    end
    return collect(range(lo, hi; length=n))
end

function min_tsd_gap_on_grid(bench::AbstractVector, port::AbstractVector, p::AbstractVector, grid::AbstractVector)
    gaps = [semivariance_lpm2(bench, t, p) - semivariance_lpm2(port, t, p) for t in grid]
    return minimum(gaps), gaps
end

function validation_metrics(bench::AbstractVector, port::AbstractVector, p::AbstractVector, grid::AbstractVector; eps_left::Float64=1e-6)
    min_gap, gaps = min_tsd_gap_on_grid(bench, port, p, grid)
    max_violation = max(0.0, -min_gap)

    tol_left = minimum(bench) + eps_left
    idx_nontrivial = findall(t -> t > tol_left, grid)

    if isempty(idx_nontrivial)
        min_nontrivial_gap = min_gap
    else
        min_nontrivial_gap = minimum(gaps[idx_nontrivial])
    end

    return (
        min_gap = min_gap,
        max_violation = max_violation,
        min_nontrivial_gap = min_nontrivial_gap,
        gaps = gaps,
    )
end

"""
提取当前 high_order 解下的 active TSD risk levels.

定义方式：
- 对 sorted(ξ_0) 的每个 threshold t，计算 g_p(t, x, ...)
- 找到达到最大值 g_bar 的那些 t
- 去重后返回其个数和具体值

这里统计的是“当前 GitHub 实现口径下的活跃风险水平个数”，
不是论文 Table 3 模型层面的固定常数。
"""
function extract_active_tsd_levels(x, ξ, ξ_0, SDorder, p_ξ, p_ξ_0; atol=1e-8)
    p = SDorder - 1.0
    sorted_ξ0 = sort(ξ_0)

    vals = [StochasticDominance.g_p(t, x, ξ, ξ_0, p, p_ξ, p_ξ_0) for t in sorted_ξ0]
    maxval = maximum(vals)

    active_idx = findall(v -> isapprox(v, maxval; atol=atol, rtol=0.0), vals)
    active_t = sorted_ξ0[active_idx]
    active_t_unique = unique(active_t)

    return (
        max_gbar = maxval,
        active_indices = active_idx,
        active_t = active_t,
        active_t_unique = active_t_unique,
        active_count = length(active_t_unique),
        all_t = sorted_ξ0,
        all_vals = vals,
    )
end

function active_levels_to_string(v::AbstractVector)
    if isempty(v)
        return ""
    end
    return join([@sprintf("%.10f", x) for x in v], "; ")
end

# =========================================================
# 4. high_order 方法
# =========================================================

function try_optimize_high_order_order(ξ, ξ_0, SDorder, p;
    eps_tol=1e-8,
    max_iter_try=800,
    n_particles_try=600)

    try
        x_sol, t_sol = optimize_max_return_SD(
            ξ,
            ξ_0,
            SDorder;
            p_ξ = p,
            p_ξ_0 = p,
            ε = eps_tol,
            max_iter = max_iter_try,
            n_particles = n_particles_try,
            verbose = true
        )
        return x_sol, t_sol, true
    catch e
        if e isa MethodError
            println("Installed StochasticDominance version does not accept max_iter / n_particles.")
            println("Falling back to the basic call with ε = ", eps_tol)
            println()

            x_sol, t_sol = optimize_max_return_SD(
                ξ,
                ξ_0,
                SDorder;
                p_ξ = p,
                p_ξ_0 = p,
                ε = eps_tol,
                verbose = true
            )
            return x_sol, t_sol, false
        else
            rethrow(e)
        end
    end
end

# =========================================================
# 5. Post-Kopa baseline（SCTSD / convex QCP）
# =========================================================

function solve_post_kopa_baseline(
    R::Matrix{Float64},
    ξ0::Vector{Float64},
    p::Vector{Float64};
    threshold_mode::Symbol = :equal25,
    n_equal_grid::Int = 25,
    ipopt_tol::Float64 = 1e-8,
    ipopt_max_iter::Int = 5000,
    ipopt_print_level::Int = 0
)
    T, K = size(R)

    perm = sortperm(ξ0)
    R_sorted = R[perm, :]
    y_sorted = ξ0[perm]
    p_sorted = p[perm]

    thresholds = build_thresholds(y_sorted; mode=threshold_mode, n_equal_grid=n_equal_grid)
    m = length(thresholds)

    E_b, S2_b, delta = compute_postkopa_tolerances(y_sorted, thresholds, p_sorted)
    μ_bench = mean_return(y_sorted, p_sorted)

    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "tol", ipopt_tol)
    set_optimizer_attribute(model, "max_iter", ipopt_max_iter)
    set_optimizer_attribute(model, "print_level", ipopt_print_level)

    @variable(model, w[1:K] >= 0)
    @variable(model, z[1:m, 1:T] >= 0)

    # 初值：等权
    w0 = fill(1.0 / K, K)
    ret0 = R_sorted * w0

    for k in 1:K
        set_start_value(w[k], w0[k])
    end
    for s in 1:m, t in 1:T
        set_start_value(z[s, t], max(thresholds[s] - ret0[t], 0.0))
    end

    @constraint(model, sum(w[k] for k in 1:K) == 1)

    @expression(model, portret[t=1:T], sum(R_sorted[t, k] * w[k] for k in 1:K))
    @expression(model, mean_port, sum(p_sorted[t] * portret[t] for t in 1:T))

    @constraint(model, mean_port >= μ_bench)
    @constraint(model, [s=1:m, t=1:T], z[s, t] + portret[t] >= thresholds[s])
    @constraint(model, [s=1:m], (1.0 + delta[s]) * sum(p_sorted[t] * z[s, t]^2 for t in 1:T) <= S2_b[s])

    @objective(model, Max, mean_port)

    optimize!(model)

    term = termination_status(model)
    pstat = primal_status(model)

    ok = term in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED, MOI.ALMOST_OPTIMAL)
    if !ok
        @warn "Post-Kopa baseline did not end with OPTIMAL/LOCALLY_SOLVED. termination_status = $term, primal_status = $pstat"
    end

    w_sol = value.(w)
    obj_sol = objective_value(model)

    model_tsd_constraint_count = m*T + m + 2

    return (
        x = w_sol,
        objective = obj_sol,
        status = term,
        primal_status = pstat,
        thresholds = thresholds,
        delta = delta,
        constraint_count = model_tsd_constraint_count,
        sorted_R = R_sorted,
        sorted_benchmark = y_sorted,
        sorted_p = p_sorted,
    )
end

# =========================================================
# 6. 先 dummy 跑 order 2，再正式跑 order 3
# =========================================================

println("===== OUR METHOD: DUMMY RUN (order = $(dummy_SDorder)) =====")
x_dummy, t_dummy, used_extended_dummy = try_optimize_high_order_order(
    ξ, ξ_0, dummy_SDorder, p;
    eps_tol = eps_tol,
    max_iter_try = max_iter_try,
    n_particles_try = n_particles_try
)
obj_dummy = expected_portfolio_return(x_dummy, ξ, p)
println("x_dummy = ", x_dummy)
println("t_dummy = ", t_dummy)
println("obj_dummy = ", obj_dummy)
println("weights_dummy = ", pretty_weights(x_dummy, asset_cols))
println("used extended kwargs? ", used_extended_dummy)
println()

println("===== OUR METHOD: FORMAL RUN (order = $(target_SDorder)) =====")
x_ours, t_ours, used_extended = try_optimize_high_order_order(
    ξ, ξ_0, target_SDorder, p;
    eps_tol = eps_tol,
    max_iter_try = max_iter_try,
    n_particles_try = n_particles_try
)
obj_ours = expected_portfolio_return(x_ours, ξ, p)
println("x_ours = ", x_ours)
println("t_ours = ", t_ours)
println("obj_ours = ", obj_ours)
println("weights_ours = ", pretty_weights(x_ours, asset_cols))
println("used extended kwargs? ", used_extended)
println()

println("===== ACTIVE TSD RISK LEVELS FOR high_order =====")
tsd_info = extract_active_tsd_levels(x_ours, ξ, ξ_0, target_SDorder, p, p; atol=active_level_atol)
println("max g_bar value = ", tsd_info.max_gbar)
println("active TSD risk levels count = ", tsd_info.active_count)
println("active TSD risk levels = ", tsd_info.active_t_unique)
println()

println("===== POST-KOPA BASELINE (SCTSD / QCP, threshold_mode = :equal25) =====")
pk = solve_post_kopa_baseline(
    R, ξ_0, p;
    threshold_mode = threshold_mode,
    n_equal_grid = n_equal_grid,
    ipopt_tol = ipopt_tol,
    ipopt_max_iter = ipopt_max_iter,
    ipopt_print_level = ipopt_print_level
)

x_pk = pk.x
obj_pk = pk.objective
println("x_pk = ", x_pk)
println("obj_pk = ", obj_pk)
println("weights_pk = ", pretty_weights(x_pk, asset_cols))
println("pk solver status = ", pk.status)
println("pk constraint count (excluding nonnegativity) = ", pk.constraint_count)
println("number of thresholds used = ", length(pk.thresholds))
println()

# =========================================================
# 7. 后验 TSD 检查
# =========================================================

perm = sortperm(ξ_0)
R_sorted = R[perm, :]
ξ0_sorted = ξ_0[perm]
p_sorted = p[perm]

ret_ours = portfolio_returns(R_sorted, x_ours)
ret_pk   = portfolio_returns(R_sorted, x_pk)

grid = dense_validation_grid(ξ0_sorted, ret_ours, ret_pk; n=validate_grid_n)

val_ours = validation_metrics(ξ0_sorted, ret_ours, p_sorted, grid; eps_left=nontrivial_eps)
val_pk   = validation_metrics(ξ0_sorted, ret_pk,   p_sorted, grid; eps_left=nontrivial_eps)

println("===== DENSE-GRID TSD VALIDATION =====")
println("OURS:")
println("  min gap              = ", val_ours.min_gap)
println("  max violation        = ", val_ours.max_violation)
println("  min nontrivial gap   = ", val_ours.min_nontrivial_gap)
println()
println("POST-KOPA:")
println("  min gap              = ", val_pk.min_gap)
println("  max violation        = ", val_pk.max_violation)
println("  min nontrivial gap   = ", val_pk.min_nontrivial_gap)
println()

# =========================================================
# 8. 汇总输出
# =========================================================

results = DataFrame(
    Method = String[],
    Objective = Float64[],
    AuxParam = String[],
    TSD_Constraint_Count = Union{Missing, Int}[],
    Active_TSD_Risk_Level_Count = Union{Missing, Int}[],
    Active_TSD_Risk_Levels = Union{Missing, String}[],
    Max_TSD_Violation = Float64[],
    Min_Nontrivial_TSD_Gap = Float64[],
    SolverStatus = String[]
)

ho_constraints = isnothing(high_order_constraints_manual) ? 2 + tsd_info.active_count : high_order_constraints_manual

push!(results, (
    Method = "high_order_after_dummy_order2",
    Objective = obj_ours,
    AuxParam = @sprintf("t_opt=%.8f", t_ours),
    TSD_Constraint_Count = ho_constraints,
    Active_TSD_Risk_Level_Count = tsd_info.active_count,
    Active_TSD_Risk_Levels = active_levels_to_string(tsd_info.active_t_unique),
    Max_TSD_Violation = val_ours.max_violation,
    Min_Nontrivial_TSD_Gap = val_ours.min_nontrivial_gap,
    SolverStatus = "package_internal"
))

push!(results, (
    Method = "PostKopa_SCTSD_threshold=equal25",
    Objective = obj_pk,
    AuxParam = @sprintf("m=%d", length(pk.thresholds)),
    TSD_Constraint_Count = pk.constraint_count,
    Active_TSD_Risk_Level_Count = missing,
    Active_TSD_Risk_Levels = missing,
    Max_TSD_Violation = val_pk.max_violation,
    Min_Nontrivial_TSD_Gap = val_pk.min_nontrivial_gap,
    SolverStatus = string(pk.status)
))

println("========== COMPARISON SUMMARY ==========")
println(results)
println()

println("========== WEIGHTS ==========")
println("high_order: ", pretty_weights(x_ours, asset_cols))
println("post_kopa: ", pretty_weights(x_pk, asset_cols))
println()

CSV.write(joinpath(@__DIR__, "compare_tsd_36m_v2_with_counts_summary.csv"), results)
CSV.write(
    joinpath(@__DIR__, "compare_tsd_36m_v2_with_counts_weights.csv"),
    DataFrame(
        Asset = asset_cols,
        high_order = x_ours,
        post_kopa = x_pk
    )
)

table3_like = DataFrame(
    Item = vcat(asset_cols, ["objective", "number of TSD constraints"]),
    high_order = vcat(x_ours, [obj_ours], [Float64(ho_constraints)]),
    post_kopa = vcat(x_pk, [obj_pk], [Float64(pk.constraint_count)])
)

println("========== TABLE 3 LIKE ==========")
println(table3_like)
println()

CSV.write(joinpath(@__DIR__, "table3_like.csv"), table3_like)