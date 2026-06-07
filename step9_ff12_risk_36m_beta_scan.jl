using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using CSV, DataFrames, Statistics, Printf, Random
using StochasticDominance

Random.seed!(1234)

# =========================
# 1. 参数设置
# =========================
sourcefile = joinpath(@__DIR__, "ff12_selected5_all_months.csv")
windowfile = joinpath(@__DIR__, "ff12_selected5_last36.csv")

n_keep = 36
SDorder = 2

# 论文 Table 1 风格：固定 second-order，比较不同 β
β_list = [0.1, 0.5, 0.8]

# 包接口需要额外指定 r
r = 2.0

# 数值参数
max_ipot_try = 300

# =========================
# 2. 读取全样本并截取最近 36 个月
# =========================
df_all = CSV.read(sourcefile, DataFrame)

if nrow(df_all) < n_keep
    error("样本不足 $(n_keep) 个时期，当前只有 $(nrow(df_all)) 行。")
end

df = last(df_all, n_keep)
CSV.write(windowfile, df)

# =========================
# 3. 固定按位置识别列
# =========================
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

println("Risk-function experiment settings:")
println("SDorder = ", SDorder)
println("β values = ", β_list)
println("r = ", r)
println()

# =========================
# 4. 转成优化器需要的格式
# =========================
R = Matrix{Float64}(df[:, asset_cols])   # scenarios × assets
ξ = permutedims(R)                       # assets × scenarios
ξ_0 = Vector{Float64}(df[:, benchmark_col])

n = size(R, 1)
p_ξ = fill(1 / n, n)
p_ξ_0 = fill(1 / n, n)

println("size(R) = ", size(R), "   (scenarios × assets)")
println("size(ξ) = ", size(ξ), "   (assets × scenarios)")
println("length(ξ_0) = ", length(ξ_0))
println()

# =========================
# 5. 小工具函数
# =========================
function pretty_weights(x, asset_cols)
    parts = String[]
    for (name, w) in zip(asset_cols, x)
        push!(parts, @sprintf("%s=%.4f", name, w))
    end
    return join(parts, ", ")
end

function try_optimize_one_beta(ξ, ξ_0, SDorder, p_ξ, p_ξ_0, β, r; max_ipot_try=100)
    try
        x_sol, q_sol, t_sol = optimize_min_riskreturn_SD(
            ξ, ξ_0, SDorder;
            p_ξ = p_ξ,
            p_ξ_0 = p_ξ_0,
            β = β,
            r = r,
            max_ipot = max_ipot_try,
            verbose = true
        )
        return x_sol, q_sol, t_sol, true
    catch e
        if e isa MethodError
            println("Current installed version does not accept max_ipot.")
            println("Falling back to the basic call.")
            println()

            x_sol, q_sol, t_sol = optimize_min_riskreturn_SD(
                ξ, ξ_0, SDorder;
                p_ξ = p_ξ,
                p_ξ_0 = p_ξ_0,
                β = β,
                r = r,
                verbose = true
            )
            return x_sol, q_sol, t_sol, false
        else
            rethrow(e)
        end
    end
end

# =========================
# 6. 主循环：固定 SDorder，扫描 β
# =========================
results = DataFrame(
    β = Float64[],
    risk_objective = Float64[],
    q_opt = Float64[],
    t_opt = Float64[],
    used_extended_kwargs = Bool[]
)

weight_store = Dict{Float64, Vector{Float64}}()

for β in β_list
    println("===== β = ", β, " =====")

    x_sol, q_sol, t_sol, used_extended = try_optimize_one_beta(
        ξ, ξ_0, SDorder, p_ξ, p_ξ_0, β, r;
        max_ipot_try = max_ipot_try
    )

    obj = riskfunction_asset_allocation(x_sol, q_sol, ξ, r, p_ξ, β)

    println("x_opt = ", x_sol)
    println("q_opt = ", q_sol)
    println("t_opt = ", t_sol)
    println("risk objective = ", obj)
    println("weights: ", pretty_weights(x_sol, asset_cols))
    println("used extended kwargs? ", used_extended)
    println()

    push!(results, (
        β = β,
        risk_objective = obj,
        q_opt = q_sol,
        t_opt = t_sol,
        used_extended_kwargs = used_extended
    ))

    weight_store[β] = x_sol
end

# =========================
# 7. 汇总输出
# =========================
println("========== SUMMARY ==========")
println(results)
println()

println("========== WEIGHTS BY BETA ==========")
for β in β_list
    println("β = ", β, ": ", pretty_weights(weight_store[β], asset_cols))
end
println()

best_idx = argmin(results.risk_objective)
println("Lowest risk objective among tested β values = ",
        results.risk_objective[best_idx],
        " at β = ", results.β[best_idx])

CSV.write(joinpath(@__DIR__, "step9_summary.csv"), results)

table1_like = DataFrame(β = β_list)

for (j, asset) in enumerate(asset_cols)
    table1_like[!, asset] = [weight_store[β][j] for β in β_list]
end

table1_like[!, :RiskFunction] = results.risk_objective

println("========== TABLE 1 LIKE ==========")
println(table1_like)
println()

CSV.write(joinpath(@__DIR__, "step9_table1_like.csv"), table1_like)