using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using CSV, DataFrames, Statistics, Printf, Random
using StochasticDominance

# 用于 ∞ 阶
using JuMP
using HiGHS

Random.seed!(1234)

sourcefile = joinpath(@__DIR__, "ff12_selected5_all_months.csv")
windowfile = joinpath(@__DIR__, "ff12_selected5_last22.csv")   

n_keep = 22                                                    
finite_orders = [2, 3, 5, 10, 15, 20]

eps_tol = 1e-8
max_iter_try = 800
n_particles_try = 600

df_all = CSV.read(sourcefile, DataFrame)
df = last(df_all, n_keep)
CSV.write(windowfile, df)

colnames = String.(names(df))
date_col = colnames[1]
asset_cols = colnames[2:end-1]
benchmark_col = colnames[end]

R = Matrix{Float64}(df[:, asset_cols])   # scenarios × assets
ξ = permutedims(R)                       # assets × scenarios
ξ_0 = Vector{Float64}(df[:, benchmark_col])

n = size(R, 1)
d = size(R, 2)

p_ξ = fill(1 / n, n)
p_ξ_0 = fill(1 / n, n)

function pretty_weights(x, asset_cols)
    join([@sprintf("%s=%.4f", name, w) for (name, w) in zip(asset_cols, x)], ", ")
end

function try_optimize_one_order(ξ, ξ_0, SDorder, p_ξ, p_ξ_0;
                                eps_tol=1e-8,
                                max_iter_try=800,
                                n_particles_try=600)
    try
        x_sol, t_sol = optimize_max_return_SD(
            ξ, ξ_0, SDorder;
            p_ξ = p_ξ,
            p_ξ_0 = p_ξ_0,
            ε = eps_tol,
            max_iter = max_iter_try,
            n_particles = n_particles_try,
            verbose = true
        )
        return x_sol, t_sol, true
    catch e
        if e isa MethodError
            x_sol, t_sol = optimize_max_return_SD(
                ξ, ξ_0, SDorder;
                p_ξ = p_ξ,
                p_ξ_0 = p_ξ_0,
                ε = eps_tol,
                verbose = true
            )
            return x_sol, t_sol, false
        else
            rethrow(e)
        end
    end
end

# ∞ 阶：maximize mean return subject to min scenario return >= min benchmark return
function optimize_infinity_order(R, ξ_0)
    n, d = size(R)
    bmin = minimum(ξ_0)

    model = Model(HiGHS.Optimizer)
    set_silent(model)

    @variable(model, x[1:d] >= 0)
    @constraint(model, sum(x[j] for j in 1:d) == 1)

    # 对每个 scenario，都要求组合收益 >= benchmark 的最小收益
    @constraint(model, [i in 1:n], sum(R[i, j] * x[j] for j in 1:d) >= bmin)

    μ = vec(mean(R, dims=1))
    @objective(model, Max, sum(μ[j] * x[j] for j in 1:d))

    optimize!(model)

    x_sol = value.(x)
    obj = objective_value(model)
    return x_sol, obj
end

results = DataFrame(
    OrderLabel = String[],
    OrderValue = Union{Missing, Float64}[],
    objective = Float64[],
    t_opt = Union{Missing, Float64}[],
    used_extended_kwargs = Union{Missing, Bool}[]
)

weight_store = Dict{String, Vector{Float64}}()

# 有限阶部分
for ord in finite_orders
    println("===== SDorder = ", ord, " =====")
    local x_sol, t_sol, used_extended

    x_sol, t_sol, used_extended = try_optimize_one_order(
        ξ, ξ_0, ord, p_ξ, p_ξ_0;
        eps_tol = eps_tol,
        max_iter_try = max_iter_try,
        n_particles_try = n_particles_try
    )

    obj = expected_portfolio_return(x_sol, ξ, p_ξ)

    println("x_opt = ", x_sol)
    println("t_opt = ", t_sol)
    println("objective = ", obj)
    println("weights: ", pretty_weights(x_sol, asset_cols))
    println()

    push!(results, (
        string(ord),
        float(ord),
        obj,
        t_sol,
        used_extended
    ))

    weight_store[string(ord)] = x_sol
end

# ∞ 阶部分
println("===== SDorder = ∞ =====")
x_inf, obj_inf = optimize_infinity_order(R, ξ_0)

println("x_inf = ", x_inf)
println("objective_inf = ", obj_inf)
println("weights: ", pretty_weights(x_inf, asset_cols))
println()

push!(results, (
    "∞",
    missing,
    obj_inf,
    missing,
    missing
))

weight_store["∞"] = x_inf

println("========== SUMMARY ==========")
println(results)

weight_df = DataFrame(OrderLabel = results.OrderLabel)
for (j, asset) in enumerate(asset_cols)
    weight_df[!, asset] = [weight_store[label][j] for label in results.OrderLabel]
end

CSV.write(joinpath(@__DIR__, "step8_summary_with_inf.csv"), results)   
CSV.write(joinpath(@__DIR__, "step8_weights_with_inf.csv"), weight_df)