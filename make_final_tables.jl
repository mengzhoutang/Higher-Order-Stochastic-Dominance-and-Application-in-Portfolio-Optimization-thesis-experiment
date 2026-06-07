using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using CSV, DataFrames, Printf

# =========================================================
# Table 1
# 来源：step9_table1_like.csv
# 保留真实资产名：NoDur, Manuf, BusEq, Shops, Hlth
# =========================================================
df1 = CSV.read("step9_table1_like.csv", DataFrame)

# 统一列名显示
rename!(df1, Dict(
    "RiskFunction" => "Risk Function"
))

# 数值保留 4 位小数
for c in names(df1)
    if c != "β"
        df1[!, c] = round.(Float64.(df1[!, c]), digits=4)
    end
end

CSV.write("Table1_final.csv", df1)
println("Saved: Table1_final.csv")

# =========================================================
# Table 3
# 来源：table3_like.csv
# 保留真实资产名：NoDur, Manuf, BusEq, Shops, Hlth
# =========================================================
df3 = CSV.read("table3_like.csv", DataFrame)

rename!(df3, Dict(
    "Item" => "asset",
    "high_order" => "Our approach",
    "post_kopa" => "Post and Kopa (2017)"
))

# 做成适合直接复制进 Word 的显示版
table3_display = DataFrame(
    asset = String[],
    our = String[],
    post = String[]
)

for i in 1:nrow(df3)
    item = String(df3[i, "asset"])
    x = Float64(df3[i, "Our approach"])
    y = Float64(df3[i, "Post and Kopa (2017)"])

    if item == "objective"
        push!(table3_display, (
            item,
            @sprintf("%.2f%%", 100x),
            @sprintf("%.2f%%", 100y)
        ))
    elseif item == "number of TSD constraints"
        push!(table3_display, (
            item,
            string(round(Int, x)),
            string(round(Int, y))
        ))
    else
        push!(table3_display, (
            item,
            @sprintf("%.4f", x),
            @sprintf("%.4f", y)
        ))
    end
end

rename!(table3_display, Dict(
    "our" => "Our approach",
    "post" => "Post and Kopa (2017)"
))

CSV.write("Table3_final.csv", table3_display)
println("Saved: Table3_final.csv")