using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using CSV, DataFrames, Statistics, Plots

default(
    dpi = 220,
    size = (1100, 700),
    legendfontsize = 10,
    guidefontsize = 11,
    tickfontsize = 9,
    titlefontsize = 13,
    linewidth = 1.5
)

function benchmark_mean(windowfile::String)
    df = CSV.read(windowfile, DataFrame)
    benchmark = Vector{Float64}(df[:, end])
    return mean(benchmark)
end

function cumulative_bands(weightsfile::String; bottom_to_top::Vector{String})
    w = CSV.read(weightsfile, DataFrame)

    x = 1:nrow(w)
    xtlbl = ["2nd", "3rd", "5th", "10th", "15th", "20th", "∞"]

    available_assets = String.(names(w)[2:end])
    for a in bottom_to_top
        if !(a in available_assets)
            error("Asset $a not found in weights file. Available assets: $(available_assets)")
        end
    end

    vals = Dict(a => 100 .* Vector{Float64}(w[:, a]) for a in available_assets)

    lower = Dict{String, Vector{Float64}}()
    upper = Dict{String, Vector{Float64}}()

    cumulative = zeros(length(x))
    for a in bottom_to_top
        lower[a] = copy(cumulative)
        cumulative .+= vals[a]
        upper[a] = copy(cumulative)
    end

    return x, xtlbl, lower, upper
end

function cumulative_area_plot(weightsfile::String; bottom_to_top::Vector{String})
    x, xtlbl, lower, upper = cumulative_bands(weightsfile; bottom_to_top = bottom_to_top)

    color_map = Dict(
        "NoDur" => :lightpink,
        "Manuf" => :lightgreen,
        "BusEq" => :thistle,
        "Shops" => :lightskyblue,
        "Hlth"  => :tan
    )

    p = plot(
        xlabel = "stochastic order",
        ylabel = "portfolio allocation (%)",
        xticks = (collect(x), xtlbl),
        ylims = (0, 100),
        legend = :outerright,
        legend_column = 1,
        grid = true,
        right_margin = 12Plots.mm
    )

    top_to_bottom = reverse(bottom_to_top)

    for a in top_to_bottom
        lab = (a == bottom_to_top[1]) ? string(a, " (bottom)") : a

        plot!(
            p,
            x,
            upper[a];
            fillrange = lower[a],
            fillalpha = 0.65,
            linealpha = 0.0,
            color = get(color_map, a, :gray),
            label = lab
        )
    end

    return p
end

function objective_plot(summaryfile::String, windowfile::String; figtitle::String="")
    s = CSV.read(summaryfile, DataFrame)

    x = 1:nrow(s)
    xtlbl = ["2nd", "3rd", "5th", "10th", "15th", "20th", "∞"]


    obj = 100 .* s.objective
    bm = 100 * benchmark_mean(windowfile)
    bm_line = fill(bm, nrow(s))

    p = plot(
        x, obj;
        marker = :circle,
        markersize = 4,
        linewidth = 1.5,
        xticks = (collect(x), xtlbl),
        xlabel = "stochastic order",
        ylabel = "portfolio return (%)",
        label = "Objective",
        title = figtitle,
        legend = :outerright,
        legend_column = 1,
        right_margin = 12Plots.mm
    )

    plot!(
        p,
        x,
        bm_line;
        linestyle = :dash,
        marker = :cross,
        markersize = 3,
        linewidth = 1.2,
        label = "Benchmark"
    )

    return p
end


function plot_figure6_style(summaryfile::String, weightsfile::String, windowfile::String, outfile::String;
                            figtitle::String="",
                            bottom_to_top::Vector{String})
    p1 = objective_plot(summaryfile, windowfile; figtitle = figtitle)
    p2 = cumulative_area_plot(weightsfile; bottom_to_top = bottom_to_top)

    p = plot(
        p1, p2;
        layout = @layout([a; b]),
        size = (1100, 760)
    )

    savefig(p, outfile)
    println("Saved figure to: ", outfile)
end


plot_figure6_style(
    "step8_summary_with_inf.csv",
    "step8_weights_with_inf.csv",
    "ff12_selected5_last22.csv",
    "fig6_style_22m_final.png";
    figtitle = "22-month experiment",
    bottom_to_top = ["Shops", "Hlth", "NoDur", "BusEq", "Manuf"]
)


plot_figure6_style(
    "step7_summary_with_inf.csv",
    "step7_weights_with_inf.csv",
    "ff12_selected5_last36.csv",
    "fig6_style_36m_final.png";
    figtitle = "36-month robustness experiment",
    bottom_to_top = ["Hlth", "Manuf", "Shops", "NoDur", "BusEq"]
)
