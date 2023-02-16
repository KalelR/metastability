using DrWatson
@quickactivate "metastability"
using CairoMakie
# using GLMakie
# GLMakie.activate!()
using DynamicalSystems, Statistics
using DifferentialEquations, CurveFit, DelimitedFiles
include("$(srcdir())/paperplottheme.jl");
include("$(scriptsdir())/utils.jl");

include("$(srcdir())/visualizations/plots.jl");
include("$(srcdir())/visualizations/paperplots.jl");
set_theme!(figure3);

generate_data = false
const c1 = :green; const c2 = :purple; const c3=:red
width = 2.5*columnsize_pt;
height = 1.0*columnsize_pt;
fig = Figure(; resolution=(width, height), figure_padding = 4)
idxcol = 1; idxrow = 1;
axs = []
fig, _axs = noisybistable(generate_data; idxcol, fig); idxcol+=1; push!(axs, _axs);
fig, _axs = heteroclinicycle(; idxcol, fig); idxcol+=1; push!(axs, _axs);
fig, _axs = amcrisis(generate_data; idxcol, fig); idxcol+=1; push!(axs, _axs);
fig, _axs = typeIintermittency(generate_data; idxcol, fig); idxcol+=1;push!(axs, _axs);
fig, _axs = chaotic_saddle(; idxcol, fig); push!(axs, _axs);
rowsize!(fig.layout, 2, Relative(0.6))
rowgap!(fig.layout, Relative(0.01))
colgap!(fig.layout, Relative(0.0))


for _axs in axs 
    _axs[1].yticks = WilkinsonTicks(2)
    _axs[1].ylabelpadding = 2.5
    _axs[3].ylabelpadding = 2.5
    _axs[1].yticksize = 3.0
    _axs[3].yticksize = 3.0
    _axs[1].xticksize = 3.0
    _axs[3].xticksize = 3.0
    _axs[1].xticks = WilkinsonTicks(2)
    _axs[3].xticks = WilkinsonTicks(2)
    hidedecorations!(_axs[2]; label=false)
    # hidespines!(_axs[2])
end

axs[1][1].yticks = [-5, 5]
axs[2][1].yticks = [0.0, 0.04]
axs[3][1].yticks = [-0.3, 0.3]
axs[3][1].xticks = [10200, 10800]
axs[4][1].yticks = [-50, 50]
axs[4][1].xticks = [800, 830]
ylims!(axs[4][1], [-50, 50])
axs[5][1].yticks = [0.0, 5.0]


fix_x_ticks_limits!(axs[4][1], [1990, 2020]; powers=false)
fix_y_ticks_limits!(axs[1][3], [1e-7, 1e-5, 1e-2])
fix_y_ticks_limits!(axs[2][3], [1e4, 1e5, 1e6])
fix_y_ticks_limits!(axs[3][3], [1e-6, 1e-4, 1e-1])
fix_y_ticks_limits!(axs[4][3], [1e-6, 1e-4, 1e-1])
fix_y_ticks_limits!(axs[5][3], [1e-8, 1e-6, 1e-5])

# safesave("$(plotsdir())/paper/mechanismsmetastability-tmp.png", fig, px_per_unit=4)
safesave("$(plotsdir())/paper/mechanismsmetastability.png", fig, px_per_unit=40)
using CairoMakie
safesave("$(plotsdir())/paper/mechanismsmetastability.pdf", fig)



# GLMakie.activate!()
# # for _axs in axs ?
# label = "A"
# label = "A\^"
#     Label(fig.layout[2,1][1, 1, TopLeft()], label,
#         textsize = 12,
#         font = :bold,
#         padding = (0, 5, 5, 0),
#         halign = :right)
# # end