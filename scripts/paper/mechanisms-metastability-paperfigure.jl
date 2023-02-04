using DrWatson
@quickactivate "metastability"
using GLMakie
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
fig = Figure(; resolution=(width, height))
idxcol = 1; idxrow = 1;
axs = []
# fig, _axs = noisybistable(generate_data; idxcol, fig); idxcol+=1; push!(axs, _axs);
# fig, _axs = heteroclinicycle(; idxcol, fig); idxcol+=1; push!(axs, _axs);
generate_data = true
# fig, _axs = amcrisis(generate_data; idxcol, fig); idxcol+=1; push!(axs, _axs);
fig, _axs = typeIintermittency(generate_data; idxcol, fig); idxcol+=1;push!(axs, _axs);
# fig, _axs = chaotic_saddle(; idxcol, fig); push!(axs, _axs);
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


# safesave("$(plotsdir())/paper/mechanismsmetastability-tmp.png", fig, px_per_unit=4)
using CairoMakie
safesave("$(plotsdir())/paper/mechanismsmetastability-tmp.pdf", fig)