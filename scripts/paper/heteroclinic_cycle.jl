using DrWatson
@quickactivate "metastability"
using GLMakie, OrdinaryDiffEq, CurveFit
include("$(scriptsdir())/utils.jl")
include("$(srcdir())/paperplottheme.jl")
include("$(srcdir())/attractor_classification/classify_fixedpoints.jl")
include("$(srcdir())/visualizations/plots.jl")
include("$(srcdir())/visualizations/animations.jl")

# -------------------------------- Rate model -------------------------------- #
include("$(srcdir())/systems/ratemodel.jl")


azim = 6.995530633326985; elev = 0.3926990816987241;
fig, axs = plot_heterocliniccycle_paper(t_plot, ss_plot, traj_colors, dwelltimes, fp_colors)
axs[1].xticks = 0:20000:tplot;
axs[2].xticks = [0.02, 0.04]; axs[2].yticks = 0:0.02:0.04; axs[2].zticks = 0:0.02:0.04;
save("$(plotsdir())/paper/heterocliniccycle.pdf", fig)