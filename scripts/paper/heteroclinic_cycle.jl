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

p = rateparams()
T = 1e6; Ttr = 1000; Δt = 1.0;
u0 = [0.5, 0.2, 0.4, 0.9, 0.5, 0.6];
tspan = (0, T);
hcgh = ODEProblem(ratemodel_rule_Z!, u0, tspan, p);
sol = solve(hcgh, Vern9(), maxiters=1e9, saveat=Ttr:Δt:T); ts = sol.t;
zs = @view sol[((1:N) .-1) .*2 .+ 1, :]; ss = map(x->z_to_s(x, Smax), zs);

tplot = 35000; ss_plot = ss[:, 1:tplot]; t_plot = sol.t[1:tplot];

fps = fixedpoints_ratemodel(p)

traj_state_idxs_all, dwelltimes = dwelltimes_heteroclinicycle(ss[1,:], ss[2,:], ss[3,:], fps; neigh_th=0.001)
fp_colors, traj_colors = color_trajectory_hc(traj_state_idxs_all[1:tplot])

azim = 6.995530633326985; elev = 0.3926990816987241;
fig, axs = plot_heterocliniccycle_paper(t_plot, ss_plot, traj_colors, dwelltimes, fp_colors)
axs[1].xticks = 0:20000:tplot;
axs[2].xticks = [0.02, 0.04]; axs[2].yticks = 0:0.02:0.04; axs[2].zticks = 0:0.02:0.04;
save("$(plotsdir())/paper/heterocliniccycle.pdf", fig)