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

neigh_th=1e-3
prob, p = heteroclinic_cycle()
Ttr = 1000; Δt = 0.25; T = 1e6;
sol = solve(prob, Vern9(), maxiters=1e9, saveat=Ttr:Δt:T); ts = sol.t;
N = 3
logvars = @view sol[((1:N) .-1) .*2 .+ 1, :];
ogvars = map(x->z_to_s(x, p.Smax), logvars);
xs, ys, zs = ogvars[1, :], ogvars[2, :], ogvars[3, :]
fps = fixedpoints_ratemodel(p)
_, dwelltimes = dwelltimes_heteroclinicycle(xs, ys, zs, fps; neigh_th=0.001)
traj_state_idxs_all = classify_points_into_fps([xs ys zs]', fps[1:2:end, :]; neigh_th);
fp_colors, traj_colors = color_trajectory_hc(traj_state_idxs_all; c1, c2, c3, trajcolor=:black)

azim=7; elev=0.3



fig = Figure()
ax2 = Axis3(fig[1, 1]; azimuth=azim, elevation=elev, ylabel="s₂", zlabel="s₃", xlabel="s₁")
lines!(ax2, xs, ys, zs, color=traj_colors, linewidth=0.5)
scatter!(ax2, fps[1:2:end, :][:,1], fps[1:2:end, :][:, 2], fps[1:2:end, :][:, 3], color=[c1, c2, c3], markersize=6)
hc = [xs ys zs]
fig, ax = plot_arrow_following_data!(fig, ax2, (0.02, 0.05, 0.02), hc)
fig, ax = plot_arrow_following_data!(fig, ax2, (0.04, 0.02, 0.02), hc)
fig, ax = plot_arrow_following_data!(fig, ax2, (0.02, 0.02, 0.05), hc)

