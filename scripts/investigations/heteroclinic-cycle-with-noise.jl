using DrWatson 
@quickactivate
using GLMakie
using DifferentialEquations
include("$(srcdir())/visualizations/animations/animations.jl")
include("$(srcdir())/visualizations/plots.jl")
include("$(srcdir())/analyses/utils.jl")
include("$(srcdir())/attractor_classification/classify_fixedpoints.jl")
include("$(srcdir())/systems/ratemodel.jl")

function plot_heteroclinic_cycle(ts, xs, ys, zs, fps; fulltrajcolor, saddle_colors, title="", azimuth=7, elevation=0.4, mksize=12, resolution=(800, 600))
    fig = Figure(; resolution); axs = [];
    
    ax = Axis3(fig[1:3, 1]; azimuth, elevation); push!(axs, ax);
    lines!(ax, xs, ys, zs, color=fulltrajcolor) 
    scatter!(ax, fps[1:2:end, :][:,1], fps[1:2:end, :][:, 2], fps[1:2:end, :][:, 3], color=saddle_colors, markersize=mksize+4)
    hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
    ax.xticks = [0.0, 0.02, 0.04]
    ax.yticks = [0.0, 0.02, 0.04]
    ax.zticks = [0.0, 0.02, 0.04]
    
    for (idx, variable) in enumerate([xs, ys, zs])
        ax = Axis(fig[idx, 2]); push!(axs, ax)
        lines!(ax, ts, variable; color=fulltrajcolor)
    end
    
    supertitle(fig, title)
    return fig, axs 
end

function run_heteroclinic_cycle_with_noise(; noise_strength=0.0, Ttr=50000, T=100000, Δt=1.0)
    u0 = [-5, 0.031, -5, 0.03, -5, 0.03]
    hcgh, params = heteroclinic_cycle(; η=noise_strength, T, rngseed=1, u0)
    solver = noise_strength == 0.0 ? Vern9() : SKenCarp()
    # solver = noise_strength == 0.0 ? Vern9() : SOSRA()
    sol = solve(hcgh, solver, maxiters=1e9, saveat=Ttr:Δt:T, abstol=1e-9, reltol=1e-9);
    # sol = solve(hcgh, solver, maxiters=1e9, saveat=Ttr:Δt:T, abstol=1e-4, reltol=1e-4);
    N = 3
    logvars = @view sol[((1:N) .- 1) .* 2 .+ 1, :];
    ogvars = map(x->z_to_s(x, params.Smax), logvars);
    ts = sol.t; xs, ys, zs = collect(eachrow(ogvars))

    fps = fixedpoints_ratemodel(params)
    c1 = :green; c2 = :purple; c3=:red; fulltrajcolor=:black; saddle_colors = [c1, c2, c3]
    traj_state_idxs_all, dwelltimes = dwelltimes_heteroclinicycle(xs, ys, zs, fps; neigh_th=0.001)
    fp_colors, traj_colors = color_trajectory_hc(traj_state_idxs_all; c1, c2, c3, trajcolor=fulltrajcolor)
    fig, axs = plot_heteroclinic_cycle(ts, xs, ys, zs, fps; fulltrajcolor=traj_colors, saddle_colors, title="noise strength = $(noise_strength)")
end


#0.01 is way too high 


# fig, axs = run_heteroclinic_cycle_with_noise(noise_strength=0.001; Ttr=0.0)  #already pretty big
# fig, axs = run_heteroclinic_cycle_with_noise(noise_strength=0.0001; Ttr=0.0)  
# fig, axs = run_heteroclinic_cycle_with_noise(noise_strength=1e-3; Ttr=0.0, T=1e5)  
fig, axs = run_heteroclinic_cycle_with_noise(noise_strength=1e-2; Ttr=0.0, T=1e5)  
# for noise_strength in range(0.0, 0.1; step=0.01)
# for noise_strength in [1e-5, 1e-4, 1e-3]
#     Ttr = 0.0; T = 100000.0;
#     fig, axs = run_heteroclinic_cycle_with_noise(; noise_strength, Ttr, T) 
#     filename = "$(plotsdir())/mechanisms/heterocliniccycle/ratemodel/dynamics-ratemodel-noise_$noise_strength-T_$T-Ttr_$(Ttr).png"
#     mkpath(dirname(filename))
#     makiesave(filename, fig)
# end