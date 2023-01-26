
using DrWatson
@quickactivate "metastability"
using GLMakie, OrdinaryDiffEq
include("$(scriptsdir())/utils.jl")


# ------------------------------------------------- Example on particle in double well (eg duffing). Following Ishii, 1986 (Physics Letters A) ------------------------------------------------ #
include("$(srcdir())/systems/duffing.jl")
include("$(srcdir())/paperplottheme.jl")
include("$(srcdir())/visualizations/plots.jl")
function plot_time_series()
    # T = 100000
    T = 1e6
    Ttr = 100000
    Δt=0.05
    tplot = 1000;

    parameters_ishii1986 = Dict(
        :a => 100, #b in their paper
        :b => 10, #-a in the paper
        :c => 0.0, #not present there (symmetric wells)
        :d => 1.0, # γ in their paper
        # :f => 0.853, #F
        :f => 0.852, #F
        # :f => 0.849, #F
        :ω => 3.5, #Ω
    )

    # Time series
    diffeq = (alg = Tsit5(), abstol=1e-8, reltol=1e-8, maxiters=1e8)
    u0 = [0.1, 0.1]
    cleft = :purple
    cright = :green
    fig, axs = subplotgrid(3, 1; resolution=(columnsize_pt, 2.0*columnsize_pt), xlabels=["t", "x", ""], ylabels=["x", "dx/dt", "PDF(τ)"], sharex=false, sharey=false)
    df = duffing_assymetric(parameters_ishii1986)
    tr = trajectory(df, T, u0; Ttr, Δt, diffeq); 
    ts = Ttr:Δt:Ttr+T; xs = tr[:,1]; ys = tr[:,2];
    time_idxs_plot = floor(Int64, T - tplot/Δt):1:ceil(Int64, T)
    ts_plot = ts[time_idxs_plot]; xs_plot = xs[time_idxs_plot]; ys_plot = ys[time_idxs_plot]
    idx_states_time_plot = [x >= 0 ? 1 : 0 for x in xs_plot]
    colors = replace(idx_states_time_plot, 1=>cright, 0=>cleft)

    ax = axs[1]
    lines!(ax, ts_plot, xs_plot, color=colors, linewidth=0.7)
    ax.xticks = [10000, 10500, 11000]
    # ylims!(ax, -0.6, 0.6)

    ax = axs[2]
    lines!(ax, xs_plot, ys_plot, color=colors)
    ax.xlabel = "x"
    rowsize!(fig.layout, 2, Relative(0.6))
    
    ax = axs[3] 
    ax.xlabel = "τ"
    idx_states_time = [x >= 0 ? 1 : 0 for x in xs]
    numbins = 10;
    for (i, key) in enumerate([0, 1])
    dwelltimes = length_samevalues_allowfluctuations(idx_states_time, 3)[1][key]
    weights, bins = histogram(dwelltimes, collect(range(0, 4000; length=numbins)))
    weights .+=1e-8
    lines!(ax, bins, weights, color=[cleft, cright][i])
    end
    ax.yscale = log10
    ax.xticks = 0:2000:4000
    save("$(plotsdir())/paper/intermittency-attractormergingcrisis.pdf", fig)
    fig
end

fig = plot_time_series()
fig