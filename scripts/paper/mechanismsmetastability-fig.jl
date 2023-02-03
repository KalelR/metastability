using DrWatson
@quickactivate "metastability"
using GLMakie
include("$(srcdir())/paperplottheme.jl")
using DynamicalSystems, Statistics
include("$(scriptsdir())/utils.jl")
include("$(srcdir())/visualizations/plots.jl")


c1 = :green; c2 = :purple; c3=:red
width = 2.0*columnsize_pt 
height = 1.0*width;
fig = Figure(; resolution=(width, height))
idxcol = 1; idxrow = 1;

# Noisy Duffing 
include("$(scriptsdir())/paper/noisybistable.jl")
prob, _ = duffing()
T = 1e4; Ttr = 0; Δt = 0.5;
sol = solve(prob, SKenCarp(), saveat=Ttr:Δt:T, maxiters=1e9, progress=true);
distribution_info = get_distribution_noisy_bistable(sol)
_fig, _axs = plot_noisy_bistable(sol, distribution_info; fig, idxrow, idxcol, markersize=2)
idxcol +=1; idxrow = 1;


# Heteroclinic cycle
include("$(srcdir())/systems/ratemodel.jl")
include("$(srcdir())/visualizations/plots.jl")

prob, p = heteroclinic_cycle()
Ttr = 1000; Δt = 1.0; T = 8e5;
sol = solve(prob, Vern9(), maxiters=1e9, saveat=Ttr:Δt:T); ts = sol.t;
N = 3
logvars = @view sol[((1:N) .-1) .*2 .+ 1, :];
ogvars = map(x->z_to_s(x, p.Smax), logvars);
xs, ys, zs = ogvars[1, :], ogvars[2, :], ogvars[3, :]
fps = fixedpoints_ratemodel(p)
_, dwelltimes = dwelltimes_heteroclinicycle(xs, ys, zs, fps; neigh_th=0.001)

restricted_idxs= 1:35000

include("$(srcdir())/visualizations/plots.jl")
fig = Figure(; resolution=(width, height))
idxcol = 1; idxrow = 1;
_fig, _axs = plot_heterocliniccycle_paper(ts[restricted_idxs], xs[restricted_idxs], ys[restricted_idxs], zs[restricted_idxs], p, dwelltimes, fps; c1, c2, c3, fig, markersize=12)

# rowsize!(fig.layout, 2, Relative(0.6))


# Attractor-merging crisis 
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

@unpack a, b, c, d, f, ω = parameters_ishii1986
prob, p = duffing(; a, b, c, d, f, ω, n₁=0, n₂=0, T, u0=[0.1, 0.1]) 
T = 1e6; Ttr = 100000; Δt = 0.05; 
sol = solve(prob, Tsit5(), saveat=Ttr:Δt:T, maxiters=1e9, progress=true, abstol=1e-8, reltol=1e-8);
ts = sol.t; xs = sol[1, :]; ys = tr[2, :];

tplot = 1000;
restricted_idxs = floor(Int64, (T - tplot)/Δt):1:ceil(Int64, T)

_fig, _axs = plot_attractormerging_crisis(ts[restricted_idxs], xs[restricted_idxs], ys[restricted_idxs]; fig, idxrow, idxcol, cleft=c2, cright=c1)
# rowsize!(fig.layout, 2, Relative(0.6))