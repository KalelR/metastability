using DrWatson
@quickactivate "metastability"
using GLMakie
include("$(srcdir())/paperplottheme.jl")
using DynamicalSystems, Statistics
using DifferentialEquations, CurveFit, DelimitedFiles
include("$(scriptsdir())/utils.jl")
include("$(srcdir())/visualizations/plots.jl")


set_theme!(figure3)


# Noisy Duffing 
include("$(srcdir())/visualizations/plots.jl")
include("$(scriptsdir())/paper/noisybistable.jl")
function noisybistable(; idxrow=1, idxcol=1, fig)
    prob, _ = duffing()
    T = 1e4; Ttr = 0; Δt = 0.5;
    sol = solve(prob, SKenCarp(), saveat=Ttr:Δt:T, maxiters=1e9, progress=true);
    distribution_info = get_distribution_noisy_bistable(sol)
    ts = sol.t; xs = sol[1,:]; xdots = sol[2,:];
    idxs_restricted = 1:floor(Int64, 5000/Δt)
    _fig, _axs = plot_noisy_bistable(ts[idxs_restricted], xs[idxs_restricted], xdots[idxs_restricted], distribution_info; fig, idxrow, idxcol, markersize=2)
    return _fig, _axs 
end


# Heteroclinic cycle
include("$(srcdir())/systems/ratemodel.jl")
function heteroclinicycle(; idxrow=1, idxcol=1, fig)
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

_fig, _axs = plot_heterocliniccycle_paper(ts[restricted_idxs], xs[restricted_idxs], ys[restricted_idxs], zs[restricted_idxs], p, dwelltimes, fps; c1, c2, c3, fig, idxrow, idxcol, markersize=12)
return _fig, _axs 
end

# rowsize!(fig.layout, 2, Relative(0.6))


# Attractor-merging crisis 
include("$(srcdir())/systems/duffing.jl")
function amcrisis(; idxrow=1, idxcol=1, fig)
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

    #first, get dwell times
    T = 1e6; Ttr = 1e4; Δt = 1.0; 
    prob, p = duffing(; a, b, c, d, f, ω, n₁=0, n₂=0, T, u0=[0.1, 0.1]) 
    sol = solve(prob, Tsit5(), saveat=Ttr:Δt:T, maxiters=1e9, progress=true, abstol=1e-8, reltol=1e-8);
    ts = sol.t; xs = sol[1, :]; ys = sol[2, :];
    dwelltimes = dwell_times_attractor_merging_crisis(xs)

    #now integrate at finer time scale to plot
    T = 1e4 + 1000; Ttr = 1e4; Δt = 0.1; 
    sol = solve(prob, Tsit5(), saveat=Ttr:Δt:T, maxiters=1e9, progress=true, abstol=1e-8, reltol=1e-8);
    ts = sol.t; xs = sol[1, :]; ys = sol[2, :];

    _fig, _axs = plot_attractormerging_crisis(ts, xs, ys, dwelltimes; fig, idxrow, idxcol, cleft=c2, cright=c1)
    return _fig, _axs 
end


# Type-I intermittency 
include("$(srcdir())/systems/lorenz.jl")
function typeIintermittency(; idxrow=1, idxcol=1, fig)
    σ = 10; β = 8/3; u0 = [0.1, 0.1, 0.1]

    Ttr = 780; T=850;
    # Ttr = 780; T=1e4; 
    Δt=0.01;
    ρ = 166.1; p = [σ, ρ, β]; #intermittency (CA);
    prob = ODEProblem(lorenz!, u0, (0, T), p)
    sol = solve(prob, Tsit5(), saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, maxiters=1e9);
    ts = sol.t; xs = sol[1, :]; ys = sol[2, :]; zs = sol[3,:]


    #get the limit cycle (before the bif) and save
    ρ = 166.06; p = [σ, ρ, β]; #LC, just before SN
    prob = ODEProblem(lorenz!, u0, (0, T), p); 

    limitcycle = solve(prob, Tsit5(), saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, maxiters=1e9);
    arewithin = mapslices(x->withinset(x, limitcycle', 2), sol[:, :], dims=1)[1, :];
    colors = [el == true ? :green : :purple for el in arewithin]
    dwelltimes = length_samevalues_allowfluctuations(arewithin, 3)[1] #1 corresponds to time in  the limit cycle (so 1=laminar)
    # dwelltimes = readdlm("$(datadir())/typeIintermittency-lorenz-dwelltimes.dat")[:, 1]

    # include("$(srcdir())/visualizations/plots.jl")
    # fig = Figure(; resolution=(width, height))
    # idxcol = 1; idxrow = 1;
    _fig, _axs = plot_typeI_intermittency(ts, xs, ys, zs, dwelltimes, colors; fig, idxrow, idxcol, c1, c2)
    return _fig, _axs 
end



#Chaotic saddle
include("$(srcdir())/systems/ikedamap.jl")
function chaotic_saddle(; idxcol=1, idxrow=1, fig)
    b = 0.9; c = 0.4; d = 6.0; a = 1.003;
    T = 2e5; Ttr = 100;
    u0 = [-2.0, -2.0]
    ik = Systems.ikedamap(u0;a, b, c, d)
    traj = trajectory(ik, T; Ttr); 
    ts = Ttr:Ttr+T; xs = traj[:,1]; ys = traj[:,2]; 

    #getting distribution of times in chaotic saddle for various ics
    #  icsperdim=50;
    icsperdim = 20
    τs = dwell_times_chaotic_saddle(ik, u0, icsperdim; T=1e6, Ttr=1e2, threshold=1)

    filename = "$(datadir())/durationinchaoticsaddle-d_0.2.dat";
    # filename = "$(datadir())/durationinchaoticsaddle-d_$d.dat";
    dwelltimes = readdlm(filename)[:,1]

    # xfit = bins; yfit = weights; A, B = CurveFit.exp_fit(xfit, yfit);

    #getting colors that distinguish chaotic saddle from fixed point
    fp = [2.97, 4.15]; dist_th = 3.0;
    bool_within_fp = map(x->iswithinneighborhood(x, fp, dist_th), traj);
    # bool_within_fp = map(x->iswithinneighborhood(x, [traj[1:1000, 1] traj[1:1000, 2]], 0.05), traj);
    colors = [el == true ? :purple : :green for el in bool_within_fp];


    # include("$(srcdir())/visualizations/plots.jl")
    # fig = Figure(; resolution=(width, height))
    _fig, _axs = plot_chaotic_saddle(ts, xs, ys, dwelltimes, colors; fig, idxrow, idxcol)
    return _fig, _axs 
end


c1 = :green; c2 = :purple; c3=:red
width = 2.0*columnsize_pt 
width = 4.0*columnsize_pt 
height = 1.0*width;
fig = Figure(; resolution=(width, height))
idxcol = 1; idxrow = 1;

fig, _axs = noisybistable(; idxcol, fig); idxcol+=1;
fig, _axs = heteroclinicycle(; idxcol, fig); idxcol+=1;
fig, _axs = amcrisis(; idxcol, fig); idxcol+=1;
fig, _axs = typeIintermittency(; idxcol, fig); idxcol+=1;
fig, _axs = chaotic_saddle(; idxcol, fig)
rowsize!(fig.layout, 2, Relative(0.6))
rowgap!(fig.layout, Relative(0.01))
colgap!(fig.layout, Relative(0.01))
