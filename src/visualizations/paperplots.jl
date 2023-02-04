include("$(srcdir())/visualizations/plots.jl")


include("$(scriptsdir())/paper/noisybistable.jl")
function noisybistable(generate_data=false; idxrow=1, idxcol=1, fig)
    @info "Running and plotting noisy bistable"
    prob, params = duffing()
    T = 1e8; Ttr = 0; Δt = 0.5;
    sol = solve(prob, SKenCarp(), saveat=Ttr:Δt:T, maxiters=1e9, progress=true, seed=2);
    n₁ = params[7]; n₂ = params[8]
    filename = "$(datadir())/noisybistable-distributioninfo-n1_$(n₁)-n2_$(n₂)-T_$T.jld2"
    distribution_data, file = produce_or_load(Dict("v"=>sol[1,:]), distribution_info_noisy_bistable; filename, force=generate_data, tag=false)
    ts = sol.t; xs = sol[1,:]; xdots = sol[2,:];
    idxs_restricted = 1:floor(Int64, 5000/Δt)
    _fig, _axs = plot_noisy_bistable(ts[idxs_restricted], xs[idxs_restricted], xdots[idxs_restricted], distribution_data; fig, idxrow, idxcol)
    return _fig, _axs 
end


include("$(srcdir())/systems/ratemodel.jl")
function heteroclinicycle(; idxrow=1, idxcol=1, fig)
    @info "Running and plotting heteroclinic cycle"
    prob, p = heteroclinic_cycle()
    Ttr = 1000; Δt = 0.25; T = 1e6;
    sol = solve(prob, Vern9(), maxiters=1e9, saveat=Ttr:Δt:T); ts = sol.t;
    N = 3
    logvars = @view sol[((1:N) .-1) .*2 .+ 1, :];
    ogvars = map(x->z_to_s(x, p.Smax), logvars);
    xs, ys, zs = ogvars[1, :], ogvars[2, :], ogvars[3, :]
    fps = fixedpoints_ratemodel(p)
    _, dwelltimes = dwelltimes_heteroclinicycle(xs, ys, zs, fps; neigh_th=0.001)

    restricted_idxs= 1:4*35000

    _fig, _axs = plot_heterocliniccycle_paper(ts[restricted_idxs], xs[restricted_idxs], ys[restricted_idxs], zs[restricted_idxs], p, dwelltimes, fps; c1, c2=c3, c3=c2, fig, idxrow, idxcol)
    return _fig, _axs 
end

include("$(srcdir())/systems/duffing.jl")
function amcrisis(generate_data=false; idxrow=1, idxcol=1, fig)
    @info "Running and plotting attractor-merging crisis"
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
    filename = "$(datadir())/attractormergingcrisis-dwelltimesinfo-T_$T.jld2"
    dwelltimes_data, file = produce_or_load(Dict("v"=>xs), dwell_times_attractor_merging_crisis; filename, force=generate_data, tag=false)

    #now integrate at finer time scale to plot
    T = 1e4 + 1000; Ttr = 1e4; Δt = 0.1; 
    sol = solve(prob, Tsit5(), saveat=Ttr:Δt:T, maxiters=1e9, progress=true, abstol=1e-8, reltol=1e-8);
    ts = sol.t; xs = sol[1, :]; ys = sol[2, :];

    _fig, _axs = plot_attractormerging_crisis(ts, xs, ys, dwelltimes_data; fig, idxrow, idxcol, cleft=c2, cright=c1)
    return _fig, _axs 
end

include("$(srcdir())/systems/lorenz.jl")
function typeIintermittency(generate_data=true; idxrow=1, idxcol=1, fig)
    @info "Running and plotting type I intermittency"
    σ = 10; β = 8/3; u0 = [0.1, 0.1, 0.1]
    ρ_chaos = 166.1; p_chaos = [σ, ρ_chaos, β]; #intermittency (CA);
    ρ_lc = 166.06; p_lc = [σ, ρ_lc, β]; #LC, just before SN
    
    #dwell times
    T = 1e6
    dwelltimes_input = Dict("T"=>T, "p_chaos"=>p_chaos, "system_func"=>lorenz!, "u0"=>u0, "p_lc"=>p_lc)
    filename = "$(datadir())/typeIintermittency-dwelltimesinfo-T_$T.jld2"
    dwelltimes_data, file = produce_or_load(dwelltimes_input, dwelltimes_info_laminar_phase; filename, force=generate_data, tag=false)
    @unpack dwelltimes, sol_lc = dwelltimes_data

    #to plot
    Ttr = 780; T=850; Δt=0.01; 
    prob = ODEProblem(lorenz!, u0, (0, T), p_chaos)
    sol = solve(prob, Tsit5(), saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, maxiters=1e9);
    ts = sol.t; xs = sol[1, :]; ys = sol[2, :]; zs = sol[3,:]
    arewithin = mapslices(x->withinset(x, sol_lc', 2), sol[:, :], dims=1)[1, :];
    colors = [el == true ? c1 : c2 for el in arewithin]
    
    _fig, _axs = plot_typeI_intermittency(ts, xs, ys, zs, dwelltimes, colors; fig, idxrow, idxcol, c1, c2)
    return _fig, _axs 
end

include("$(srcdir())/systems/ikedamap.jl")
function chaotic_saddle(; idxcol=1, idxrow=1, fig)
    @info "Running and plotting chaotic saddle"
    b = 0.9; c = 0.4; d = 6.0; a = 1.003;
    T = 2e5; Ttr = 100;
    u0 = [-2.0, -2.0]
    ik = Systems.ikedamap(u0;a, b, c, d)
    traj = trajectory(ik, T; Ttr); 
    ts = Ttr:Ttr+T; xs = traj[:,1]; ys = traj[:,2]; 

    #getting distribution of times in chaotic saddle for various ics
    #  icsperdim=50;
    icsperdim = 20
    fp = [2.97, 4.15]; dist_th = 3.0;
    # τs = dwell_times_chaotic_saddle(ik, u0, fp; icsperdim; T=1e6, Ttr=1e2, threshold=1)

    filename = "$(datadir())/durationinchaoticsaddle-d_0.2.dat";
    # filename = "$(datadir())/durationinchaoticsaddle-d_$d.dat";
    dwelltimes = readdlm(filename)[:,1]

    # xfit = bins; yfit = weights; A, B = CurveFit.exp_fit(xfit, yfit);

    #getting colors that distinguish chaotic saddle from fixed point
    bool_within_fp = map(x->iswithinneighborhood(x, fp, dist_th), traj);
    # bool_within_fp = map(x->iswithinneighborhood(x, [traj[1:1000, 1] traj[1:1000, 2]], 0.05), traj);
    colors = [el == true ? c1 : c2 for el in bool_within_fp];


    # include("$(srcdir())/visualizations/plots.jl")
    # fig = Figure(; resolution=(width, height))
    _fig, _axs = plot_chaotic_saddle(ts, xs, ys, dwelltimes, colors; fig, idxrow, idxcol, color_pdf=c2)
    return _fig, _axs 
end

