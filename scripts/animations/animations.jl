using DrWatson 
@quickactivate
using GLMakie
using DifferentialEquations
include("$(srcdir())/visualizations/animations/animations.jl")
include("$(srcdir())/visualizations/plots.jl")
include("$(srcdir())/analyses/utils.jl")
include("$(srcdir())/attractor_classification/classify_fixedpoints.jl")

set_theme!(animationtheme)

BALL_COLOR = :orange 
METASTABLE_STATE_1_COLOR = :green
METASTABLE_STATE_2_COLOR = :purple
METASTABLE_STATE_3_COLOR = :red

include("$(srcdir())/systems/ratemodel.jl")
function run_animation_shc()
    Ttr = 50000; T = 100000; Δt= 1.0;
    p = rateparams(); N = 3;
    u0 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    hcgh = ODEProblem(ratemodel_rule_Z!, u0, (0, T), p);
    sol = solve(hcgh, Vern9(), maxiters=1e9, saveat=Ttr:Δt:T);
    logvars = @view sol[((1:N) .-1) .*2 .+ 1, :];
    ogvars =map(x->z_to_s(x, p.Smax), logvars);

    fps = fixedpoints_ratemodel(p)
    az = 6.995530633326985; el = 0.3926990816987241
    framerate = 25; numframes = 400
    filename = "$(plotsdir())/mechanisms/heterocliniccycle/ratemodel/animation-ratemodel-logvars-T_$T-Ttr_$(Ttr)-framerate_$framerate.mp4"
    ts = sol.t .- sol.t[1] #remove transient (improve visualization in animation)
    animate_heterolinic_cycle(filename, ts, ogvars[1, :], ogvars[2, :], ogvars[3, :], fps; framerate, numframes, az, el, ball_color = BALL_COLOR, c1=METASTABLE_STATE_1_COLOR, c2=METASTABLE_STATE_2_COLOR, c3=METASTABLE_STATE_3_COLOR, resolution=_RESOLUTION)
end



include("$(srcdir())/systems/duffing.jl")
function run_animation_doublewell()
    T = 150; Ttr = 0; Δt = 0.1
    u0 = [-3.0, 0]
    
    n₁ = n₂ = n = 0.25; 
    rngseed = 5 #one transition in the middle
    rngseed = 7 #one transition in the middle
    prob_duffing, params = duffing(; n₁, n₂, T, u0, rngseed)
    sol = solve(prob_duffing, SKenCarp(), saveat=Ttr:Δt:T, maxiters=1e9, seed=0); 

    ts = sol.t .- sol.t[1];
    xs = sol[1, :]

    framerate = 40
    numframes = 600

    xs_well = -6.5:0.1:6.5
    a, b, c = params[1:3]

    filename = "$(plotsdir())/mechanisms/noisebistable/duffing/animation-duffing-T_$T-Ttr_$(Ttr)-framerate_$framerate-numframes_$numframes-rngseed_$rngseed.mp4"
    animate_double_well(filename, xs, ts, xs_well, x->doublewell(x, a, b, c); ball_color = BALL_COLOR, ball_size=22, c1=METASTABLE_STATE_1_COLOR, c2=METASTABLE_STATE_2_COLOR, framerate, numframes, resolution=_RESOLUTION)
end

function run_animation_chaoticsaddle()
    b = 0.9; c = 0.4; d = 6.0; a = 1.003;
    T = 110000; Ttr = 100;
    u0 = [-2.0, -2.0]
    ik = Systems.ikedamap(u0;a, b, c, d)
    traj = trajectory(ik, T; Ttr); 
    ts = Ttr:Ttr+T; xs = traj[:,1]; ys = traj[:,2]; 

    ts = ts .- ts[1];

    framerate = 25 
    numframes = 400


    fp_idx_1 = 101251 - 10 #hacky way to 
    bool_within_fp = zeros(Bool, length(traj))
    bool_within_fp[fp_idx_1:end] .= true
    c1 = METASTABLE_STATE_1_COLOR; c2 = METASTABLE_STATE_2_COLOR
    colors_state_space = [el == true ? c1 : c2 for el in bool_within_fp];
    
    #reduce density of points in saddle (to make it clearer), decrease their markersizes compared to the points going to the focus
    transition_idx = findfirst(x->x != bool_within_fp[1], bool_within_fp)
    idxs_redu = vcat( collect(1:30:transition_idx), collect(transition_idx+1:length(bool_within_fp)))
    mksize_cs = 0.8
    mksize_nd = 3.4
    markersizes = [el == true ? mksize_nd : mksize_cs for el in bool_within_fp];
    xs_fixedstructures = xs[idxs_redu]; ys_fixedstructures = ys[idxs_redu]; colors_fixedstructures = colors_state_space[idxs_redu]; markersizes_fixedstructures = markersizes[idxs_redu];
    
    idxs_animation = transition_idx-250:transition_idx+100
    ts_animation = ts[idxs_animation]; xs_animation = xs[idxs_animation]; ys_animation = ys[idxs_animation];
    colors_animation = colors_state_space[idxs_animation]
    
    filename = "$(plotsdir())/mechanisms/chaoticsaddle/ikedamap/animation-ikedamap-a_$(a)-b_$(b)-c_$(c)-d_$(d)-$framerate-$numframes-$T.mp4"
    animate_chaotic_saddle(filename, xs_animation, ys_animation, ts_animation, xs_fixedstructures, ys_fixedstructures, colors_fixedstructures, markersizes_fixedstructures, colors_animation; ball_color = BALL_COLOR, ball_size=22, c1=METASTABLE_STATE_1_COLOR, c2=METASTABLE_STATE_2_COLOR, framerate, numframes, resolution=(800, 600))
end


include("$(srcdir())/systems/lorenz.jl")
"""
Integrates the Lorenz system for parameters beefore the SN, giving a LC, and store that LC. Then integrate for after the SN, with chaos. Compares the chaotic trajecotry with the LC, and considers that the chaotic trajecotry is in the laminar period when its distance to the nearest point in the LC is below a threshold. 
"""
function run_animate_type_I_intermittency()
    σ = 10; β = 8/3; u0 = [0.1, 0.1, 0.1]
    
    #get the LC
    ρ_lc = 166.06; p_lc = [σ, ρ_lc, β]; 
    Ttr_lc = 1e3; T_lc = Ttr_lc + 3; Δt_lc = 0.0025;
    prob = ODEProblem(lorenz!, u0, (0, T_lc), p_lc); sol = solve(prob, Tsit5(), saveat=Ttr_lc:Δt_lc:T_lc, abstol=1e-8, reltol=1e-8, maxiters=1e9); t = sol.t; traj = sol[:,:]'
    limitcycle = deepcopy(traj)

    #now the chaotic trajecotry
    # Ttr_ch = 1980; T_ch = 2030; Δt_ch = 0.0025;
    Ttr_ch = 1985; T_ch = 2025; Δt_ch = 0.0025;

    ρ_ch = 166.1; p_ch = [σ, ρ_ch, β]
    prob = ODEProblem(lorenz!, u0, (0, T_ch), p_ch)
    sol = solve(prob, Tsit5(), saveat=Ttr_ch:Δt_ch:T_ch, abstol=1e-8, reltol=1e-8, maxiters=1e9); ts = sol.t .- sol.t[1]; traj = sol[:,:]'
    
    # threshold = 2; #better than 3
    threshold = 2.5; #better than 3
    arewithin = mapslices(x->withinset(x, limitcycle, threshold), traj, dims=2)[:,1]; #to distinguish between laminar and chaotic periods
    
    xs = traj[:,1]; ys = traj[:,2]; zs = traj[:,3];
    c1 = METASTABLE_STATE_1_COLOR; c2 = METASTABLE_STATE_2_COLOR
    alpha = 0.8
    colors_state_space = [el == 1 ? (c1, alpha) : (c2, alpha) for el ∈ arewithin]
    
    framerate = 25 
    numframes = 800
    
    filename = "$(plotsdir())/mechanisms/type_I_intermittency/lorenz/animation-lorenz-rho_$(ρ_ch)-$framerate-$numframes-$T_ch-$Ttr_ch.mp4"
    animate_type_I_intermittency(filename, xs, ys, zs, ts, colors_state_space, xs, ys, zs; framerate, numframes, ball_color=BALL_COLOR, ball_size=22, resolution=_RESOLUTION )
end

function run_animate_interior_crisis()
    T = 2e3; Ttr = 1e3
    framerate = 25 
    numframes = 800
    numframes = 1000
    idxs_t_anim = 600:1100
    fps_5 = [-0.149036 -0.111942; 0.6901164331310736 -0.07534279207228838; 0.7081550024831971 0.6107258394579234; 0.08793801081244812 -0.377787118392552; 1.0755007518118098 -0.25769927892208205]
    
    innercore = []
    for d in [7.22, 7.3]
        a = 0.84; b = 0.9; c = 0.4; 
        ik = Systems.ikedamap(; a, b, c, d)
        traj = trajectory(ik, T; Ttr); ts = Ttr:Ttr+T
        xs = traj[:, 1]; ys = traj[:, 2]
        filename = "$(plotsdir())/mechanisms/interior_crisis/ikeda/animation-ikedamap-d_$(d)-$framerate-$numframes-$T-$Ttr.mp4"
        
        if d <= 7.22 #
            innercore = deepcopy(Matrix(traj))
            colors_state_space = [:black for _ in ts]
        else 
            c1 = METASTABLE_STATE_1_COLOR 
            c2 = METASTABLE_STATE_2_COLOR
            # alpha = 0.6
            alpha = 1.0
            maxdist_threshold = 0.04
            # maxdist_threshold = 0.1
            arewithin = mapslices(x->withinset(x, innercore, maxdist_threshold), Matrix(traj), dims=2)[:,1]; 
            colors_state_space = [el == 1  ? (c1, alpha) : (c2, alpha) for el in arewithin]
        end
        
        xs_anim = xs[idxs_t_anim]; ys_anim = ys[idxs_t_anim]; colors_state_space_anim = colors_state_space[idxs_t_anim]; ts_anim = ts[idxs_t_anim] .- ts[1]
        animate_interior_crisis(filename, xs_anim, ys_anim, ts_anim, colors_state_space_anim, xs, ys, colors_state_space, fps_5; framerate, numframes, ball_color=BALL_COLOR, ball_size=22, resolution=(600, 600) )
    end
end

# run_animation_shc()
# run_animation_doublewell()
# run_animation_chaoticsaddle()
# run_animate_type_I_intermittency()
run_animate_interior_crisis()
    
