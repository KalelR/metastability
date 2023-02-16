using DynamicalSystemsBase:CDS

@inbounds function duffing_rule!(du, u, p, t)
    a,b,c,d,f, ω, n₁, n₂ = p
    x, v = u
    du[1] = v
    du[2] = (-a*x^3 + b*x - c) - d*v + f*cos(ω*t)
end

@inbounds function duffing_noise_rule!(du,u,p,t)
    a,b,c,d,f, ω, n₁, n₂ = p
    du[1] = n₁
    du[2] = n₂
end

"""
This function establishes the default parameters used in the simulations and then creates the structure in the code that is needed to numerically integrate.
U(x) = ax^4/4 - bx 2/2 + cx (c: assymetry paramter); gamma: dissipation, f forcing omega forcing freq, n noise strength
Alternativly, mdx2/dt2 = -νdx/dt - dV/dx + f*sen(ωt)
Default parameters
d = 1 # = ν/m = γ/m; m =1
b = 10 #-a in the paper
a = 100 # = b in the paper
c = 0
ω = 3.5 # = Ω
"""
function duffing(; a=0.5, b=8.0, c=0.0, d = 0.2, f = 0.0, ω = 1.0, n₁ = 0.18, n₂ = n₁, T = 1e4, u0 = [0.0, 0.0])
    p =  [a, b, c, d, f, ω, n₁, n₂];
    if n₁ != 0 || n₂ != 0
        prob_duffing = SDEProblem(duffing_rule!, duffing_noise_rule!, u0, (0, T), p; seed=0)
    else
        prob_duffing = ODEProblem(duffing_rule!, u0, (0, T), p; seed=0)
    end
    return prob_duffing, p
end

function distribution_times_up_below_threshold(d::Dict, numbins, side=0; Δt=1.0) #0 is left, 1 is right
    v = d["v"]
    return distribution_times_up_below_threshold(v, numbins, side)
end

function distribution_times_up_below_threshold(v, numbins, side=0; Δt=1.0) #0 is left, 1 is right
    dwelltimes_int, _  = length_samevalues_allowfluctuations(bistable_laminarperiods(v), side)
    dwelltimes = multiply_dict_number(dwelltimes_int, Δt)
    ws, bins = histogram(dwelltimes[1], numbins);
end


U(x, a, b, c) = (a/4)*x^4 - (b/2) * x^2 + c*x #double well

# function integrate_and_get_distribution(; T=1e4)
#     a=0.5; b=8.0; c=0.0; d = 0.2; f = 0.0; ω = 1.0; n₁ = n₂ = n = 0.18;
#     Ttr = 0; Δt = 0.5
#     u0 = [0.0, 0];
#     prob_duffing, p = duffing(; a, b, c, d, f, ω, n₁, n₂, T, u0) 
#     sol = solve(prob_duffing, SKenCarp(), saveat=Ttr:Δt:T, maxiters=1e9, progress=true);
# end

function integrate_noisy_bistable(d::Dict)
    @unpack T, Ttr, Δt = d
    prob, params = duffing(; T)
    sol = solve(prob, SKenCarp(), saveat=Ttr:Δt:T, maxiters=1e9, progress=true, seed=2);
    return @strdict sol
end 


function integrate_and_get_distribution_info_noisy_bistable(d::Dict)
    @unpack T, Ttr, Δt, generate_data, numbins = d
    filename = "$(datadir())/noisybistable-solution-T_$T-Ttr_$Ttr-Δt_$Δt.jld2"
    sol_data, file = produce_or_load(d, integrate_noisy_bistable; filename, force=generate_data, tag=false)
    @unpack sol = sol_data; 
    distribution_data = distribution_info_noisy_bistable(sol[1,:], numbins; Δt)
end
    

function distribution_info_noisy_bistable(d::Dict)
    @unpack xs, numbins, Δt, numbins = d
    return distribution_info_noisy_bistable(xs, numbins; Δt)
end

function distribution_info_noisy_bistable(xs, numbins=10; Δt=1.0)
    distribution_data = Dict()
    for side in [0, 1]
        weights, bins = distribution_times_up_below_threshold(xs, numbins, side; Δt)
        if minimum(weights) == 0 weights .+= 1e-8 end #just make sure that log and fitting don't error because of weights=0
        xfit = bins; yfit = deepcopy(weights) 
        A, B = CurveFit.exp_fit(xfit, yfit);
        yfit = A .* exp.(B .* xfit)
        data = @strdict weights bins xfit yfit A B
        distribution_data["$side"] = data
    end
    return distribution_data
end


function dwell_times_attractor_merging_crisis(d::Dict)
    @unpack xs, numbins, Δt = d
    return dwell_times_attractor_merging_crisis(xs, numbins; Δt)
end

function dwell_times_attractor_merging_crisis(xs, numbins=10; Δt=1.0)
    idx_states_time = [x >= 0 ? 1 : 0 for x in xs]
    dwelltimes = Dict()
    for (i, key) in enumerate([0, 1])
        dwelltimes["$key"] = Δt .* length_samevalues_allowfluctuations(idx_states_time, 3)[1][key]
    end
    return dwelltimes
end