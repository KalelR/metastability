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

function distribution_times_noisywell(xs, numbins, side=0) #0 is left, 1 is right
    τs, _  = length_samevalues_allowfluctuations(bistable_laminarperiods(xs_sol), side)
    ws, bins = histogram(τs[1], numbins);
end

U(x, a, b, c) = (a/4)*x^4 - (b/2) * x^2 + c*x #double well

function  integrate_and_get_distribution(; T=1e4)
    a=0.5; b=8.0; c=0.0; d = 0.2; f = 0.0; ω = 1.0; n₁ = n₂ = n = 0.18;
    Ttr = 0; Δt = 0.5
    u0 = [0.0, 0];
    prob_duffing, p = duffing(; a, b, c, d, f, ω, n₁, n₂, T, u0) 
    sol = solve(prob_duffing, SKenCarp(), saveat=Ttr:Δt:T, maxiters=1e9, progress=true);
end

function get_distribution_noisy_bistable(sol, numbins=15; obtain_dwell_time="readfromfile")
    if obtain_dwell_time == "readfromfile"
        τs_l = readdlm("$(datadir())/noisybistable-dwelltimes-left-n_0.18-T_1.0e7.dat")[:,1];
    else
        τs_l, _  = length_samevalues_allowfluctuations(bistable_laminarperiods(sol[1,:]), 0)
        τs_l = τs_l[1];
        τs_r, _  = length_samevalues_allowfluctuations(bistable_laminarperiods(sol[1,:]), 1)
        τs_r = τs_r[1];
    end

    weightsl, binsl =  histogram(τs_l, numbins);
    xfitl = binsl; yfitl = weightsl; Al, Bl = CurveFit.exp_fit(xfitl, yfitl .+ 1e-7);
    yfitl = Al .* exp.(Bl .* xfit_l)
    yfitl .+= 1e-8; #small value so log10 doesnt break
    weightsl .+= 1e-8; #small value so log10 doesnt break
    @info "Exponents of exp fit are $Al and $Bl"
    Dict("bins_l"=>binsl, "weights_l"=>weightsl, "xfit_l"=>xfitl, "yfit_l"=>yfitl, "A_l"=>Al, "B_l"=>Bl)
end