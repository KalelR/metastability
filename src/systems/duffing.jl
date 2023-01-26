using DynamicalSystemsBase:CDS

"""
U(x) = ax^4/4 - bx 2/2 + cx (c: assymetry paramter); gamma: dissipation, f forcing omega forcing freq, n noise strength
Alternativly, mdx2/dt2 = -νdx/dt - dV/dx + f*sen(ωt)
Default parameters
d = 1 # = ν/m = γ/m; m =1
b = 10 #-a in the paper
a = 100 # = b in the paper
c = 0
ω = 3.5 # = Ω
"""
function duffing_assymetric(u0=[0.1, 0.1]; a = 100, b=10, c=0, d=1, ω=3.5, f=0.849)
    @show a, b, c, d, f, ω
	return CDS(duffing_rule, u0, [a,b,c,d,f,ω], (y,w,z)->nothing)
end

function duffing_assymetric(pvals)
    @unpack a, b, c, d, f, ω = pvals
    df = duffing_assymetric(; a, b, c, d, f, ω)
end

@inbounds function duffing_rule(u, p, t)
    a,b,c,d,f, ω = p
    x, v = u
    du1 = v
    du2 = (-a*x^3 + b*x - c) - d*v + f*cos(ω*t)
    return SVector(du1, du2)
end

function noise_duffing(du,u,p,t)
    a,b,c,d,f, ω, n₁, n₂ = p
    du[1] = n₁
    du[2] = n₂
end

function distribution_times_noisywell(xs, numbins, side=0) #0 is left, 1 is right
    τs, _  = length_samevalues_allowfluctuations(bistable_laminarperiods(xs_sol), side)
    ws, bins = histogram(τs[1], numbins);
end

U(x, a, b, c) = (a/4)*x^4 - (b/2) * x^2 + c*x #double well


function noisy_duffing(; a=0.5, b=8.0, c=0.0, d = 0.2, f = 0.0, ω = 1.0, n₁ = n₂ = n = 0.18, T = 1e4, u0 = [0.0, 0.0])
    p =  [a, b, c, d, f, ω, n₁, n₂];
    prob_duffing = SDEProblem(duffing_rule, noise_duffing, u0, (0, T), p; seed=0)
end