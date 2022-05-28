using GLMakie, DynamicalSystems, Random

function duffing_assymetric(u0 = [0.1, 0.25]; a = 1, b = 1, c = 0.1, d=0.1, ω = 2.2, f = 27.0,  n =0)
    J = zeros(eltype(u0), 2, 2)
    J[1,2] = 1
    # return ContinuousDynamicalSystem(duffing_assymetric_rule, u0, [ω, f, d, β, a, n], duffing_assymetric_jacob, J)
    return ContinuousDynamicalSystem(duffing_assymetric_rule, u0, [a, b, c, d, f, ω, n])
end
# """
# U(x) = ax^4/4 - bx 2/2 + cx (c: assymetry paramter); gamma: dissipation, f forcing omega forcing freq, n noise strength
# """
@inbounds function duffing_assymetric_rule(u, p, t)
    a,b,c,d,f, ω, n = p
    x, v = u
    dx1 = v
    dx2 = (-a*x^3 + b*x - c) - d*v + f*cos(ω*t) + n * (2*rand()-1)
    return SVector(dx1, dx2)
end
# @inbounds function duffing_assymetric_jacob(u, p, t)
#     ω, f, d, β, n = p
#     return @SMatrix [0 1 ;
#     (-β - 3u[1]^2) -d]
# end

a=0.5
b=8.0
c=0.0

d = 0.2

f = 0.0
ω = 1.0

n = 0.01

T = 30000
Δt=  1

Random.seed!(0)
u0 = [0.0, 0]
duf = duffing_assymetric(u0; a,b,c,d,f, ω, n)
tr = trajectory(duf, T; Δt)
t = 0:Δt:T


fig = Figure()
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[2, 1])
lines!(ax1, t, tr[:,1])

U(x, a, b, c) = (a/4)*x^4 - (b/2) * x^2 + c*x
x = -10:0.1:10
Us = U.(x, a, b, c)
lines!(ax2, x, Us)
lines!(ax2, tr[:,1], U.(tr[:,1], a, b, c))



#using SDE 
using DifferentialEquations
@inbounds function duffing_assymetric_rule(du, u, p, t)
    a,b,c,d,f, ω, n = p
    x, v = u
    du[1] = v
    du[2] = (-a*x^3 + b*x - c) - d*v + f*cos(ω*t)
end

function noise_duffing(du,u,p,t)
    a,b,c,d,f, ω, n = p
    du[1] = 0 
    du[2] = n
end
p =  [a, b, c, d, f, ω, 0.01]
tspan = (0, 10000.0)
prob_duffing = SDEProblem(duffing_assymetric_rule, noise_duffing, u0, tspan, p)
sol = solve(prob_duffing); t = sol.t

fig = Figure()
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[2, 1])
lines!(ax1, t, sol[1,:])

U(x, a, b, c) = (a/4)*x^4 - (b/2) * x^2 + c*x
x = -8:0.1:8
Us = U.(x, a, b, c)
lines!(ax2, x, Us)
lines!(ax2, sol[1,:], U.(sol[1,:], a, b, c))