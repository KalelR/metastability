using CairoMakie,  DifferentialEquations


# """
# U(x) = ax^4/4 - bx 2/2 + cx (c: assymetry paramter); gamma: dissipation, f forcing omega forcing freq, n noise strength
# """
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

a=0.5
b=8.0
c=0.0

d = 0.2

f = 0.0
ω = 1.0

n = 2.5

T = 30000.

u0 = [0.0, 0]

p =  [a, b, c, d, f, ω, n]
tspan = (0, T)
prob_duffing = SDEProblem(duffing_assymetric_rule, noise_duffing, u0, tspan, p; seed=0)
sol = solve(prob_duffing); t = sol.t

#animation of double well
framerate = 20
datastep = 2
dataduration = datastep *1000
dataduration_2 = datastep *1000
videotstep= 1


energy(x, a, b, c) = (a/4)*x^4 - (b/2) * x^2 + c*x
timestamps = range(1, length(t), step=1)
time = Observable(1)
x = @lift(sol[1, $time])
E = @lift(energy(sol[1, $time], a,b,c))
fig = scatter([x], [E], color = :black)
filename = "$animation-double-well-datastep_$(datastep)-framerate_$(framerate)-videotend_$(dataduration).mp4"; mkpath(dirname(filename))
record(fig, filename, timestamps; framerate = framerate) do t
    time[] = t
end
