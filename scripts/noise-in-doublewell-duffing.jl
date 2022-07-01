using GLMakie,  DifferentialEquations


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

T = 1000
Δt = 0.5

u0 = [0.0, 0]

p =  [a, b, c, d, f, ω, n]
tspan = (0, T)
prob_duffing = SDEProblem(duffing_assymetric_rule, noise_duffing, u0, tspan, p; seed=0)
# sol = solve(prob_duffing, SOSRA(), saveat=0:Δt:T); t = sol.t #nonstiff, didnt test
sol = solve(prob_duffing, SKenCarp(), saveat=0:Δt:T); t = sol.t #stiff, worked well but slow

fig = Figure()
ax1 = Axis(fig[1, 1], ylabel="x", xlabel="t")
ax2 = Axis(fig[2, 1], ylabel="U(x)", xlabel="x")
lines!(ax1, t, sol[1,:])

U(x, a, b, c) = (a/4)*x^4 - (b/2) * x^2 + c*x
x = -6.5:0.1:6.5
Us = U.(x, a, b, c)
lines!(ax2, x, Us, color=:black, linewidth=4)
scatter!(ax2, sol[1,:], U.(sol[1,:], a, b, c), color=(:orange, 0.5))

save("duffing-doublwell-a_$(a)-b_$(b)-c_$(c)-d_$(d)-f_$(f)-n_$(n)-solver_SKenCarp.png", fig)


