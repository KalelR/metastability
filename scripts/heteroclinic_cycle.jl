using GLMakie,  DifferentialEquations


@inbounds function heteroclinic_cycle_gh(du, u, p, t)
	μ, a, b, c = p
    du[1] = μ*u[1] + u[1]*(a*u[1]^2 + b*u[2]^2 + c*u[3]^2)
    du[2] = μ*u[2] + u[2]*(a*u[2]^2 + b*u[3]^2 + c*u[1]^2)
    du[3] = μ*u[3] + u[3]*(a*u[3]^2 + b*u[1]^2 + c*u[2]^2)
end

function gh_jac(u, p)
    x,y,z = u
    l, a, b, c = p
    return [(l + 3*a*x^2 + b*y^2 + c*z^2) (2*b*x*y) (2*c*x*z); (2*c*y*x) (l+3*a*y^2 + b*z^2 + c*x^2) (2*b*y*z); (2*b*z*x) (2*c*z*y) (l+3*a*z^2 + b*x^2 +c*y^2)]
end



T = 5000
# T = 200000
μ = 1.0 

##GH - stable HC!!!!!
a = -0.2
b = -0.1
# a = -0.3 
# b = -0.2
# a = -0.5
# b = -0.01
c = -1 - a - b

## stabilitizes and becomes node apparently
# T = 1500
# a = -1
# b = -1.5
# c = (2a - b) - 0.1 #unstable ## 
# c = (2a - b) + 0.1#not even a HC 


## parameters from scholarpedia (adapted so c>a>b); unstable
# c = -0.55 
# a = -1.0
# b = -1.5

## directly from scholarpedia; unstable
# a = -1.0 
# b = -0.55 
# c = -1.5

@warn c > a > b #Krupa exists if true 
@warn 2a > b + c #Krupa stable if true
@warn a + b + c == -1 #GH
@warn -1/3 < a < 0 #GH 
@warn c < a < b < 0 #GH


p = [μ, a, b, c]


xpcoor(μ,p) = sqrt(-μ/p)
fps_x = [[xpcoor(μ, a), 0, 0], [xpcoor(μ, b), 0, 0], [xpcoor(μ, c), 0, 0]]; fps_x = [fps_x; -1 .* fps_x]
fps_y = [[0, xpcoor(μ, a), 0], [0, xpcoor(μ, b), 0], [0, xpcoor(μ, c), 0]]; fps_y = [fps_y; -1 .* fps_y]
fps_z = [[0, 0, xpcoor(μ, a)], [0, 0, xpcoor(μ, b)], [0, 0, xpcoor(μ, c)]]; fps_z = [fps_z; -1 .* fps_z]
fps = hcat([fps_x; fps_y; fps_z]...)'

# u0 = fps_y[3] .+ [0.1, -0.5, 0.1] #still converges to node, this fp is repeller
# u0 = fps_y[1] .+ [0.1, +0.5, 0.1]
u0 = [0.1, 0.2, 0.3] #off the coordinate plane or axis should be attracted towards HC
# u0 = 5 .* rand(3)
# u0 = [5.3, 5.2, 5.4]
# u0 = [-0.1, -0.2, -0.3] 
# u0 = [10,10,10] #goes to a node
# u0 = [10,-10,10] #goes to a node

tspan = (0, T)
hcgh = ODEProblem(heteroclinic_cycle_gh, u0, tspan, p)
# sol = solve(hcgh, saveat=0:0.01:T, abstol=1e-16, reltol=1e-16); t = sol.t;
# sol = solve(hcgh, AutoTsit5(Rosenbrock23()), saveat=0:0.01:T, abstol=1e-8, reltol=1e-8); t = sol.t;
sol = solve(hcgh, Rodas5(), saveat=0:0.01:T, abstol=1e-10, reltol=1e-10, maxiters=1e9); t = sol.t;
# sol = solve(hcgh, RK4(), dt=1/1000); t = sol.t;
# sol = solve(hcgh, saveat=0:0.01:T); t = sol.t;


fig = Figure()
ax1 = Axis(fig[1, 1])
lines!(ax1, t, sol[1,:])
lines!(ax1, t, sol[2,:])
lines!(ax1, t, sol[3,:])

# fig2 = Figure()
ax2 = Axis3(fig[1,2])
# scatter(sol[1,:], sol[2,:], sol[3,:])
lines!(ax2, sol[1,:], sol[2,:], sol[3,:], color=t)
scatter!(ax2, fps[:,1], fps[:,2], fps[:,3], color=:red, markersize=5000)
scatter!(ax2, [0], [0], [0], color=:blue, markersize=8000)




## Stability analysis
using LinearAlgebra
xp = sol[end]
xp = fps_y[2]
# xp = fps_y[2]
# xp = fps_y[3]
# xp = [sqrt(-1/a), 0, 0]
J = gh_jac(xp, p)
evals, evecs = eigen(J)


# With noise
@inbounds function heteroclinic_cycle_gh_noise(du, u, p, t)
	μ, a, b, c, n = p
    du[1] = μ*u[1] + u[1]*(a*u[1]^2 + b*u[2]^2 + c*u[3]^2)
    du[2] = μ*u[2] + u[2]*(a*u[2]^2 + b*u[3]^2 + c*u[1]^2)
    du[3] = μ*u[3] + u[3]*(a*u[3]^2 + b*u[1]^2 + c*u[2]^2)
end
function noise_hc_gh(du,u,p,t)
	μ, a, b, c, n = p
    du[1] = abs(n) 
    du[2] = abs(n)
    du[3] = abs(n)
end
Δt = 0.01
n = 1e-8
p = [μ, a, b, c, n]
prob_hc = SDEProblem(heteroclinic_cycle_gh_noise, noise_hc_gh, u0, tspan, p; seed=0)
sol = solve(prob_hc, SOSRA(), saveat=0:Δt:T); t = sol.t #nonstiff, didnt test
# sol = solve(prob_hc, SKenCarp(), saveat=0:Δt:T); t = sol.t #stiff, worked well but slow

fig = Figure()
ax1 = Axis(fig[1, 1])
lines!(ax1, t, sol[1,:])
lines!(ax1, t, sol[2,:])
lines!(ax1, t, sol[3,:])

# fig2 = Figure()
ax2 = Axis3(fig[1,2])
# scatter(sol[1,:], sol[2,:], sol[3,:])
lines!(ax2, sol[1,:], sol[2,:], sol[3,:], color=t)

scatter!(ax2, fps[:,1], fps[:,2], fps[:,3], color=:red, markersize=5000)
scatter!(ax2, [0], [0], [0], color=:blue, markersize=8000)