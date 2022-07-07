using GLMakie,  DifferentialEquations


@inbounds function heteroclinic_cycle_gh(du, u, p, t)
	μ, a, b, c = p
    du[1] = μ*u[1] + u[1]*(a*u[1]^2 + b*u[2]^2 + c*u[3]^2)
    du[2] = μ*u[2] + u[2]*(a*u[2]^2 + b*u[3]^2 + c*u[1]^2)
    du[3] = μ*u[3] + u[3]*(a*u[3]^2 + b*u[1]^2 + c*u[2]^2)
end

T = 5000
μ = 1.0 

##GH - stable HC!!!!!
a = -0.2
b = -0.1
a = -0.3 
b = -0.2
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
u0 = [0.1, 0.2, 0.3] #off the coordinate plane or axis should be attracted towards HC

u0 = 2 .* rand(3)
# u0 = [10,10,10] #goes to a node
# u0 = [10,-10,10] #goes to a node

tspan = (0, T)
hcgh = ODEProblem(heteroclinic_cycle_gh, u0, tspan, p)
sol = solve(hcgh, saveat=0:0.01:T, abstol=1e-16, reltol=1e-16); t = sol.t;
# sol = solve(hcgh, AutoTsit5(Rosenbrock23()), saveat=0:0.01:T, abstol=1e-8, reltol=1e-8); t = sol.t;
# sol = solve(hcgh, saveat=0:T); t = sol.t;


fig = Figure()
ax1 = Axis(fig[1, 1])
lines!(ax1, t, sol[1,:])
lines!(ax1, t, sol[2,:])
lines!(ax1, t, sol[3,:])

# fig2 = Figure()
ax2 = Axis3(fig[1,2])
# scatter(sol[1,:], sol[2,:], sol[3,:])
lines!(ax2, sol[1,:], sol[2,:], sol[3,:], color=t)
scatter!(ax2, [sqrt(-μ/a), 0, 0], [0, sqrt(-μ/b), 0], [0, 0, sqrt(-μ/c)], color=:red, markersize=1000)