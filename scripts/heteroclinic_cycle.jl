using GLMakie,  DifferentialEquations


@inbounds function heteroclinic_cycle_gh(du, u, p, t)
	μ, a, b, c = p
    du[1] = μ*u[1] + u[1]*(a*u[1]^2 + b*u[2]^2 + c*u[3]^2)
    du[2] = μ*u[2] + u[2]*(a*u[2]^2 + b*u[3]^2 + c*u[1]^2)
    du[3] = μ*u[3] + u[3]*(a*u[3]^2 + b*u[1]^2 + c*u[2]^2)
end

T = 10000
μ = 1.0 
a = 1.0 
b = 0.55
c = 1.5

μ = 1.0
a = -0.2
b = -0.1
c = -1 - (a+b)


p = [μ, a, b, c]
u0 = rand(3)
tspan = (0, T)
hcgh = ODEProblem(heteroclinic_cycle_gh, u0, tspan, p)
sol = solve(hcgh); t = sol.t


fig = Figure()
ax1 = Axis(fig[1, 1])
lines!(ax1, t, sol[1,:])
lines!(ax1, t, sol[2,:])
lines!(ax1, t, sol[3,:])

fig2 = Figure()
# scatter(sol[1,:], sol[2,:], sol[3,:])
lines(sol[1,:], sol[2,:], sol[3,:], color=1:T)
fp = sqrt(μ)
scatter!([fp, 0,0], [0, fp, 0], [0,0, fp], color=:red, markersize=5)