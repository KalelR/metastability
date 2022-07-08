using GLMakie,  DifferentialEquations


gh_coupling(γ, x, y, z) = γ*x*(y^2+z^2)

@inbounds function heteroclinic_chaos_lorenz(du, u, p, t)
	# x, y, z = u
	x = u[1:3]
	y = u[4:6]
	z = u[7:9]
	μ, σ, p, β = p
	p_lor = [σ, ρ, β]
    du[1] = Systems.loop(x, p_lor, t) + gh_coupling(μ, x, y, z)
    du[1] = Systems.loop(y, p_lor, t) + gh_coupling(μ, y, x, z)
    du[1] = Systems.loop(z, p_lor, t) + gh_coupling(μ, z, x, y)
end

T = 5000
# T = 200000
μ = -2 
σ = 15
ρ = 58
β = 2.4
u01 = [10, -11, 30]
u02 = [10, -13, 20]
u03 = [10, -12, 30]
u0 = [u01; u02; u03]
p = [μ, σ, ρ, β]

hcgh = ODEProblem(heteroclinic_chaos_lorenz, u0, tspan, p)
sol = solve(hcgh, Rodas5(), saveat=0:0.01:T, abstol=1e-10, reltol=1e-10, maxiters=1e9); t = sol.t;


fig = Figure()
ax1 = Axis(fig[1, 1])
lines!(ax1, t, sol[1,:])
lines!(ax1, t, sol[4,:])
lines!(ax1, t, sol[7,:])




# fig2 = Figure()
ax2 = Axis3(fig[1,2])
# scatter(sol[1,:], sol[2,:], sol[3,:])
lines!(ax2, sol[1,:], sol[2,:], sol[3,:], color=t)
scatter!(ax2, fps[:,1], fps[:,2], fps[:,3], color=:red, markersize=5000)
scatter!(ax2, [0], [0], [0], color=:blue, markersize=8000)