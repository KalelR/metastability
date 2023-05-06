using GLMakie,  DifferentialEquations, StaticArrays

# function lorenz_stat(u, p, t)
#     @inbounds begin
#         σ = p[1]; ρ = p[2]; β = p[3]
#         du1 = σ*(u[2]-u[1])
#         du2 = u[1]*(ρ-u[3]) - u[2]
#         du3 = u[1]*u[2] - β*u[3]
#         return SVector{3}(du1, du2, du3)
#     end
# end
# function lorenz_op(u, p, t)
#     @inbounds begin
#         σ = p[1]; ρ = p[2]; β = p[3]
#         du1 = σ*(u[2]-u[1])
#         du2 = u[1]*(ρ-u[3]) - u[2]
#         du3 = u[1]*u[2] - β*u[3]
#         return [du1, du2, du3]
#     end
# end
# @time a=lorenz_op(u, p, t)

# gh_coupling(γ, x, y, z) = γ*x*(y^2+z^2)
# normsq(v) = sum(v[1]^2 + v[2]^2 + v[3]^2)
# @inbounds function heteroclinic_chaos_lorenz!(du, u, p, t)
# 	x = @view u[1:3]
# 	y = @view u[4:6]
# 	z = @view u[7:9]
# 	p_lor = @view p[2:end]
# 	μ = p[1]
#     # du[1:3] .= lorenz_op(x, p_lor, t) .+ μ .* x .* normsq(y)  #3 alloc
#     # du[1:3] = lorenz_op(x, p_lor, t) .+ μ .* x .* normsq(y)  #3 alloc
#     # du[1:3] .= lorenz_op(x, p_lor, t) #2 aloc
#     # du[1:3] .=  μ .* x .* normsq(y)  #3 alloc
#     du[1:3] .=  [1,2,3]  #2 alloc
#     # du[1:3] =  [1,2,3]  #1 alloc
#     # du[1] =  1  #3 alloc
#     # du[4:6] .= lorenz_op(y, p_lor, t) .+ μ .* y .* normsq(z) 
#     # du[7:9] .= lorenz_op(z, p_lor, t) .+ μ .* z .* normsq(x) 
# end



# @inbounds function heteroclinic_chaos_lorenz(u, p, t)
# 	x = @view u[1:3]
# 	y = @view u[4:6]
# 	z = @view u[7:9]
# 	p_lor = @view p[2:end]
# 	μ = p[1]
#     du1 = lorenz_op(x, p_lor, t) .+ μ .* x .* normsq(y) 
#     du2 = lorenz_op(y, p_lor, t) .+ μ .* y .* normsq(z) 
#     du3 = lorenz_op(z, p_lor, t) .+ μ .* z .* normsq(x) 
# 	[du1; du2; du3]
# end

# @inbounds function heteroclinic_chaos_lorenz_stat(u, p, t)
# 	x = @view u[1:3]
# 	y = @view u[4:6]
# 	z = @view u[7:9]
# 	p_lor = @view p[2:end]
# 	μ = p[1]
#     du1 = lorenz_stat(x, p_lor, t) .+ μ .* x .* normsq(y) 
#     du2 = lorenz_stat(y, p_lor, t) .+ μ .* y .* normsq(z) 
#     du3 = lorenz_stat(z, p_lor, t) .+ μ .* z .* normsq(x) 
# 	SA[du1; du2; du3]
# end

# du = zeros(9)
# u = zeros(9)
# p = [1,1,1,1]
# t = 0.0
# @time heteroclinic_chaos_lorenz!(du, u, p, t)
# @time heteroclinic_chaos_lorenz_loop!(du, u, p, t)
# @time heteroclinic_chaos_lorenz(u, p, t)
# @time a=heteroclinic_chaos_lorenz_stat(u, p, t)
# @code_warntype heteroclinic_chaos_lorenz!(du, u, p, t)

"""
Optimized code.. maybe could look better
"""
@inbounds function heteroclinic_chaos_lorenz!(du, u, p, t)
	μ, σ, ρ, β = p
	for i=1:3
		X = @view u[(i-1)*3+1:i*3]
		if i == 3 Y = @view u[1:3] else Y = @view u[i*3+1:(i+1)*3] end
		normvec = normsq(Y);
		du[(i-1)*3 + 1] = σ*(X[2] - X[1])      + μ * X[1] * normvec
		du[(i-1)*3 + 2] = X[1]*(ρ-X[3]) - X[2] + μ * X[2] * normvec
		du[(i-1)*3 + 3] = X[1]*X[2] - β*X[3]   + μ * X[3] * normvec
	end
end

# T = 5000
T = 25000
tspan = (0, T)
# T = 200000
μ = -0.025
# μ = 0.0
σ = 15
ρ = 58
β = 2.4

# σ = 10
# ρ = 28 
# β = 8/3

u01 = [10, -11, 30]
u02 = [10, -13, 20]
u03 = [10, -12, 30]
u0 = [u01; u02; u03]
p = [μ, σ, ρ, β]

hcgh = ODEProblem(heteroclinic_chaos_lorenz!, u0, tspan, p)
# sol = solve(hcgh, Rodas5(), saveat=0:0.01:T, abstol=1e-10, reltol=1e-10, maxiters=1e9); t = sol.t;
# sol = solve(hcgh, Rodas5(), saveat=0:0.01:T, abstol=1e-8, reltol=1e-8, maxiters=1e8); solver_s = "rodas5"; t = sol.t;
# sol = solve(hcgh, RK4(), dt=1/10000); t = sol.t;
sol = solve(hcgh, RK4(), dt=1/10); solver_s = "RK4"; t = sol.t;

c1 = (:red , 1.0)
c2 = (:purple , 1.0)
c3 = (:green, 1.0)
colors = [c1, c2, c3]
axs = []
fig = Figure()
for i=1:3 #each cell
	x = sol[(i-1)*3+1,:];
	y = sol[(i-1)*3+2,:];
	z = sol[(i-1)*3+3,:];
	ax = Axis(fig[i, 1], ylabel="Cell $(i) \n x"); push!(axs, ax);
	lines!(ax, t, x, color=colors[i])
	ax = Axis(fig[i, 2], ylabel="Cell $(i) \n y"); push!(axs, ax);
	lines!(ax, t, y, color=colors[i])
	ax = Axis(fig[i, 3], ylabel="Cell $(i) \n z", xlabel="t"); push!(axs, ax);
	lines!(ax, t, z, color=colors[i])
	ax = Axis3(fig[i,4])
	lines!(ax, x, y, z, color=colors[i]); 
end
linkxaxes!(axs...)

save("../plots/heteroclinicchaos/heteroclinichaos-lorenz-sigma_$(σ)-rho_$(ρ)-beta-$(β)-mu_$(μ)-solver_$(solver_s).png", fig)









# fig2 = Figure()
ax2 = Axis3(fig[1,2])
# scatter(sol[1,:], sol[2,:], sol[3,:])
lines!(ax2, sol[1,:], sol[2,:], sol[3,:], color=t)
scatter!(ax2, fps[:,1], fps[:,2], fps[:,3], color=:red, markersize=5000)
scatter!(ax2, [0], [0], [0], color=:blue, markersize=8000)

fig = Figure()