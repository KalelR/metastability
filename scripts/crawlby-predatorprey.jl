using DrWatson
@quickactivate "metastability"

using GLMakie, DynamicalSystems

""" prey = N; predator = P """
@inbounds @inline function predator_prey(u, p, t)
	α, γ, ϵ, ν, h, K, m = p
	N, P = u
    du1 = α*N*(1 - N/K) - γ*N*P / (N+h)
    du2 = ϵ * (  ν*γ*N*P/(N+h) - m*P    )
	return SVector{2}(du1, du2)
end

T = 300
Ttr=100
Δt = 0.1

γ = 2.5
h = 1 
ν = 0.5 
m = 0.4

#normal LC  
α = 1.5 
K = 2.2 
ϵ = 1.0

##boring, wont do

# -------------------------------------------------------------------- crawl-by -------------------------------------------------------------------- #
ϵ = 1.0
##state space 
# α =1.5 
# K = 10
## time-series 
α = 0.8 
K = 15 

u0 = rand(2)
pp = ContinuousDynamicalSystem(predator_prey, u0, [α, γ, ϵ, ν, h, K, m])
tr = trajectory(pp, T; Δt, Ttr); t = 0:Δt:T


fig = Figure(resolution=(4000, 2000), fontsize = 20)
ax1 = Axis(fig[1, 1], xlabel="t", ylabel="Prey")
ax2 = Axis(fig[1, 1], xlabel="t", ylabel="Predator")
ax3 = Axis(fig[2, :], xlabel="Prey", ylabel="Predator")
scatterlines!(ax3, tr[:,1], tr[:, 2], linewidth=1, color=t, marker=:circle)
lines!(ax1, t, tr[:, 1], linewidth=1, color=:black)
lines!(ax2, t, tr[:, 2], linewidth=1, color=:blue)
save("$(plotsdir())/crawlby/predatorprey-crawlby-α_$(α)-K_$(K).png", fig)
# xlims!(ax1, Ttr, T)
# xlims!(ax2, Ttr, T)
# xlims!(ax3, Ttr, T)

using Symbolics
@variables N P 
@variables α γ ϵ ν h K m
J = Symbolics.jacobian([α*N*(1 - N/K) - γ*N*P / (N+h), ϵ * (  ν*γ*N*P/(N+h) - m*P    )], [N, P]; simplify=true)
J_fp = substitute.(J, (Dict(α=>1.5, γ=>2.5, ϵ=>1.0, ν=>0.5, h=>1.0, K=>10., m=>0.4, N=>0., P=>0.),))


# -------------------------------------------------------------------- fast-slow ------------------------------------------------------------------- #
T = 1000
Ttr=500
Δt = 0.1
α = 1.5 
K = 2.2 
ϵ = 0.01


u0 = rand(2)
pp = ContinuousDynamicalSystem(predator_prey, u0, [α, γ, ϵ, ν, h, K, m])
tr = trajectory(pp, T; Δt, Ttr); t = 0:Δt:T


fig = Figure(resolution=(4000, 2000), fontsize = 20)
ax1 = Axis(fig[1, 1], xlabel="t", ylabel="Prey-Predator")
# ax2 = Axis(fig[1, 1], xlabel="t", ylabel="Predator")
ax3 = Axis(fig[2, :], xlabel="Prey", ylabel="Predator")
scatterlines!(ax3, tr[:,1], tr[:, 2], linewidth=1, color=t, marker=:circle)
lines!(ax1, t, tr[:, 1], linewidth=1, color=:black, label="Prey")
lines!(ax1, t, tr[:, 2], linewidth=1, color=:blue, label="Predator")
fig[1,2] = Legend(fig, ax1)
save("$(plotsdir())/fastslow/predatorprey-fastslow-α_$(α)-K_$(K)-ϵ_$(ϵ).png", fig)
# xlims!(ax1, Ttr, T)
# xlims!(ax2, Ttr, T)
# xlims!(ax3, Ttr, T)