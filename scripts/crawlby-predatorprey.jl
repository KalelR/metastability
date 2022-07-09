using DrWatson
@quickactivate "metastability"

using GLMakie, DynamicalSystems

# """ prey = N; predator = P """
@inbounds @inline function predator_prey(u, p, t)
	α, γ, ϵ, ν, h, K, m = p
	N, P = u
    du1 = α*N*(1 - N/K) - γ*N*P / (N+h)
    du2 = ϵ * (  ν*γ*N*P/(N+h) - m*P    )
	return SVector{2}(du1, du2)
end

using DataStructures
"""
for each point in trajectory, round it to number digits (sort of equivalent of getting the eps-neighborhood of each point, eps being 10^-numdigits), count number of occurrences of rounded point and divided by total amount of points to give the measure of each point. 
"""
function histmeasure(tr, numdigits)
	v = [[tr[i,1], tr[i,2]] for i=1:length(tr)]
	vround = [round.(el, digits=numdigits) for el in v]
	c = counter(vround)
	measure = [c[elround] for elround in vround ]
	measure ./ length(v)
end




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
T = 1e5
Ttr= 100
Δt = 0.1


ϵ = 1.0
##state space 
α =1.5 
K = 10
## time-series 
# α = 0.8 
# K = 15 

u0 = rand(2)
pp = ContinuousDynamicalSystem(predator_prey, u0, [α, γ, ϵ, ν, h, K, m])
tr = trajectory(pp, T; Δt, Ttr); t = 0:Δt:T
numdigits = 3
ms = histmeasure(tr, numdigits)
tplot = Int64(300/Δt)
t = t[1:tplot]
tr = tr[1:tplot]

c2 = :blue
# fig = Figure(resolution=(4000, 2000), fontsize = 20)
fig = Figure()
ax1 = Axis(fig[1, 1], xlabel="t", ylabel="Prey", title="Crawl by in predator-prey system")
lines!(ax1, t, tr[:, 1], linewidth=1, color=:black)
ax2 = Axis(fig[1, 1], xlabel="t", ylabel="Predator", yticklabelcolor = c2, ylabelcolor=c2, yaxisposition = :right)
lines!(ax2, t, tr[:, 2], linewidth=1, color=c2)
hidespines!(ax2)
hidexdecorations!(ax2)
linkxaxes!(ax1, ax2); xlims!(0, 300)

ax3 = Axis(fig[2, 1][1,1], xlabel="Prey", ylabel="Predator")
s=scatterlines!(ax3, tr[:,1], tr[:, 2], linewidth=1, color=log.(ms), marker=:circle)
lines!(ax3, tr[:,1], tr[:, 2], linewidth=1, color=log.(ms))
Colorbar(fig[2,1][1,2], limits=(minimum(ms), maximum(ms)), label="Histogram measure, log scale")

ax4 = Axis(fig[3, 1], xlabel="measure", ylabel="PDF(measure)")
hist!(ax4, ms, bins=20, normalization=:pdf)

save("$(plotsdir())/crawlby/predatorprey-crawlby-α_$(α)-K_$(K)-numdigitsmeasure_$(numdigits)-pdf.png", fig)
# xlims!(ax1, Ttr, T)
# xlims!(ax2, Ttr, T)
# xlims!(ax3, Ttr, T)

using Symbolics
@variables N P 
@variables α γ ϵ ν h K m
J = Symbolics.jacobian([α*N*(1 - N/K) - γ*N*P / (N+h), ϵ * (  ν*γ*N*P/(N+h) - m*P    )], [N, P]; simplify=true)
J_fp = substitute.(J, (Dict(α=>1.5, γ=>2.5, ϵ=>1.0, ν=>0.5, h=>1.0, K=>10., m=>0.4, N=>0., P=>0.),))



# -------------------------------------------------------------------- fast-slow ------------------------------------------------------------------- #
T = 1e5
Ttr=500
Δt = 0.1
α = 1.5 
K = 2.2 
ϵ = 0.01


u0 = rand(2)
pp = ContinuousDynamicalSystem(predator_prey, u0, [α, γ, ϵ, ν, h, K, m])
tr = trajectory(pp, T; Δt, Ttr); t = 0:Δt:T
numdigits = 1
ms = histmeasure(tr, numdigits);

c2 = :blue
fig = Figure(resolution=(4000, 2000), fontsize = 20)
ax1 = Axis(fig[1, 1], xlabel="t", ylabel="Prey", title="Crawl by in predator-prey system")
lines!(ax1, t, tr[:, 1], linewidth=1, color=:black)
ax2 = Axis(fig[1, 1], xlabel="t", ylabel="Predator", yticklabelcolor = c2, ylabelcolor=c2, yaxisposition = :right)
lines!(ax2, t, tr[:, 2], linewidth=1, color=c2)
hidespines!(ax2)
hidexdecorations!(ax2)
linkxaxes!(ax1, ax2); xlims!(0, 1000)

ax3 = Axis(fig[2, 1][1,1], xlabel="Prey", ylabel="Predator")
s=scatterlines!(ax3, tr[:,1], tr[:, 2], linewidth=1, color=log.(ms), marker=:circle)
Colorbar(fig[2,1][1,2], limits=(minimum(ms), maximum(ms)), label="Histogram measure, log scale")
# save("$(plotsdir())/fastslow/predatorprey-fastslow-α_$(α)-K_$(K)-ϵ_$(ϵ)-numdigitsmeasure_$(numdigits).png", fig)

ax4 = Axis(fig[3, 1], xlabel="measure", ylabel="PDF(measure)")
hist!(ax4, ms, bins=20, normalization=:pdf)
