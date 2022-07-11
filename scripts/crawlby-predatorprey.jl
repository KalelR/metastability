# using DrWatson
# @quickactivate "metastability"

# using GLMakie, DynamicalSystems
using GLMakie, OrdinaryDiffEq
plotsdir() = "../plots/"

# """ prey = N; predator = P """
@inbounds @inline function predator_prey!(du, u, p, t)
	α, γ, ϵ, ν, h, K, m = p
	N, P = u
    du[1] = α*N*(1 - N/K) - γ*N*P / (N+h)
    du[2] = ϵ * (  ν*γ*N*P/(N+h) - m*P    )
	# return SVector{2}(du1, du2)
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


timederivative(sol) = [sol(t, Val{1}) for t ∈ sol.t]
norm(v) = sum(v.^2)
γ = 2.5
h = 1 
ν = 0.5 
m = 0.4



#crawlby (fast-slow dynamics) (due to saddle point)
typename = "crawlby"
ϵ = 1.0
T = 1e5
Ttr= 100
Δt = 0.01
Tplot = 300
##state space  fig
α =1.5 
K = 10
## time-series figs
# α = 0.8 
# K = 15 
#

#fast-slow system 
typename = "fastslowsystem"
T = 1e5
Ttr=500
Δt = 0.01
α = 1.5 
K = 2.2 
ϵ = 0.01
Tplot = 1000
#

#regular LC
# typename = "regularLC"
# ϵ = 1.0
# T = 1e5
# Ttr= 100
# Δt = 0.01
# α = 1.5 
# K = 2.2 
# Tplot=300
#

u0 = [0.1, 0.1]
p = [α, γ, ϵ, ν, h, K, m]
# pp = ContinuousDynamicalSystem(predator_prey, u0, p)
# tr = trajectory(pp, T; Δt, Ttr); t = 0:Δt:T;
solver = Tsit5()
tspan = (0, T)
prob = ODEProblem(predator_prey!, u0, tspan, p)
sol = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8)
t = sol.t 
tr = sol[:,:]'

speeds = norm.(timederivative(sol))

numdigits = 5
ms = histmeasure(sol', numdigits);
tplot = Int64(Tplot/Δt);
t_plot = t[1:tplot];
tr_plot = tr[1:tplot, :];
ms_plot = ms[1:tplot];
speeds_plot = speeds[1:tplot]

c2 = :blue
# fig = Figure(resolution=(4000, 2000), fontsize = 20)
fig = Figure()
ax1 = Axis(fig[1, 1], xlabel="t", ylabel="Prey", title="$(typename) in predator-prey system")
lines!(ax1, t_plot, tr_plot[:, 1], linewidth=1, color=:black)
ax2 = Axis(fig[1, 1], xlabel="t", ylabel="Predator", yticklabelcolor = c2, ylabelcolor=c2, yaxisposition = :right)
lines!(ax2, t_plot, tr_plot[:, 2], linewidth=1, color=c2)
hidespines!(ax2)
hidexdecorations!(ax2)
linkxaxes!(ax1, ax2); xlims!(Ttr, Tplot)


ax3 = Axis(fig[2, 1][1,1], xlabel="Prey", ylabel="Predator")
# s=scatter!(ax3, tr_plot[:,1], tr_plot[:, 2], linewidth=1, color=log.(ms_plot), marker=:circle)
s=scatter!(ax3, tr_plot[:,1], tr_plot[:, 2], linewidth=1, color=ms_plot, marker=:circle)
Colorbar(fig[2,1][1,2], limits=(minimum(ms_plot), maximum(ms_plot)), label="Measure, log scale")

ax4 = Axis(fig[3, 1][1,1], xlabel="Prey", ylabel="Predator")
# s=scatterlines!(ax4, tr_plot[:,1], tr_plot[:, 2], linewidth=1, color=speeds_plot, marker=:circle)
s=scatter!(ax4, tr_plot[:,1], tr_plot[:, 2], linewidth=1, color=log.(speeds_plot), marker=:circle)
Colorbar(fig[3,1][1,2], limits=(minimum(speeds_plot), maximum(speeds_plot)), label="|xdot|, log scale")

ax5 = Axis(fig[4, 1], xlabel="measure", ylabel="PDF(measure)")
# hist!(ax4, ms, bins=20, normalization=:pdf)
hist!(ax5, ms, bins=20, normalization=:pdf)

save("$(plotsdir())/$(typename)/predatorprey-$(typename)-α_$(α)-K_$(K)-numdigitsmeasure_$(numdigits)-intstep_$(Δt)-pdf.png", fig)
# xlims!(ax1, Ttr, T)
# xlims!(ax2, Ttr, T)
# xlims!(ax3, Ttr, T)

using Symbolics
@variables N P 
@variables α γ ϵ ν h K m
J = Symbolics.jacobian([α*N*(1 - N/K) - γ*N*P / (N+h), ϵ * (  ν*γ*N*P/(N+h) - m*P    )], [N, P]; simplify=true)
J_fp = substitute.(J, (Dict(α=>1.5, γ=>2.5, ϵ=>1.0, ν=>0.5, h=>1.0, K=>10., m=>0.4, N=>0., P=>0.),))



