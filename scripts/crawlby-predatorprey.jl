using DrWatson
@quickactivate "metastability"
using GLMakie, OrdinaryDiffEq
include("$(scriptsdir())/utils.jl")

# """ prey = N; predator = P """
@inbounds @inline function predator_prey!(du, u, p, t)
	α, γ, ϵ, ν, h, K, m = p
	N, P = u
    du[1] = α*N*(1 - N/K) - γ*N*P / (N+h)
    du[2] = ϵ * (  ν*γ*N*P/(N+h) - m*P    )
	# return SVector{2}(du1, du2)
end


γ = 2.5
h = 1 
ν = 0.5 
m = 0.4


#crawlby (fast-slow dynamics) (due to saddle point)
typename = "crawlby"
ϵ = 1.0
# T = 1e5
T = 2e3
Ttr= 100
Ttr=1e4; T = 15000 #for α=0.1
Δt = 0.01
# Tplot = 300
Tplot=1000
#state space  fig
α = 1.5 
α = 1.0
α = 0.5
α = 0.3
α = 0.1
K = 10
numbins = 80
# time-series figs
# α = 0.8 
# K = 15 
dirname=typename
#

#fast-slow system 
# typename = "fastslowsystem"
# T = 1e5
# Ttr=500
# Δt = 0.01
# α = 1.5 
# K = 2.2 
# ϵ = 0.01
# Tplot = 2000
# numbins = 40
# dirname=typename
# #

# #regular LC
# typename = "regularLC"
# dirname = "not_metastable/$(typename)"
# ϵ = 1.0
# # T = 1e5
# T = 300
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
fig = Figure(resolution=(800, 800), fontsize = 20)
# ax1 = Axis(fig[1, 1:2], xlabel="t", ylabel="Prey", title="$(typename) in predator-prey system, min = ($(round(minimum(tr[:,1]), digits=6)), $(round(minimum(tr[:,2]), digits=6)))")
ax1 = Axis(fig[1, 1:2], xlabel="t", ylabel="Prey", title="$(typename) in predator-prey system")
lines!(ax1, t_plot, tr_plot[:, 1], linewidth=1, color=:black)
ax2 = Axis(fig[1, 1:2], xlabel="t", ylabel="Predator", yticklabelcolor = c2, ylabelcolor=c2, yaxisposition = :right)
lines!(ax2, t_plot, tr_plot[:, 2], linewidth=1, color=c2)
hidespines!(ax2)
hidexdecorations!(ax2)
linkxaxes!(ax1, ax2); 
xlims!(ax1, Ttr, Ttr+Tplot)

ax3 = Axis(fig[2, 1][1,1], xlabel="Prey", ylabel="Predator")
s=scatter!(ax3, tr_plot[:,1], tr_plot[:, 2], linewidth=1, color=log10.(ms_plot), marker=:circle)
Colorbar(fig[2,1][1,2], s, label="log10(ρ)")

ax4 = Axis(fig[3, 1][1,1], xlabel="Prey", ylabel="Predator")
s=scatter!(ax4, tr_plot[:,1], tr_plot[:, 2], linewidth=1, color=log10.(1 ./ speeds_plot), marker=:circle)
Colorbar(fig[3,1][1,2], s, label="log10(1/|xdot|)")

# ax5 = Axis(fig[2, 2], xlabel="density of points ρ", ylabel="PDF(ρ)")
# hist!(ax5, log.(ms), bins=numbins, normalization=:pdf, color=:black)
ax5 = Axis(fig[2, 2], xlabel="density of points ρ", ylabel="PDF(ρ)", yscale=log10, xscale=log10)
# hist!(ax5, ms, bins=bins, normalization=:pdf, color=:black, offset=1e-8)
weights, bins = histogram(ms, numbins)
weights .+= 1
lines!(ax5, bins[1:end-1], weights, color=:black)
ax6 = Axis(fig[3, 2], xlabel="1/|xdot|", ylabel="PDF(1/|xdot|)", yscale=log10, xscale=log10)
# ax6 = Axis(fig[3, 2], xlabel="1/|xdot|", ylabel="PDF(1/|xdot|)",  xscale=log10)
# hist!(ax6, log.(1 ./ speeds), bins=numbins, normalization=:pdf, color=:black)
weights, bins = histogram(1 ./ speeds, numbins)
lines!(ax6, bins[1:end-1], weights, color=:black)

save("$(plotsdir())/$(dirname)/predatorprey-$(typename)-α_$(α)-K_$(K)-numdigitsmeasure_$(numdigits)-intstep_$(Δt)-pdf.png", fig)
# save("$(plotsdir())/$(typename)/homoclinicbifurcation/predatorprey-$(typename)-α_$(α)-K_$(K)-numdigitsmeasure_$(numdigits)-intstep_$(Δt)-pdf.png", fig)
# xlims!(ax1, Ttr, T)
# xlims!(ax2, Ttr, T)
# xlims!(ax3, Ttr, T)


# Looking for homoc bif
using Statistics
# Ks = [5, 10, 15, 18, 19, 20, 21, 22, 25, 30]
logrange(x1, x2; length) = (10^y for y in range(log10(x1), log10(x2), length=length))
# Ks = collect(logrange(5, 30; length=20))
Ks = range(2, 30; length=100)
T=5e4
speedstats = zeros(Float64, (length(Ks), 2))
msstats = zeros(Float64, (length(Ks), 2))
posstats = zeros(Float64, (length(Ks), 4))
for (idx, K) ∈ enumerate(Ks)
	u0 = [0.1, 0.1]
	p = [α, γ, ϵ, ν, h, K, m]
	solver = Tsit5()
	tspan = (0, T)
	prob = ODEProblem(predator_prey!, u0, tspan, p)
	sol = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8)
	t = sol.t 
	tr = sol[:,:]'
	
	speeds = norm.(timederivative(sol))
	
	numdigits = 5
	ms = histmeasure(sol', numdigits);
	speedstats[idx, :] = [mean(1 ./ speeds), maximum(1 ./ speeds)]
	msstats[idx, :] = [mean(ms), maximum(ms)]
	posstats[idx, :] = [minimum(tr[:,1]), maximum(tr[:,1]), minimum(tr[:,2]), maximum(tr[:,2])]
end

yscale=log 
c1 = :black
c2 = :purple
c3 = :green
c4 = :cyan
fig=Figure(resolution=(800, 800), fontsize = 20) 
ax = Axis(fig[1,1], yscale=log10, ylabel="mean(ρ)")
scatterlines!(ax, Ks, msstats[:,1], color=c1)
ax = Axis(fig[1,2], yscale=log10, ylabel="max(ρ)")
scatterlines!(ax, Ks, msstats[:,2], color=c1)
ax = Axis(fig[2,1], yscale=log10, ylabel="mean(1/|xdot|)", xlabel="K")
scatterlines!(ax, Ks, speedstats[:,1], color=c2)
ax = Axis(fig[2,2], yscale=log10, ylabel="mean(1/|xdot|)", xlabel="K")
scatterlines!(ax, Ks, speedstats[:,2], color=c2)
ax = Axis(fig[3,1], yscale=log10, ylabel="minimum(x(t))", xlabel="K")
scatterlines!(ax, Ks, abs.(posstats[:,1]), color=c3)
ax = Axis(fig[3,2], yscale=log10, ylabel="maximum(x(t))", xlabel="K")
scatterlines!(ax, Ks, abs.(posstats[:,2]), color=c3)
ax = Axis(fig[4,1], yscale=log10, ylabel="minimum(y(t))", xlabel="K")
scatterlines!(ax, Ks, abs.(posstats[:,3]), color=c4)
ax = Axis(fig[4,2], yscale=log10, ylabel="maximum(y(t))", xlabel="K")
scatterlines!(ax, Ks, abs.(posstats[:,4]), color=c4)

save("$(plotsdir())/$(typename)/homoclinicbifurcation/predatorprey-$(typename)-scaling-α_$(α)-numdigitsmeasure_$(numdigits)-intstep_$(Δt).png", fig)

#-----------------Analyss
fp1 = [0,0] #origin

Nfp = h*m/(ν*γ - m)
Pfp = (Nfp+h)*α*(1-Nfp/K)/γ
fp2 = [Nfp, Pfp]

fp3 = [K, 0]

using Symbolics, LinearAlgebra
function get_jacobian_predatorprey(fp, p)
	@variables N P 
	@variables α γ ϵ ν h K m
	J = Symbolics.jacobian([α*N*(1 - N/K) - γ*N*P / (N+h), ϵ * (  ν*γ*N*P/(N+h) - m*P    )], [N, P]; simplify=true)
	J_fp = substitute.(J, (Dict(α=>p[1], γ=>p[2], ϵ=>p[3], ν=>p[4], h=>p[5], K=>p[6], m=>p[7], N=>fp[1], P=>fp[2]),))
	J_fp = Symbolics.value.(J_fp)
end

##origin
J = get_jacobian_predatorprey(fp1, p)
eigvals, eigvecs = eigen(J)

J = get_jacobian_predatorprey(fp3, p)
eigvals, eigvecs = eigen(J)

J = get_jacobian_predatorprey(fp2, p)
eigvals, eigvecs = eigen(J)







# ---------------------------- frequency analysis ---------------------------- #
include("$(scriptsdir())/utils.jl")
using FFTW
using GLMakie
typename = "crawlby"
ϵ = 1.0
T = 2e2
Ttr= 100
Δt = 0.01
Tplot = 300
#state space  fig
α =1.5 
K = 10
numbins = 80
# time-series figs
# α = 0.8 
# K = 15 
dirname=typename
u0 = [0.1, 0.1]
p = [α, γ, ϵ, ν, h, K, m]
# pp = ContinuousDynamicalSystem(predator_prey, u0, p)
# tr = trajectory(pp, T; Δt, Ttr); t = 0:Δt:T;
solver = Tsit5()
tspan = (0, T)
prob = ODEProblem(predator_prey!, u0, tspan, p)
sol = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8)
ts = sol.t 
tr = sol[:,:]'
Ts = Δt

axs = []
fig = Figure()
for i=1:2
	signal = tr[:,i]
	variable = ["x", "y"][i]
	F, freqs, periods, maxperiod, Fmax = getspectrum2(signal, ts, Ts)
	ax1 = Axis(fig[1,i], title="Tmax = $maxperiod", ylabel="$variable", xlabel="t")
	time_domain = lines!(ax1, ts, signal, color=:black)
	# F .+= 1
	# ax2 = Axis(fig[2,i], ylabel="FFT($variable)", xlabel="Frequency", yscale=log10)
	ax2 = Axis(fig[2,i], ylabel="FFT($variable)", xlabel="Frequency")
	freq_domain = lines!(ax2, freqs, abs.(F), color=:black) 
	scatter!(ax2, [1/maxperiod], [Fmax], color=:red)
	push!(axs, [ax1, ax2]...)
	# ylims!(ax2, 1e-4, 1.3*Fmax)
end
linkxaxes!(axs[1], axs[3])
linkaxes!(axs[2], axs[4])
save("$(plotsdir())/$(savedir)/crawlby-spectrumanalysis.png")







#---------- MOVIE 
# u0 = [0.1, 0.1]
using ColorSchemes
u0 = [5, 0.1]
# α=0.1 #very slow, ̢o~-1
α=0.1 #very slow, ̢o~-1
p = [α, γ, ϵ, ν, h, K, m]
solver = Tsit5();
tspan = (0, T);
Ttr=1e4; T = Ttr + 200 #for α=0.1
Δt = 0.01
prob = ODEProblem(predator_prey!, u0, tspan, p)
sol = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8);
t = sol.t; tr = sol[:,:]';

Tplot = 44.1
Tplot = 200
# framerate=450
framerate=900
Δtanimation = 0.01
# Δtanimation = 1
t_plot, tr_plot, speeds_plot, frames = animationdata(sol, Tplot, Δt)
colors = pointspeed_as_colors(speeds_plot);

points = Observable(Point2f[(tr[1,1], tr[1,2])])
colors_ob = Observable([colors[1]])
time = Observable(t_plot[1])
fig = Figure(resolution=(800, 600))
ax = Axis(fig[1,1], ylabel="Predator", xlabel="Prey",  title= @lift("t = $((round($time, digits=0)))"));
# fig, ax = scatter(points, color=colors_obg)
scatter!(ax, points, color=colors_ob)
hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
limits!(ax, -0.3, 10.5, -0.15, 3.5)
hidespines!(ax, :t, :r) 
record(fig, "$(plotsdir())/$(dirname)/crawlby_animation-alpha_$(α).mp4", frames;
        framerate) do frame
    new_point = Point2f(tr_plot[frame,1], tr_plot[frame,2])
    points[] = push!(points[], new_point)
	colors_ob[] = push!(colors_ob[], colors[frame])
    time[] = t_plot[frame]
end





### ------for limit cycle 
using ColorSchemes
T=113
Ttr=100
Δt = 0.01
u0 = [5, 0.1]
p = [α, γ, ϵ, ν, h, K, m]
solver = Tsit5()
tspan = (0, T)
prob = ODEProblem(predator_prey!, u0, tspan, p)
sol = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8)
t = sol.t 
tr = sol[:,:]'
speeds = norm.(timederivative(sol))
Tplot= T-Ttr
framerate=200
tplot = Int64(Tplot/Δt);
Δt_plot = Int64(Δt/Δt)
t_plot = t[1:Δt_plot:tplot];
tr_plot = tr[1:Δt_plot:tplot, :];
speeds_plot = log10.(1 ./ speeds[1:Δt_plot:tplot])
frames = 2:length(t_plot)

#transform the vector with the info for the colors onto a Int vector going from 1 to 264; this is used to index the colormap (wihch has 264 colors); basically transforming it into an vector of indices
v = (speeds_plot .- minimum(speeds_plot)) ./ (maximum(speeds_plot) .- minimum(speeds_plot)) .* 255 .+ 1
v = round.(Int, v )
colors = ColorSchemes.viridis[v]

points = Observable(Point2f[(tr[1,1], tr[1,2])])
colors_ob = Observable([colors[1]])
t_title = Observable(t_plot[1])
fig = Figure(resolution=(800, 600))
ax = Axis(fig[1,1], ylabel="Predator", xlabel="Prey", title= @lift("t = $($t_title)"))
# fig, ax = scatter(points, color=colors_obg)
scatter!(ax, points, color=colors_ob)
hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
limits!(ax, 0, 1.25, 0.25, 1.1)
hidespines!(ax, :t, :r) 
record(fig, "$(plotsdir())$(dirname)/$(typename)_animation.mp4", frames;
        framerate) do frame
	t_title[] = t_plot[frame]
    new_point = Point2f(tr_plot[frame,1], tr_plot[frame,2])
    points[] = push!(points[], new_point)
	colors_ob[] = push!(colors_ob[], colors[frame])
end


### -----------fast-slow system
using ColorSchemes
T=820
Ttr=500
Δt = 0.05
# u0 = [5, 0.1]
p = [α, γ, ϵ, ν, h, K, m]
solver = Tsit5()
tspan = (0, T)
prob = ODEProblem(predator_prey!, u0, tspan, p)
sol = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8)
t = sol.t 
tr = sol[:,:]'
speeds = norm.(timederivative(sol))
Tplot= T-Ttr
framerate=350
tplot = Int64(Tplot/Δt);
Δt_plot = Int64(Δt/Δt)
# Δt_plot = Int64(1.0/Δt)
t_plot = t[1:Δt_plot:tplot];
tr_plot = tr[1:Δt_plot:tplot, :];
speeds_plot = log10.(1 ./ speeds[1:Δt_plot:tplot])
frames = 2:length(t_plot)

#transform the vector with the info for the colors onto a Int vector going from 1 to 264; this is used to index the colormap (wihch has 264 colors); basically transforming it into an vector of indices
v = (speeds_plot .- minimum(speeds_plot)) ./ (maximum(speeds_plot) .- minimum(speeds_plot)) .* 255 .+ 1
v = round.(Int, v )
colors = ColorSchemes.viridis[v]

points = Observable(Point2f[(tr[1,1], tr[1,2])])
colors_ob = Observable([colors[1]])
t_title = Observable(t_plot[1])
fig = Figure(resolution=(800, 600))
ax = Axis(fig[1,1], ylabel="Predator", xlabel="Prey", title= @lift("t = $($t_title)"))
# fig, ax = scatter(points, color=colors_obg)
scatter!(ax, points, color=colors_ob)
hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
limits!(ax, -0.02, 1.42, 0.5, 0.73)
hidespines!(ax, :t, :r) 
record(fig, "$(plotsdir())$(dirname)/$(typename)_animation.mp4", frames;
        framerate) do frame
	t_title[] = t_plot[frame]
    new_point = Point2f(tr_plot[frame,1], tr_plot[frame,2])
    points[] = push!(points[], new_point)
	colors_ob[] = push!(colors_ob[], colors[frame])
end
