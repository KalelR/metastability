using DrWatson
@quickactivate "metastability"
using GLMakie,  DifferentialEquations

include("$(scriptsdir())/utils.jl")

# """
# U(x) = ax^4/4 - bx 2/2 + cx (c: assymetry paramter); gamma: dissipation, f forcing omega forcing freq, n noise strength
# """
@inbounds function duffing_assymetric_rule(du, u, p, t)
    a,b,c,d,f, ω, n₁, n₂ = p
    x, v = u
    du[1] = v
    du[2] = (-a*x^3 + b*x - c) - d*v + f*cos(ω*t)
end

function noise_duffing(du,u,p,t)
    a,b,c,d,f, ω, n₁, n₂ = p
    du[1] = n₁
    du[2] = n₂
end

a=0.5
b=8.0
c=0.0
d = 0.2
f = 0.0
ω = 1.0
n₁ = n₂ = n = 0.18

T = 5000 
Ttr=0
# T = 15000
Δt = 0.5

u0 = [0.0, 0]

p =  [a, b, c, d, f, ω, n₁, n₂]
tspan = (0, T)
prob_duffing = SDEProblem(duffing_assymetric_rule, noise_duffing, u0, tspan, p; seed=0)
# sol = @time solve(prob_duffing, SOSRA(), saveat=0:Δt:T); t = sol.t #nonstiff, didnt test# vert slow
sol = @time solve(prob_duffing, SKenCarp(), saveat=0:Δt:T); t = sol.t #stiff, worked well but slow, 68secs **BEST
# sol = @time solve(prob_duffing, ImplicitEM(), saveat=0:Δt:T); t = sol.t #stiff, worked well but slow 59secs
# sol = @time solve(prob_duffing, STrapezoid(), saveat=0:Δt:T); t = sol.t #stiff, very fast but didnt have transitions
# sol = @time solve(prob_duffing, ImplicitEulerHeun(), saveat=0:Δt:T); t = sol.t #stiff, 
# sol = @time solve(prob_duffing, SROCK1(), dt=0.01, saveat=0:Δt:T); t = sol.t #stiff, 
# sol = @time solve(prob_duffing, RKMil(), saveat=0:Δt:T); t = sol.t #nonstiff, 
# choice_function(integrator) = (Int(integrator.dt<0.001) + 1)
# alg_switch = StochasticCompositeAlgorithm((EM(), SKenCarp()),choice_function)
# sol = @time solve(prob_duffing, alg_switch, saveat=0:Δt:T); t = sol.t #nonstiff, 


# fig = Figure()
# ax1 = Axis(fig[1, 1], ylabel="x", xlabel="t")
# ax2 = Axis(fig[2, 1], ylabel="U(x)", xlabel="x")
# lines!(ax1, t, sol[1,:], color=[el > 0 ? :red : :black for el in sol[1,:]])

# U(x, a, b, c) = (a/4)*x^4 - (b/2) * x^2 + c*x
# x = -6.5:0.1:6.5
# Us = U.(x, a, b, c)
# lines!(ax2, x, Us, color=:black, linewidth=4)
# scatter!(ax2, sol[1,:], U.(sol[1,:], a, b, c), color=(:orange, 0.5))

# save("duffing-doublwell-a_$(a)-b_$(b)-c_$(c)-d_$(d)-f_$(f)-n_$(n)-solver_SKenCarp-noiseonboth.png", fig)


#just the time series
c1 = "#440154FF"
c2 = "#FDE725FF"
fig = Figure(resolution=(800,250))
ax1 = Axis(fig[1, 1], ylabel="x", xlabel="t")
lines!(ax1, t, sol[1,:], color=[el > 0 ? c1 : c2 for el in sol[1,:]])
save("$(plotsdir())/noise/duffing-timeseries-a_$(a)-b_$(b)-c_$(c)-d_$(d)-f_$(f)-n_$(n)-solver_SKenCarp-viridiscolors.png", fig)

#verifying kramers 
U_peak = U(0, a, b, c)
Uwell, idx = findmin(Us)
Eb = U_peak - U_well
ϵ = 4*n^2
U²(x, a, b, c) = 3a*x^2 - b
xwell = x[idx]
prefactor = 2π/(sqrt(U²(xwell, a, b, c) * abs(U²(0, a, b, c)) ) )
exponent = Eb/ϵ

#------------------------------------------------------------------SCALING
using CurveFit 
laminarperiods(v) = v .> 0 #1 is positive, 0 is negative
T = 1e7
Δt = 0.5
u0 = [0.0, 0]
p =  [a, b, c, d, f, ω, n₁, n₂]
tspan = (0, T)
prob_duffing = SDEProblem(duffing_assymetric_rule, noise_duffing, u0, tspan, p; seed=0)
# u0s = [rand(2) for i=1:10]
# all_τs = []
# @Threads.threads for idx_u0 in 1:length(u0s)
# u0 = u0s[idx_u0]
prob_duffing = SDEProblem(duffing_assymetric_rule, noise_duffing, u0, tspan, p; seed=0)
sol = solve(prob_duffing, SKenCarp(), saveat=0:Δt:T, maxiters=1e9, progress=true); t = sol.t #stiff, worked well but slow
# sol = solve(prob_duffing, EM(), dt=0.01, saveat=0:Δt:T, maxiters=1e9, progress=true); t = sol.t #stiff, worked well but slow
# sol = @time solve(prob_duffing, SROCK1(), dt=0.01, saveat=0:Δt:T, maxiters=1e9); t = sol.t #stiff, 
x = sol[1,:]
τs, _ = length_samevalues_allowfluctuations(laminarperiods(x), 1)
τs_2, _ = length_samevalues_allowfluctuations(laminarperiods(x), 0)
numbins=15
ws, bins = histogram(τs[1], numbins); 


# ensembleprob = EnsembleProblem(prob_duffing)
# sols = solve(ensembleprob, SKenCarp(), EnsembleThreads(), trajectories=1, maxiters=1e9)
# τs_all = Int64[]
# for sol in sols 
#     x = sol[1,:]
#     τs, _ = length_samevalues_allowfluctuations(laminarperiods(x), 1)
#     push!(τs_all, τs[1]...)
# end
# ws, bins = histogram(τs_all, 50); 


# ws, bins = histogram(τs[1], 50); ws .+= 1
# ws2, bins2 = histogram(τs_2[1], 50); ws2 .+= 1
fig = Figure(resolution=(800, 600), fontsize=30,figure_padding=(5, 35, 5, 30))
ax = Axis(fig[1,1], yscale=log10, ylabel="P(τ)", xlabel="τ")
scatterlines!(ax, bins[1:end-1], ws, color=:black)
# scatterlines!(ax, bins2[1:end-1], ws2, color=:orange)

xfit = bins[1:end-1]; yfit = ws[1:end]; a, b = exp_fit(xfit, yfit)
l=lines!(ax, xfit, a .* exp.(b .* xfit), color=:red, label="P(τ) = $(round(a, digits=4)) exp($(round(b, digits=4)) τ)" )
axislegend(ax)

save("$(plotsdir())/noise/duffing-doublwell-a_$(a)-b_$(b)-c_$(c)-d_$(d)-f_$(f)-n_$(n)-solver_SKenCarp-numbins_$(numbins).png", fig)
# save("../plots/noise/duffing-doublwell-a_$(a)-b_$(b)-c_$(c)-d_$(d)-f_$(f)-n_$(n)-solver_SKenCarp-noiseonx.png", fig)


#---------------------------------------------------- Recurrence 
include("utils.jl")
T = 1e3
Δt = 5.0
Δt = 0.5
u0 = [0.0, 0]
p =  [a, b, c, d, f, ω, n]
tspan = (0, T)
prob_duffing = SDEProblem(duffing_assymetric_rule, noise_duffing, u0, tspan, p; seed=0)
sol = solve(prob_duffing, SKenCarp(), saveat=0:Δt:T, maxiters=1e9); t = sol.t #stiff, worked well but slow
traj = Dataset(sol[:,:]')
fig = Figure(resolution = (1000,500))
plot_RM!(fig, sol.t, traj, 0.1; tsmode="lines")
RM = RecurrenceMatrix(traj, ϵ)
rqa(RM)
meanrecurrencetime(RM)
a = recurrencestructures(RM, recurrencetimes=true)

# ---------------------------- Frequency analysis ---------------------------- #
include("utils.jl")
T = 2e3
T=400
Δt = 0.5
u0 = [0.0, 0]
p =  [a, b, c, d, f, ω, n]
tspan = (0, T)
prob_duffing = SDEProblem(duffing_assymetric_rule, noise_duffing, u0, tspan, p; seed=0)
sol = solve(prob_duffing, SKenCarp(), saveat=0:Δt:T, maxiters=1e9); t = sol.t #stiff, worked well but slow
traj = sol[:,:]'; ts = sol.t

signal = traj[140:end,1]; ts = ts[140:end]
Ts = Δt
F, freqs, periods, maxperiod, Fmax = getspectrum2(signal, ts, Ts)

using GLMakie
fig = Figure()
ax = Axis(fig[1,1], title="Tmax = $maxperiod", ylabel="x", xlabel="t")
time_domain = lines!(ax, ts, signal, color=:black)
ax = Axis(fig[2,1], ylabel="FFT(x)", xlabel="Frequency")
freq_domain = lines!(ax, freqs, abs.(F), color=:black) 
scatter!(ax, [1/maxperiod], [Fmax], color=:red)
save("$(plotsdir())/$(savedir)/noisydoublewell-spectrumanalysis-onewell.png")


# τmean = (kt/eb)^2 * ν


#--------------------------------------------------- ANIMATION ---------------------------------
include("$(scriptsdir())/utils.jl")
# plotsdir() = "../plots/"
T = 200
Ttr=0
# T = 800
# T = 900
# Δt = 0.5
Δt = 0.1
u0 = [-3.0, 0]
n₁ = n₂ = n = 0.25
# n₁ = n₂ = n = 0.2
# n₁ = n₂ = n = 0.18
# n₁ = n₂ = n = 0.16
p =  [a, b, c, d, f, ω, n₁, n₂]
tspan = (0, T)
prob_duffing = SDEProblem(duffing_assymetric_rule, noise_duffing, u0, tspan, p; seed=0)
sol = solve(prob_duffing, SKenCarp(), saveat=0:Δt:T, maxiters=1e9, seed=0); t = sol.t #stiff, worked well but slow
t = sol.t ;
traj = sol[:,:]';
speeds = norm.(timederivative(sol));

Tplot= T-Ttr
# framerate=25
# framerate=100
# framerate=150
framerate=75
tplot = Int64(Tplot/Δt);
Δt_plot = round(Int64, Δt/Δt);
t_plot = t[1:Δt_plot:tplot];
tr_plot = traj[1:Δt_plot:tplot, :];
speeds_plot = log10.(1 ./ speeds[1:Δt_plot:tplot]);
frames = 2:length(t_plot)
duration = length(frames) / framerate

#transform the vector with the info for the colors onto a Int vector going from 1 to 264; this is used to index the colormap (wihch has 264 colors); basically transforming it into an vector of indices
v = (speeds_plot .- minimum(speeds_plot)) ./ (maximum(speeds_plot) .- minimum(speeds_plot)) .* 255 .+ 1;
v = round.(Int, v );
colors = ColorSchemes.viridis[v];
doublewell(x, a, b, c) = (a/4)*x^4 - (b/2) * x^2 + c*x
x = -6.5:0.1:6.5
Us = doublewell.(x, a, b, c)


c1 = "#440154FF"
c2 = "#FDE725FF"
cball = "#FF1400"
points = Observable(Point2f[(tr_plot[1,1], doublewell.(tr_plot[1,1], a, b, c))])
colors_ob = Observable([colors[1]])
time = Observable(t_plot[1])
fig = Figure(resolution=(800, 600))
ax = Axis(fig[1,1])
lines!(ax, x, Us, linewidth=4, color=[el > 0 ? c1 : c2 for el in x])
s=scatter!(ax, points, color=cball, markersize=30)
# hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
# limits!(ax, -50, 50, -100, 100, 50, 250)
record(fig, "$(plotsdir())/noise/duffing-doublewell-animation-singlepoint-a_$(a)-b_$(b)-c_$(c)-d_$(d)-f_$(f)-n_$(n)-solver_SKenCarp-faster-viridiscolors.mp4", frames;
        framerate) do frame
    new_point = Point2f(tr_plot[frame,1], doublewell.(tr_plot[frame,1], a, b, c))
    points[] = [new_point]
	# colors_ob[] = push!(colors_ob[], colors[frame])
end


# --------------------------- Well and time-series --------------------------- #

c1 = "#440154FF"
c2 = "#FDE725FF"
cball = "#FF1400"
points = Observable(Point2f[(tr_plot[1,1], doublewell.(tr_plot[1,1], a, b, c))])
# colors_ob = Observable([colors[1]])
time = Observable(t_plot[1])
points2 = Observable(Point2f[(t_plot[1], tr_plot[1,1])])
fig = Figure(resolution=(800, 600))
ax = Axis(fig[1,1], title= @lift("t = $($time)"))
lines!(ax, x, Us, linewidth=4, color=[el > 0 ? c1 : c2 for el in x])
s=scatter!(ax, points, color=cball, markersize=20)


ax = Axis(fig[2, 1], ylabel="x", xlabel="t")
lines!(ax, t_plot, tr_plot[:,1], color=[el > 0 ? c1 : c2 for el in tr_plot[:,1]])
scatter!(ax, points2, color=cball, markersize=20)

# hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
# limits!(ax, -50, 50, -100, 100, 50, 250)
record(fig, "$(plotsdir())/noise/duffing-doublewell-animation-timeseries-singlepoint-a_$(a)-b_$(b)-c_$(c)-d_$(d)-f_$(f)-n_$(n)-solver_SKenCarp-faster-viridiscolors-framerate_$(framerate).mp4", frames;
        framerate) do frame
    time[] = t_plot[frame]
    new_point = Point2f(tr_plot[frame,1], doublewell.(tr_plot[frame,1], a, b, c))
    new_point2 = Point2f(t_plot[frame,1], tr_plot[frame,1])
    points[] = [new_point]
	points2[] = [new_point2]
	# colors_ob[] = push!(colors_ob[], colors[frame])
end