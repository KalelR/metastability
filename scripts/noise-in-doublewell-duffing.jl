using DrWatson
@quickactivate "metastability"
using GLMakie,  DifferentialEquations

include("$(scriptsdir())/utils.jl")

# """
# U(x) = ax^4/4 - bx 2/2 + cx (c: assymetry paramter); gamma: dissipation, f forcing omega forcing freq, n noise strength
# """
@inbounds function duffing_assymetric_rule(du, u, p, t)
    a,b,c,d,f, ω, n = p
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
n₁ = n₂ = n = 0.2

T = 5000 #on v
Ttr=0
# T = 15000
Δt = 0.5

u0 = [0.0, 0]

p =  [a, b, c, d, f, ω, n₁, n₂]
tspan = (0, T)
prob_duffing = SDEProblem(duffing_assymetric_rule, noise_duffing, u0, tspan, p; seed=0)
# sol = solve(prob_duffing, SOSRA(), saveat=0:Δt:T); t = sol.t #nonstiff, didnt test
sol = solve(prob_duffing, SKenCarp(), saveat=0:Δt:T); t = sol.t #stiff, worked well but slow

fig = Figure()
ax1 = Axis(fig[1, 1], ylabel="x", xlabel="t")
ax2 = Axis(fig[2, 1], ylabel="U(x)", xlabel="x")
lines!(ax1, t, sol[1,:], color=[el > 0 ? :red : :black for el in sol[1,:]])

U(x, a, b, c) = (a/4)*x^4 - (b/2) * x^2 + c*x
x = -6.5:0.1:6.5
Us = U.(x, a, b, c)
lines!(ax2, x, Us, color=:black, linewidth=4)
scatter!(ax2, sol[1,:], U.(sol[1,:], a, b, c), color=(:orange, 0.5))

save("duffing-doublwell-a_$(a)-b_$(b)-c_$(c)-d_$(d)-f_$(f)-n_$(n)-solver_SKenCarp-noiseonboth.png", fig)


#just the time series
c1 = "#440154FF"
c2 = "#FDE725FF"
fig = Figure(resolution=(800,250))
ax1 = Axis(fig[1, 1], ylabel="x", xlabel="t")
lines!(ax1, t, sol[1,:], color=[el > 0 ? c1 : c2 for el in sol[1,:]])
save("../plots/noise/duffing-timeseries-a_$(a)-b_$(b)-c_$(c)-d_$(d)-f_$(f)-n_$(n)-solver_SKenCarp-viridiscolors.png", fig)

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
laminarperiods(v) = v .> 0 #1 is positive, 0 is negative
T = 1e8
Δt = 0.5
u0 = [0.0, 0]
p =  [a, b, c, d, f, ω, n₁, n₂]
tspan = (0, T)
# u0s = [rand(2) for i=1:10]
# all_τs = []
# @Threads.threads for idx_u0 in 1:length(u0s)
# u0 = u0s[idx_u0]
prob_duffing = SDEProblem(duffing_assymetric_rule, noise_duffing, u0, tspan, p; seed=0)
sol = solve(prob_duffing, SKenCarp(), saveat=0:Δt:T, maxiters=1e9, progress=true); t = sol.t #stiff, worked well but slow

x = sol[1,:]
τs, _ = length_samevalues_allowfluctuations(laminarperiods(x), 1)
τs_2, _ = length_samevalues_allowfluctuations(laminarperiods(x), 0)

ws, bins = histogram(τs[1], 50); ws .+= 1
ws2, bins2 = histogram(τs_2[1], 50); ws2 .+= 1
fig = Figure()
ax = Axis(fig[1,1], yscale=log10, ylabel="PDF(τ)", xlabel="τ")
scatterlines!(ax, bins[1:end-1], ws, color=:black)
# scatterlines!(ax, bins2[1:end-1], ws2, color=:orange)


using CurveFit 
xfit = bins[2:end-1]
yfit = ws[2:end]
a, b = exp_fit(xfit, yfit)
l=lines!(ax, xfit, a .* exp.(b .* xfit), color=:red, label="y = $(a) exp($(b) x)" )
fig[2,1] = Legend(fig, ax, orientation= :horizontal, tellwidth = false, tellheight = true)

save("../plots/noise/duffing-doublwell-a_$(a)-b_$(b)-c_$(c)-d_$(d)-f_$(f)-n_$(n)-solver_SKenCarp.png", fig)
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
tr = Dataset(sol[:,:]')
fig = Figure(resolution = (1000,500))
plot_RM!(fig, sol.t, tr, 0.1; tsmode="lines")
RM = RecurrenceMatrix(tr, ϵ)
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
tr = sol[:,:]'; ts = sol.t

signal = tr[140:end,1]; ts = ts[140:end]
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
include("utils.jl")
plotsdir() = "../plots/"
T = 400
# Δt = 0.5
Δt = 0.1
u0 = [-3.0, 0]
p =  [a, b, c, d, f, ω, n]
tspan = (0, T)
prob_duffing = SDEProblem(duffing_assymetric_rule, noise_duffing, u0, tspan, p; seed=0)
sol = solve(prob_duffing, SKenCarp(), saveat=0:Δt:T, maxiters=1e9); t = sol.t #stiff, worked well but slow
t = sol.t ;
tr = sol[:,:]';
speeds = norm.(timederivative(sol));

Tplot= T-Ttr
framerate=25
# framerate=50
tplot = Int64(Tplot/Δt);
Δt_plot = round(Int64, Δt/Δt);
t_plot = t[1:Δt_plot:tplot];
tr_plot = tr[1:Δt_plot:tplot, :];
speeds_plot = log10.(1 ./ speeds[1:Δt_plot:tplot]);
frames = 2:length(t_plot)
duration = length(frames) / framerate

#transform the vector with the info for the colors onto a Int vector going from 1 to 264; this is used to index the colormap (wihch has 264 colors); basically transforming it into an vector of indices
v = (speeds_plot .- minimum(speeds_plot)) ./ (maximum(speeds_plot) .- minimum(speeds_plot)) .* 255 .+ 1;
v = round.(Int, v );
colors = ColorSchemes.viridis[v];

# points = Observable(Point2f[(tr[1,1], tr[1,2])])
# colors_ob = Observable([colors[1]])
# fig = Figure(resolution=(800, 600))
# ax = Axis(fig[1,1])
# scatter!(ax, points, color=colors_ob)
# hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
# # limits!(ax, -50, 50, -100, 100, 50, 250)
# record(fig, "../plots/noise/duffing-doublewell-animation-statespace-a_$(a)-b_$(b)-c_$(c)-d_$(d)-f_$(f)-n_$(n)-solver_SKenCarp.mp4", frames;
#         framerate) do frame
#     new_point = Point2f(tr_plot[frame,1], tr_plot[frame,2])
#     points[] = push!(points[], new_point)
# 	colors_ob[] = push!(colors_ob[], colors[frame])
# end

# scatter(tr_plot[:,1], tr_plot[:,2])

doublewell(x, a, b, c) = (a/4)*x^4 - (b/2) * x^2 + c*x
x = -6.5:0.1:6.5
Us = doublewell.(x, a, b, c)


# #all points
# points = Observable(Point2f[(tr_plot[1,1], doublewell.(tr_plot[1,1], a, b, c))])
# colors_ob = Observable([colors[1]])
# fig = Figure(resolution=(800, 600))
# ax = Axis(fig[1,1])
# lines!(ax, x, Us, color=:black, linewidth=4)
# scatter!(ax, points, color=colors_ob)
# # hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
# # limits!(ax, -50, 50, -100, 100, 50, 250)
# record(fig, "../plots/noise/duffing-doublewell-animation-statespace-a_$(a)-b_$(b)-c_$(c)-d_$(d)-f_$(f)-n_$(n)-solver_SKenCarp.mp4", frames;
#         framerate) do frame
#     new_point = Point2f(tr_plot[frame,1], doublewell.(tr_plot[frame,1], a, b, c))
#     points[] = push!(points[], new_point)
# 	colors_ob[] = push!(colors_ob[], colors[frame])
# end
# lines!(ax, t, x, color=[el > 0 ? :red : :black for el in x], label="initial condition on the left")

c1 = "#440154FF"
c2 = "#FDE725FF"
cball = "#FF1400"
points = Observable(Point2f[(tr_plot[1,1], doublewell.(tr_plot[1,1], a, b, c))])
colors_ob = Observable([colors[1]])
time = Observable(t_plot[1])
fig = Figure(resolution=(800, 600))
ax = Axis(fig[1,1], title= @lift("t = $($time)"))
lines!(ax, x, Us, linewidth=4, color=[el > 0 ? c1 : c2 for el in x])
s=scatter!(ax, points, color=cball, markersize=35)
# hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
# limits!(ax, -50, 50, -100, 100, 50, 250)
record(fig, "../plots/noise/duffing-doublewell-animation-singlepoint-a_$(a)-b_$(b)-c_$(c)-d_$(d)-f_$(f)-n_$(n)-solver_SKenCarp-faster-viridiscolors.mp4", frames;
        framerate) do frame
    time[] = frame
    new_point = Point2f(tr_plot[frame,1], doublewell.(tr_plot[frame,1], a, b, c))
    points[] = [new_point]
	# colors_ob[] = push!(colors_ob[], colors[frame])
end