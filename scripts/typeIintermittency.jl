using DrWatson
@quickactivate "metastability"

# using GLMakie
using DynamicalSystems, Statistics
include("$(scriptsdir())/utils.jl")

# ---------------------------------------------------------------------------- #
#                                   LOGISTIC                                   #
# ---------------------------------------------------------------------------- #
reduced_r(r, rc) = abs(r-rc)/rc
reduced_to_normal(μ, rc) = rc*(1-μ)
rc = 1+√8
savedir="typeI"
# ------------------------------- Phenomenology ------------------------------ #
# Cobweb interactive
# the second range is a convenience for intermittency example of logistic
# rrange = 1:0.001:4.0
# rrange = (rc = 1 + sqrt(8); [rc, rc - 1e-5, rc - 1e-3])
# rrange = 3.8248:0.00001:3.8249
# rrange = 3.7:0.0001:3.845
# using InteractiveDynamics,
# rrange = 3.8:0.00001:3.84
# lo = Systems.logistic(0.4; r = rrange[1])
# interactive_cobweb(lo, rrange, 5)

#Time series
T = 1000
T = 5000
Ttr = 0
r = 3.8282
lo = Systems.logistic(0.4; r); t = Ttr:Ttr+T
tr = trajectory(lo, T; Ttr)
fig = Figure(resolution=(800, 300), fontsize=30)
ax = Axis(fig[1,1], ylabel="x", xlabel="t")
scatterlines!(ax, t, tr, color=:black, markersize=6)
ylims!(0, 1)
xlims!(0, 150)
save("$(plotsdir())/$(savedir)/logistic-timeseries-r_$(r)-uncoloured.png", fig)

# ------------------ Estimate dwell time on laminar periods ------------------ #
"""
Finds laminar periods by identifying the phases in which the windowed standard deviation
of the f^3 map is small (er than a threshold). In the laminar period, f^3 is a fp,
constant, so std is very small. In chaos, it's nonconstant so higher. Works quite well,
but not perfectly, and needs some fine tuning of the parameters.
"""
function laminarperiods(tr; win_size=4, std_th)
    std_tr3 = moving_std(tr, 3*win_size, 3);
    laminar_period_bool = std_tr3 .< std_th #1s are laminar, 0s are nonlaminar
end
# # bool_fp_order = repetition_every_order([4, 1, 2, 3, 1,2,3, 1,2, 3, 3, 2, 1, 3, 2, 1], 3)
# a,b=length_samevalues_allowfluctuations(bool_fp_order, 0)
"""
Finds laminar periods by identifying the phases that repeat every order=3 iterations
"""
function laminarperiods(tr, order; kwargs...)
    bool_fp_order = repetition_every_order(tr, order; kwargs...)
    T, bool_laminarperiods = length_samevalues_allowfluctuations(bool_fp_order, order)
    return T[1], bool_laminarperiods
end

function estimate_laminarperiod_duration(tr; win_size=4, std_th, num_allowed_fluctuations=3)
    laminar_period_bool = laminarperiods(tr; win_size, std_th)
    T, corr_lpb = length_samevalues_allowfluctuations(laminar_period_bool, num_allowed_fluctuations)
    return T[1], corr_lpb ##laminar period only for 1
end

function chaoticperiods(tr, order; kwargs...)
    bool_fp_order = repetition_every_order(tr, order; kwargs...)
    T, bool_laminarperiods = length_samevalues_allowfluctuations(bool_fp_order, order)
    return T[0], bool_laminarperiods
end

T = 1000
T = 15000
Ttr = 0

#coloring the laminar and chaotic periods; not very good though
# win_size = 2;# win_size = 5
# std_th = 0.05;

r = 3.8284
# r = rc-1e-7
lo = Systems.logistic(0.4; r); t = Ttr:Ttr+T
traj = trajectory(lo, T; Ttr);
# lam_periods_bool = laminarperiods(tr; win_size, std_th); tr = tr[1:end-3*win_size+1]; t = t[1:end-3*win_size+1];
# Ts, corr_lam_periods_bool = estimate_laminarperiod_duration(tr; win_size, std_th, num_allowed_fluctuations=5)

Ts, corr_lam_periods_bool= laminarperiods(traj, 3; atol=0.0, rtol=0.02)
t = t[1:length(corr_lam_periods_bool)];
traj = traj[1:length(corr_lam_periods_bool)]
fig = Figure(resolution=(800, 300), fontsize=30, figure_padding=(5, 30, 5, 30))
ax = Axis(fig[1,1], ylabel="x", xlabel="t")
# colors = [el > 0 ? :red : :black for el in lam_periods_bool];
colors = [el > 0 ? :green : :purple for el in corr_lam_periods_bool];
scatterlines!(ax, t, traj, color=colors)
ylims!(0, 1)
xlims!(0, 500)
# save("../plots/typeI/logistic-timeseries-r_$(r)-coloured-winsize_$(win_size)-stdth_$(std_th).png", fig)
# save("$(plotsdir())/$(savedir)/logistic-timeseries-r_$(r)-colouredphases-patternmatching.png", fig, px_per_unit=3)
save("$(plotsdir())/$(savedir)/logistic-timeseries-r_$(r)-colouredphases-patternmatching-stricter-rtol_$(0.02).png", fig, px_per_unit=3)
# save("../plots/typeI/logistic-timeseries-r_$(r)-coloured-winsize_$(win_size)-stdth_$(std_th)-examplesmalllaminarperiodsinsidechaos.png", fig)




# --------------------------------------------------------- Distribution of laminar period --------------------------------------------------------- #
Ttr = 500
T = 5e6
# win_size = 4; std_th = 0.05;

r = rc-1e-7
for r ∈ [rc-1e-3, rc-1e-4, rc-1e-5, rc-1e-6, rc-1e-7, rc-1e-8]
    lo = Systems.logistic(0.4; r)
    μ = reduced_r(r, rc)
    est_cutoff = μ^(-1/2)

    u0s = sort(rand(10))
    Ts_all = []
    Ts_chaos_all = []
    for (idx,u0) ∈ enumerate(u0s)
        tr = trajectory(lo, T, u0; Ttr); t_tr = Ttr:1:Ttr+T
        Ts, corr_lam_periods_bool= laminarperiods(tr, 3; atol=0.0, rtol=0.02)
        Ts_chaos, corr_lam_periods_bool= chaoticperiods(tr, 3; atol=0.0, rtol=0.02)
        Ts_all = vcat(Ts_all, Ts)
        Ts_chaos_all = vcat(Ts_chaos_all, Ts_chaos)
    end
    fig = Figure(resolution=(500, 300), fontsize=30, figure_padding=(5, 30, 5, 30))
    ax = Axis(fig[1,1], title="|r-rc| = $(round(abs(r-rc), digits=10)), μ^(-1/2)=$(floor(Int, est_cutoff))", ylabel="P(τ)", xlabel="τ")
    hist!(ax, Ts_all; bins=500, color=:black)
    save("$(plotsdir())/$(savedir)/distributions-logistic/logistic-distribution-laminartimes-r_$(r)-estimatelaminarperiodviadirectcomparison.png", fig)

    fig = Figure(resolution=(500, 300), fontsize=30, figure_padding=(5, 30, 5, 30))
    ax = Axis(fig[1,1], title="|r-rc| = $(round(abs(r-rc), digits=10))", ylabel="P(τchaos)", xlabel="τchaos")
    hist!(ax, Ts_chaos_all; bins=500, color=:black)
    save("$(plotsdir())/$(savedir)/distributions-logistic/logistic-distribution-chaotictimes-r_$(r)-estimatelaminarperiodviadirectcomparison.png", fig)
end






# -------------------------------------------------------- Scaling of laminar period with r -------------------------------------------------------- #
logrange(x1, x2; length) = (10^y for y in range(log10(x1), log10(x2), length=length))

Ttr = 500
T = 1000000
# win_size = 4
# std_th = 0.001;
μs = collect(logrange(1e-9, 1e-5; length=30))
rs = reduced_to_normal.(μs, rc)

u0 = rand();
T_means = zeros(Float64, length(rs)); T_chaos_means = zeros(Float64, length(rs));
T_maxs = zeros(Float64, length(rs)); T_chaos_maxs = zeros(Float64, length(rs));
for (i, r) ∈ enumerate(rs)
    lo = Systems.logistic(u0; r)
    tr = trajectory(lo, T, u0; Ttr); t_tr = Ttr:1:Ttr+T
    Ts, corr_lam_periods_bool= laminarperiods(tr, 3; atol=0.0, rtol=0.02)
    Ts_chaos, corr_lam_periods_bool= chaoticperiods(tr, 3; atol=0.0, rtol=0.02)
    T_means[i] = mean(Ts);  T_chaos_means[i] = mean(Ts_chaos)
    T_maxs[i] = maximum(Ts);  T_chaos_maxs[i] = maximum(Ts_chaos)
end

using CurveFit
offset, exponent = linear_fit(log10.(μs), log10.(T_means))
fig = Figure(resolution=(500, 300), fontsize=30, figure_padding=(5, 30, 5, 30))
ax = Axis(fig[1,1], yscale=log10, xscale=log10, title = "<τ> ~ μ^{$(round(exponent, digits=3))}")
scatter!(ax, μs, T_means, color=:black)
lines!(ax, μs, μs .^(exponent) .* 10^offset, color=:red )
ax.ylabel="<τ>"; ax.xlabel="μ"
save("$(plotsdir())/$(savedir)/logistic-scaling-with-r-estimatelaminarperiodviadirectcomparison.png", fig)

fig = Figure(resolution=(500, 300), fontsize=30, figure_padding=(5, 30, 5, 30))
ax = Axis(fig[1,1], xscale=log10)
scatter!(ax, μs, T_chaos_means, color=:black)
ax.ylabel="<τchaos>"; ax.xlabel="μ"
save("$(plotsdir())/$(savedir)/logistic-scaling-with-r-estimatechaoticperiodviadirectcomparison.png", fig)

fig = Figure(resolution=(500, 300), fontsize=30, figure_padding=(5, 30, 5, 30))
ax = Axis(fig[1,1], xscale=log10, yscale=log10)
scatter!(ax, μs, T_means./T_chaos_means, color=:black)
ax.ylabel="<τ>/<τchaos>"; ax.xlabel="μ"
save("$(plotsdir())/$(savedir)/logistic-scaling-with-r-ratio-laminar-to-chaotic.png", fig)


fig = Figure(resolution=(500, 300), fontsize=30, figure_padding=(5, 30, 5, 30))
ax = Axis(fig[1,1], xscale=log10, yscale=log10)
scatter!(ax, μs, T_maxs./T_chaos_maxs, color=:black)
ax.ylabel="τmax/τchaosmax"; ax.xlabel="μ"
save("$(plotsdir())/$(savedir)/logistic-scaling-with-r-ratio-maximumduration-laminar-to-chaotic.png", fig)

# length_samevalues([1,1,1,0,0,0,1,1,1,0,1,0,1,1,0,0,0,0])
# length_samevalues_allowfluctuations([1,1,1,0,0,0,1,1,1,0,1,0,1,1,0,0,0,0,0,0,0], 1)


# ---------- Below: tried a few different ways: through finite-time lyapunovs (tried mean and min, its hard because they fluctuate a lot); --------- #
#  and finally with std of f^3. This worked best.
#=

fig = Figure()
ax = Axis(fig[1,1])
lines!(t, λt, color=λt .> 0, colormap=:bluesreds)
ax2 = Axis(fig[2,1])
lines!(t_tr, tr)

using DynamicalSystemsBase:DDS
function lyapunovspectrum_instantaneous(ds::DDS{false, T, 1}, N; Ttr = 0) where {T}
    x = get_state(ds); f = ds.f
    p = ds.p; t0 = ds.t0
    if Ttr > 0
        for i in t0:(Ttr+t0)
            x = f(x, p, i)
        end
    end
    t = (t0+Ttr):(t0+Ttr+N)
    λs = zeros(T, length(t));
    for idx = 1:length(t)
        i = t[idx]
        x = f(x, p, i)
        @inbounds λs[idx] = log(abs(ds.jacobian(x, p, i)))
    end
    t = collect(t)
    return λs, t
end


"""
moving average in window of size ws, walking through each element (so offset is 1)
"""
function moving_average(v, ws)
    mv_average = zeros(length(v)-(ws-1))
    for i=1:length(mv_average)
        mv_average[i] = mean(v[i:i+ws-1])
    end
    return mv_average
end


function moving_maximum(v, ws)
    mv_average = zeros(length(v)-(ws-1))
    for i=1:length(mv_average)
        mv_average[i] = maximum(v[i:i+ws-1])
    end
    return mv_average
end

tr = trajectory(lo, T; Ttr); t_tr = Ttr:1:Ttr+T
λt, t = lyapunovspectrum_instantaneous(lo, T; Ttr)
win_size = 12
λt_fs = moving_average(λt, win_size);
# λt_fs = moving_maximum(λt, win_size);
t_fs = t[1:end-(win_size-1)]
t_tr = t_tr[1:end-(win_size-1)]
tr = tr[1:end-(win_size-1)]

# moving_average([1,1,1,2,2,2,3,1,1,1,1], 3)
λ_th = 0.18

fig = Figure()
ax = Axis(fig[1,1])
# lines!(t, λt, color=λt .> 0, colormap=:bluesreds)
lines!(ax, t_fs, λt_fs, color=λt_fs .> λ_th, colormap=:bluesreds)
hlines!(ax, λ_th, color=:black, linestyle="--")
ax2 = Axis(fig[2,1])
# scatterlines!(t_tr, tr, color=λt_fs .> 0, colormap=:bluesreds)
# scatterlines!(t_tr, tr, color=λt_fs .> 0, colormap=[:red, :blue])
lines!(ax2, t_tr, tr, color=λt_fs .> λ_th, colormap=[:blue, :red])
linkxaxes!(ax, ax2)

## trying with f3

# tr3 = tr[1:3:end]
# t_tr3 = t[1:3:end]
lo = Systems.logistic(0.4; rs[end])
Ttr = 500
T = 5000
win_size = 4
tr = trajectory(lo, T; Ttr); t_tr = Ttr:1:Ttr+T;
std_tr3 = moving_std(tr, 3*win_size);
t_stdtr3 = t_tr[1:end-3*win_size+1];
t_tr = t_tr[1:end-(3*win_size-1)];
tr = tr[1:end-(3*win_size-1)];

# std_th = 0.0006; winsize = 2
std_th = 0.001;
fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, t_stdtr3, std_tr3, color=std_tr3 .< std_th, colormap=[:red, :blue])
hlines!(ax, std_th, color=:black, linestyle="--")

ax2 = Axis(fig[2,1])
lines!(ax2, t_tr, tr, color=std_tr3 .< std_th, colormap=[:red, :blue])
linkxaxes!(ax, ax2)

laminar_period_bool = std_tr3 .< std_th #1s are laminar, 0s are nonlaminar
T = length_samevalues_allowfluctuations(laminar_period_bool, 10)[1] #laminar period only for 0


# scaling of laminar period
r = 3.83
x = 0:0.01:1
logistic_derivative = lo.jacobian.(x, r, 1)
lines(x, logistic_derivative)
=#

#---------------------phase space discretization and count density of points --------------------------------------------------------------
include("utils.jl")
Δx = 0.01
statespace_boxes = 0:Δx:1
Ttr = 500
T = 100000
# μ = 1e-8 #mainly period
μ = 1e-5 #mainly period
r = reduced_to_normal(μ, rc)
u0 = 0.3;
lo = Systems.logistic(u0; r)
tr = trajectory(lo, T, u0; Ttr); ts = Ttr:1:Ttr+T

#won't really work for logistic
# box_durs, boxes_ext = durationinboxes(tr, collect(statespace_boxes))
# mean_box_durs = mean.(box_durs)
#
Δx = 0.001
statespace_boxes = 0:Δx:1
# ρ = densitypointsintraj(tr, Int(abs(log10(Δx))), statespace_boxes)
ρ_pts = densitypointsintraj_inbins(tr, Int(abs(log10(Δx))), statespace_boxes)
ρ_pts = densitypointsintraj_inpoints(tr, Int(abs(log10(Δx))))
colors = vector_to_colors(ρ_pts)

fig = Figure(resolution=(800, 400))
ax = Axis(fig[1,1], ylabel="ρ(x)", xlabel="x")
scatterlines!(statespace_boxes, ρ_pts, color=:black)
ax = Axis(fig[2,1], ylabel="x", xlabel="t")
scatterlines!(ts, tr, color=colors)
xlims!(ax, Ttr, Ttr+500)


# ------------------- recurrence analysis ------------------------
include("$(scriptsdir())/utils.jl")
using CairoMakie
Ttr = 550
T = 200
# μ = 1e-8 #mainly period
μ = 1e-5 #mainly period
r = reduced_to_normal(μ, rc)
u0 = 0.3;
lo = Systems.logistic(u0; r)
traj = trajectory(lo, T, u0; Ttr); ts = Ttr:1:Ttr+T
ϵ=0.1
ϵ=0.05
fig = Figure(resolution=(1000, 500), fontsize=30, figure_padding=(5, 30, 5, 30))
plot_RM!(fig, ts, traj, ϵ)
fig
save("$(plotsdir())/$(savedir)/logistic-recurrenceplot-r_$(r).png", fig, px_per_unit=3)

# statistics of recurrence matrix
RM = RecurrenceMatrix(traj, ϵ)
a = RM[1:10, 1:10]
rs = recurrencestructures(a; lmin=1, theiler=0)

tchaos = 98:145
RM_chaos = RM[tchaos, tchaos]
rs = recurrencestructures(RM_chaos; lmin=1, theiler=0)
hist = rs["recurrencetimes"]
fig = Figure()
ax = Axis(fig[1,1])
scatterlines!(ax, 1:length(hist), hist, color=:black)
hlines!(ax, mean(hist))
fig

# ---------------------------- frequency analysis ---------------------------- #
include("utils.jl")
using FFTW
Ttr = 550
T = 200
Δt = 1
# μ = 1e-8 #mainly period
μ = 1e-5 #mainly period
r = reduced_to_normal(μ, rc)
u0 = 0.3;
lo = Systems.logistic(u0; r)
tr = trajectory(lo, T, u0; Ttr); ts = Ttr:Δt:Ttr+T
signal = tr
Ts = Δt = 1

F, freqs, periods, maxperiod, Fmax = getspectrum2(signal, ts, Ts)


using GLMakie
fig = Figure()
ax = Axis(fig[1,1], title="Tmax = $maxperiod", ylabel="x", xlabel="t")
time_domain = lines!(ax, ts, signal, color=:black)
ax = Axis(fig[2,1], ylabel="FFT(x)", xlabel="Frequency")
freq_domain = lines!(ax, freqs, abs.(F), color=:black)
scatter!(ax, [1/maxperiod], [Fmax], color=:red)
save("$(plotsdir())/$(savedir)/logistic-spectrumanalysis-r_$(r).png")



# ---------------------------------------------------------------------------- #
#                                    LORENZ                                    #
# ---------------------------------------------------------------------------- #
using OrdinaryDiffEq
include("utils.jl")
plotsdir() = "../plots/"

@inbounds @inline function lorenz!(du, u, p, t)
    σ, ρ, β = p
    du[1] = σ*(u[2]-u[1])
    du[2] = u[1]*(ρ-u[3]) - u[2]
    du[3] = u[1]*u[2] - β*u[3]
end

T = 5000; Ttr = 100; Δt = 0.05 #rho = 166.3
# T = 10000; Ttr = 100; Δt = 0.05 #rho = 166.062

u0 = [0.1, 0.1, 0.1]
σ = 10
β = 8/3
# ρ = 166 #period (LC)
# ρ = 166.06 #LC
# ρ = 166.062 #intermittency (CA) GOOD
# ρ = 166.063 #intermittency (CA) GOOD
ρ = 166.1 #intermittency (CA)
# ρ = 166.2 #intermittency (CA)
# ρ = 166.3 #intermittency (CA) GOOD
# ρ = 166.3 #intermittency (CA)
p = [σ, ρ, β]
solver = Tsit5()
tspan = (0, T)
typename = "typeI/lorenz"
prob = ODEProblem(lorenz!, u0, tspan, p)
sol = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, maxiters=1e9)
t = sol.t
tr = sol[:,:]'
speeds = norm.(timederivative(sol))


fig = Figure(resolution=(800,600), fontsize = 20)
ax1 = Axis(fig[1, 1],  ylabel="x", xlabel="t")
# ax2 = Axis(fig[1, 2])
# ax3 = Axis(fig[1, 3])
ax4 = Axis3(fig[2, :])
lines!(ax1, t, tr[:, 1], linewidth=1, color=:black)
# lines!(ax2, t, tr[:, 2], linewidth=1, color=:black)
# lines!(ax3, t, tr[:, 3], linewidth=1, color=:black)
# xlims = (990, 1000) #LC, sharing 3 cols
xlims = (900, 1000) #intermittency at 166.062
xlims!(ax1, xlims[1], xlims[2])
# xlims!(ax2, xlims[1], xlims[2])
# xlims!(ax3, xlims[1], xlims[2])
scatter!(ax4, tr[:,1], tr[:,2], tr[:,3], color= log10.(1 ./ speeds), markersize=5)
save("$(plotsdir())$(typename)/lorenz-x-attractor-rho_$(ρ).png", fig)


##just time series
fig = Figure(resolution=(800,250))
ax1 = Axis(fig[1, 1],  ylabel="x", xlabel="t")
lines!(ax1, t, tr[:, 1], linewidth=1, color=:black)
xlims!(ax1, 780, 840) #for 166.1
save("$(plotsdir())$(typename)/lorenz-timeseriesx-rho_$(ρ).png", fig)




#---------------------------------------movie --------------
Tplot= T-Ttr
framerate=500 #166.3
# framerate=1000 #166.3
# framerate=2000 #166.3
tplot = Int64(Tplot/Δt);
Δt_plot = Int64(1/Δt);
t_plot = t[1:Δt_plot:tplot];
tr_plot = tr[1:Δt_plot:tplot, :];
speeds_plot = log10.(1 ./ speeds[1:Δt_plot:tplot]);
frames = 2:length(t_plot);

az = -1.142477796076938
el = 0.08906585039886639

#transform the vector with the info for the colors onto a Int vector going from 1 to 264; this is used to index the colormap (wihch has 264 colors); basically transforming it into an vector of indices
v = (speeds_plot .- minimum(speeds_plot)) ./ (maximum(speeds_plot) .- minimum(speeds_plot)) .* 255 .+ 1;
v = round.(Int, v );
colors = ColorSchemes.viridis[v];

points = Observable(Point3f[(tr[1,1], tr[1,2], tr[1,3])])
colors_ob = Observable([colors[1]])
fig = Figure(resolution=(800, 600))
ax = Axis3(fig[1,1], azimuth = az, elevation = el)
scatter!(ax, points, color=colors_ob)
hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
limits!(ax, -50, 50, -100, 100, 50, 250)
# hidespines!(ax, :t, :r)
record(fig, "$(plotsdir())$(typename)/lorenz-rho_$(ρ)-faster2.mp4", frames;
        framerate) do frame
    new_point = Point3f(tr_plot[frame,1], tr_plot[frame,2], tr_plot[frame,3])
    points[] = push!(points[], new_point)
	colors_ob[] = push!(colors_ob[], colors[frame])
end


#---------------------------------------movie with grey background--------------
ρ = 166.06 #LC, just before SN
p = [σ, ρ, β]; prob = ODEProblem(lorenz!, u0, tspan, p); sol = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, maxiters=1e9); t = sol.t; traj = sol[:,:]'
limitcycle=deepcopy(traj)

# ρ = 166.1 #intermittency (CA)
ρ = 166.1 #intermittency (CA); a bit too slow
# ρ = 166.2 #intermittency (CA); a bit too slow
# ρ = 166.3 #intermittency (CA)
# ρ = 166.5 #intermittency (CA)
p = [σ, ρ, β]
Ttr = 780; T=850
tspan = (0, T)
Δt = 0.01
prob = ODEProblem(lorenz!, u0, tspan, p)
sol = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, maxiters=1e9); t = sol.t; traj = sol[:,:]'
arewithin = mapslices(x->withinset(x, limitcycle, 2), traj, dims=2)[:,1]; #a bit expensive, have to compare each point to all other points in the set
Tplot= T-Ttr
# framerate=500 #166.3
framerate=100 #166.3
t_plot, tr_plot, speeds_plot, frames = animationdata(sol, Tplot, Δt, Δt)
# colors = pointspeed_as_colors(speeds_plot);
arewithin_plot = mapslices(x->withinset(x, limitcycle, 2), tr_plot, dims=2)[:,1]; #a bit expensive, have to compare each point to all other points in the set

az = -1.142477796076938
el = 0.08906585039886639

points = Observable(Point3f[(tr_plot[1,1], tr_plot[1,2], tr_plot[1,3])])
points2 = Observable(Point2f[(t_plot[1], tr_plot[1,2])])
colors_ob = Observable([colors[1]])
tanim = Observable(t_plot[1])

fig = Figure(resolution=(800, 600))
ax = Axis3(fig[1:2,1], azimuth = az, elevation = el, title= @lift("t = $((round($tanim; digits=0)))"))
scatter!(ax, sol[1,:], sol[2,:], sol[3,:], color=[el == 1 ? (:green,0.8) : (:purple,0.8) for el ∈ arewithin], markersize=3)
scatter!(ax, points, color=(:orange, 1), markersize=15)
# scatter!(ax, tr_plot[:,1], tr_plot[:,2], tr_plot[:,3], color=(:black, 0.2), markersize=4)

hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
limits!(ax, -50, 50, -100, 100, 50, 250)

ax = Axis(fig[3,1], ylabel="x", xlabel="t")
lines!(ax, t_plot, tr_plot[:,1], color=[el == 1 ? :green : :purple for el ∈ arewithin_plot])
scatter!(ax, points2, color=(:orange, 1), markersize=12)
hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
hidespines!(ax, :t, :r)

record(fig, "$(plotsdir())/$(typename)/lorenz-rho_$(ρ)-colorghost.mp4", frames;
        framerate) do frame
    tanim[] = t_plot[frame]
    new_point = Point3f(tr_plot[frame,1], tr_plot[frame,2], tr_plot[frame,3])
    new_point2 = Point2f(t_plot[frame,1], tr_plot[frame,1])
    points[] = [new_point]
    points2[] = [new_point2]
end

# -------- attempt at using lyapunovs exponents to identify the ghost -------- #
using DynamicalSystems
ds = Systems.lorenz(u0; ρ=ρ, σ, β)
tinteg = tangent_integrator(ds, 3)
λs_vec, traj_vec, t= lyapunovspectrum_convergence(tinteg, Int(T/0.01), 0.01, Ttr)
λs = reduce(hcat, λs_vec)'
traj = reduce(hcat, traj_vec)'
# mean(λs, dims=1)
λmaxs = [maximum(λ) for λ in λs_vec]

fig = Figure(resolution=(800, 600))
ax = Axis3(fig[1:2,1], azimuth = az, elevation = el)
# scatter!(ax, traj[:,1], traj[:,2], traj[:,3], color=λs[:,1], markersize=4)
scatter!(ax, traj[:,1], traj[:,2], traj[:,3], color=λmaxs, markersize=4)

ax2 = Axis(fig[3,1], ylabel="x", xlabel="t")
# lines!(ax2, t, traj[:,1], color=λs[:,1])
lines!(ax2, t, traj[:,1], color=λmaxs)
ax3 = Axis(fig[4,1], ylabel="λmax", xlabel="t")
lines!(ax3, t, λs[:,1], color=(:black, 0.4))
lines!(ax3, t, λs[:,2], color=(:black, 0.4))
lines!(ax3, t, λs[:,3], color=(:black, 0.4))
# lines!(ax3, t, λmaxs, color=(:black, 0.4))
linkxaxes!(ax2, ax3)
save("$(plotsdir())/$(typename)/lorenz-rho_$(ρ)-lyapunovs.png", fig)

# -------------- trying to compare direcltty to points in ghost -------------- #
Ttr = 780; T=850
tspan = (0, T)
Δt = 0.01
ρ = 166.06
p = [σ, ρ, β]
prob = ODEProblem(lorenz!, u0, tspan, p)
sol = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, maxiters=1e9)
t = sol.t; traj = sol[:,:]'

fig = Figure(resolution=(800,600), fontsize = 20)
ax = Axis3(fig[1, 1])
scatter!(ax, traj[:,1], traj[:,2], traj[:,3], markersize=5)
limitcycle=deepcopy(traj)

ρ = 166.1
p = [σ, ρ, β]
prob = ODEProblem(lorenz!, u0, tspan, p)
sol = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, maxiters=1e9)
t = sol.t; traj = sol[:,:]'
scatter!(ax, traj[:,1], traj[:,2], traj[:,3], markersize=5)

arewithin = mapslices(x->withinset(x, limitcycle, 2), traj, dims=2)[:,1]
fig = Figure(resolution=(800,600), fontsize = 20)
ax = Axis3(fig[1, 1])
scatter!(ax, traj[:,1], traj[:,2], traj[:,3], markersize=5, color=[el == 1 ? :red : :blue for el ∈ arewithin])