using DrWatson
@quickactivate "metastability"

using GLMakie
include("$(srcdir())/paperplottheme.jl")
using DynamicalSystems, Statistics
include("$(scriptsdir())/utils.jl")
# using CairoMakie;

# ---------------------------------------------------------------------------- #
#                                    LORENZ                                    #
# ---------------------------------------------------------------------------- #

using OrdinaryDiffEq, DelimitedFiles
T = 5000; Ttr = 100; Δt = 0.01;
σ = 10; β = 8/3;
u0 = [0.1, 0.1, 0.1]
solver = Tsit5()

#get the limit cycle (before the bif) and save
ρ = 166.06; p = [σ, ρ, β]; #LC, just before SN
prob = ODEProblem(lorenz!, u0, (0, T), p); sol = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, maxiters=1e9);
sol = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, maxiters=1e9);
limitcycle=deepcopy(sol[:,:]);

#solve with chaotic atttractor
Ttr = 780; T=850;
ρ = 166.1; p = [σ, ρ, β]; #intermittency (CA);
prob = ODEProblem(lorenz!, u0, (0, T), p)
sol = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, maxiters=1e9);
arewithin = mapslices(x->withinset(x, limitcycle', 2), sol[:, :], dims=1)[1, :];
colors = [el == true ? :green : :purple for el in arewithin]

#get the pre-calculated distribution of dwell times in ghost of limit cycle
dwelltimes = readdlm("$(datadir())/typeIintermittency-lorenz-dwelltimes.dat")[:, 1]
numbins = 15;
weights, bins = histogram(dwelltimes, numbins)

azim = 5.475530633326973;
elev = 0.050203673205103286;

fig = Figure(resolution=(columnsize_pt, 1.0*width_pt))
ax1 = Axis(fig[1, 1],  ylabel="x", xlabel="t")
lines!(ax1, sol.t, sol[1, :], color=colors )
xlims!(ax1, 790, 848)
ax1.xticks = 800:20:840
# ax2 = Axis3(fig[2, :], azimuth=azim, elevation=elev, xlabeloffset=20, ylabeloffset=20, zlabeloffset=30, protrusions=0, viewmode=:fit)
# ax2 = Axis3(fig[2, :], azimuth=azim, elevation=elev, xlabeloffset=20, ylabeloffset=20, zlabeloffset=30, protrusions=0, viewmode=:fitzoom) #def
ax2 = Axis3(fig[2, :], azimuth=azim, elevation=elev, xlabeloffset=5, ylabeloffset=5, zlabeloffset=5, protrusions=5, viewmode=:stretch)
# hidedecorations!(ax2); hidespines!(ax2);
hidedecorations!(ax2, label=false, grid=false);
scatter!(ax2, sol[1,:], sol[2,:], sol[3,:], markersize=2, color=colors)
ax3 = Axis(fig[3, 1],  ylabel="PDF(τ)", xlabel="τ", yscale=log10)
# weights .+= 1e-8;
lines!(ax3, bins, weights, color=:green)
rowsize!(fig.layout, 2, Relative(0.6))
ax3.yticks = [1e-3, 1e-2]
# rowgap!(fig.layout, -60)
# save("$(plotsdir())/paper/intermittency-typeI-lorenz.png", fig)
save("$(plotsdir())/paper/intermittency-typeI-lorenz.pdf", fig)


#Get distribution of times
#=
T = 515; Ttr = 500; Δt = 0.01;
ρ = 166.06; p = [σ, ρ, β]; #LC, just before SN
prob = ODEProblem(lorenz!, u0, (0, T), p); sol = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, maxiters=1e9);
sol = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, maxiters=1e9);
limitcycle=deepcopy(sol[:,:]);

Ttr = 300; T=1e5; Δt=0.1;
ρ = 166.1; p = [σ, ρ, β]; #intermittency (CA);
prob = ODEProblem(lorenz!, u0, (0, T), p)
sol = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, maxiters=1e9);
arewithin = mapslices(x->withinset(x, limitcycle', 2), sol[:, :], dims=1)[1, :];
dwelltimes = length_samevalues_allowfluctuations(arewithin, 3)[1][1]
writedlm("$(datadir())/typeIintermittency-lorenz-dwelltimes.dat", dwelltimes)
=#


σ = 10; β = 8/3; u0 = [0.1, 0.1, 0.1]
ρ_chaos = 166.1; p_chaos = [σ, ρ_chaos, β]; #intermittency (CA);
ρ_lc = 166.06; p_lc = [σ, ρ_lc, β]; #LC, just before SN

#dwell times
T = 1e6
T = 1025
dwelltimes_input = Dict("T"=>T, "p_chaos"=>p_chaos, "system_func"=>lorenz!, "u0"=>u0, "p_lc"=>p_lc)
dwelltimes_input = Dict("T"=>T, "p_chaos"=>p_chaos, "system_func"=>lorenz, "u0"=>u0, "p_lc"=>p_lc)
@unpack T, p_chaos, system_func, u0, p_lc = dwelltimes_input

Ttr = 1e3; Δt_chaos = 0.05; #like this to look good
prob = ODEProblem(system_func, u0, (0., T), p_chaos)
u0 =  SVector{3}([0.1, 0.1, 0.1])
sol = solve(prob, Tsit5(), saveat=Ttr:Δt_chaos:T, abstol=1e-8, reltol=1e-8, maxiters=1e9);

#get the limit cycle (before the bif) and save
T = Ttr + 3; Δt = 0.005;
prob = ODEProblem(system_func, u0, (0, T), p_lc); 
sol_lc = solve(prob, Tsit5(), saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, maxiters=1e9);

arewithin = mapslices(x->withinset(x, sol_lc', 10), sol[:, :], dims=1)[1, :];
dwelltimes_int = length_samevalues_allowfluctuations(arewithin, 3)[1] #1 corresponds to time in  the limit cycle (so 1=laminar)
dwelltimes = multiply_dict_number(dwelltimes_int, Δt_chaos)
colors = [el == true ? :green : :red for el in arewithin]

numbins = 10
maximum_time = maximum(maximum.(values(dwelltimes)))
bins = collect(range(0, maximum_time, length=numbins))



fig = Figure()
# for i=1:2
i=1
ax = Axis(fig[i,1])
lines!(ax, sol.t, sol[i, :], color=colors)
lines!(ax, sol_lc.t, sol_lc[i,:])
# end
ax = Axis(fig[i+1, 1])
for (i, (withinlc, dwelltimes_one)) in enumerate(dwelltimes)
weights, _bins = histogram(dwelltimes_one, bins)
if minimum(weights) <= 0 weights .+= 1e-8 end
@show _bins
scatterlines!(ax, _bins, weights, color=[(:red, 1.0), (:green, 1.0)][i], markersize=4, linewidth=0.7)
end


# ---------------------------------------------------------------------------- #
#                                   LOGISTIC                                   #
# ---------------------------------------------------------------------------- #
reduced_r(r, rc) = abs(r-rc)/rc
reduced_to_normal(μ, rc) = rc*(1-μ)
rc = 1+√8
# savedir="typeI"

#Time series
function plot_time_series_logistic()
    T = 5e6
    Ttr = 500
    # r = 3.8282
    r = rc-1e-7
    lo = Systems.logistic(0.4; r); t = Ttr:Ttr+T
    traj_full = trajectory(lo, T; Ttr)
    Tplot=15000
    traj = traj_full[1:Tplot]; t = Ttr:Tplot+Ttr

    Ts, _= logistic_laminarperiods(traj_full, 3; atol=0.0, rtol=0.02);
    hist, edges = histogram(Ts, 500);
    _, corr_lam_periods_bool= logistic_laminarperiods(traj, 3; atol=0.0, rtol=0.02); #for colors
    # c1 = :green; c2 = :purple
    c1 = :purple; c2 = :green
    colors = [el > 0 ? c1 : c2 for el in corr_lam_periods_bool];
    t = t[1:length(corr_lam_periods_bool)]; traj = traj[1:length(corr_lam_periods_bool)];

    fig = Figure(resolution=(columnsize_pt, 0.80*width_pt))
    ax = Axis(fig[1,1], ylabel="x", xlabel="t")
    scatterlines!(ax, t, traj, color=colors, markersize=5, linewidth=1)
    ylims!(ax, 0, 1)
    xlims!(ax, 10700, 10950)
    ax.xticks = [10700, 10800, 10900, 11000]
    # ax2 = Axis(fig[2,1], ylabel="f(x)", xlabel="x")
    # scatter!(ax2, traj_full[1:100:end], lo.f.(traj_full[1:100:end], r, 1), color=:black, markersize=6)
    #Distribution
    # Ts, corr_lam_periods_bool= logistic_laminarperiods(traj_full, 3; atol=0.0, rtol=0.02)
    # Ts_chaos, corr_lam_periods_bool= logistic_chaoticperiods(traj_full, 3; atol=0.0, rtol=0.02)
    # ax3 = Axis(fig[3,1], ylabel="PDF(τ)", xlabel="τ", yscale=log10)
    # hist!(ax3, Ts; bins=500, color=:green)
    # lines!(ax3, edges, hist.+1e-5, color=:green)
    save("$(plotsdir())/mechanisms/intermittency-typeI/logistic-intermittency-typeI.png", fig, px_per_unit=3)
end



#examples showing correctness of my the laminar durations
# fig = Figure(resolution=(columnsize_pt, 0.80*width_pt))
# ax = Axis(fig[1,1], ylabel="x", xlabel="t")
# scatterlines!(ax, t, traj, color=colors, markersize=4)
# xlims!(47000, 47500); save("$(plotsdir())/$(savedir)/logistic-example-shortlaminarperiods-r_$(r).png", fig)
# xlims!(7000, 8000); save("$(plotsdir())/$(savedir)/logistic-example-chaotiburst-between-longlaminarperiods-r_$(r).png", fig)