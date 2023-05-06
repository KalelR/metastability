using DrWatson
@quickactivate

using DifferentialEquations

function fitzhugh_nagumo!(u = 0.5ones(2); a=-0.1, b=0.01, c=0.02, I=0.0, n=0.0)
    ds = ContinuousDynamicalSystem(fitzhugh_nagumo_rule!, u, [a, b, c, I, n])
end
function fitzhugh_nagumo_rule!(du, u, p, t)
    V, w = u
    a, b, c, I, n = p
    du[1] = V*(a-V)*(V-1) -w + I
    du[2] = b*V - c*w
end

# using DynamicalSystems
T = 300; Ttr = 100; Δt = 0.1
u0 = [0.5, 0.5]
a = -0.5; b = 0.1; c = 0.2; I = 0; n = 0 #nice tonic spiking
p = [a, b, c, I, n]
# ds = fitzhugh_nagumo!(u0; a, b, c, I, n)
# tr = trajectory(ds, T; Ttr, Δt); t=Ttr:Δt:T+Ttr
prob = ODEProblem(fitzhugh_nagumo_rule!,  u0, (0, T), p; )
sol = solve(prob, Tsit5(), maxiters=1e9, saveat=Ttr:Δt:T, seed=0);
tr = sol[:,:]'; t = sol.t


using CairoMakie
fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, t, tr[:,1], color=:black)
display(fig)
ax = Axis(fig[2,1])
lines!(ax, tr[:,1], tr[:,2], color=:black)
display(fig)
save("$(plotsdir())/fitzhughnagumo/timeseries-statespace-a_$a-b_$b-c_$c-I_$I-n_$I.png", fig, px_per_unit=3)





# ---------------------------------------------------------------------------- #
#                                 Noisy version                                #
# ---------------------------------------------------------------------------- #
function noise_fitzhugh_nagumo!(du, u,p,t)
    du .= p[5]
end
T = 150; Ttr = 0; Δt = 0.1
u0 = [0.1, 0.1]
# a = -0.5; b = 0.1; c = 0.2; I = -0.3; n = 0.5 #nice tonic spiking
a = -0.5; b = 0.1; c = 0.2; I = -0.2; n = 0.1 #tonic spiking due to noise (almost vertical w-nullcline intersection with minimum of v-nullcline)
# a = -0.0; b = 0.1; c = 0.01; I = 0.5; n = 0.0 #tonic spiking due to LC (not very nice spike though); clearly shows Hopf bif.

p = [a, b, c, I, n]
prob = SDEProblem(fitzhugh_nagumo_rule!, noise_fitzhugh_nagumo!, u0, (0, T), p; saveat=Ttr:Δt:T, seed=0)
# prob = ODEProblem(fitzhugh_nagumo_rule!,  u0, (0, T), p; saveat=Ttr:Δt:T, seed=0)
sol = solve(prob, SKenCarp(), saveat=0:Δt:T, maxiters=1e9);
# sol = solve(prob, Tsit5(), saveat=0:Δt:T, maxiters=1e9);
t = sol.t; tr = sol[:,:]'
fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, t, tr[:,1])
display(fig)
# xlims!(0,100)

# -------------------------------- Recurrences ------------------------------- #
include("$(scriptsdir())/utils.jl")
include("$(srcdir())/recurrences.jl")
T = 250; Ttr = 100; Δt = 0.01; a = -0.5; b = 0.1; c = 0.2; I = 0; n = 0 #nice tonic spiking
prob = ODEProblem(fitzhugh_nagumo_rule!,  u0, (0, T), p; )
sol = solve(prob, Tsit5(), maxiters=1e9, saveat=Ttr:Δt:T, seed=0);
tr = sol[:,:]'; t = sol.t;
fig = Figure(resolution = (1000,500))
axs=  plot_RM!(fig, t, Dataset(tr), 0.1; tsmode="lines", recurrencetimes=true,  logy=false, Δt)
fig
save("$(plotsdir())/fitzhughnagumo/timeseries-recurrence-a_$a-b_$b-c_$c-I_$I-n_$I.png", fig, px_per_unit=3)
# xlims!(axs[3], 0, 10)