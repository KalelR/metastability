using DrWatson
@quickactivate "metastability"
using GLMakie
include("$(srcdir())/paperplottheme.jl")
using DynamicalSystems, Statistics
# include("$(scriptsdir())/utils.jl")

sols = []

# Noisy Duffing 
prob = noisy_duffing()
Ttr = 0; Δt = 0.5
sol = solve(prob_duffing, SKenCarp(), saveat=Ttr:Δt:T, maxiters=1e9, progress=true);
push!(sols, sol)


# Heteroclinic cycle
prob = heteroclinic_cycle()
Ttr = 1000; Δt = 1.0;
sol = solve(hcgh, Vern9(), maxiters=1e9, saveat=Ttr:Δt:T); ts = sol.t;
push!(sols, sol)

