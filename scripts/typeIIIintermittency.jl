using GLMakie, DynamicalSystems


using DynamicalSystemsBase:DDS
"""
intermittency following subcritical period doubling at μ=0.0 (so for μ>0) (See Argyris)
"""
function typeIIIintermittency(u0 = 0.25; μ=0.01)
    return DDS(typeIIIintermittency_rule, u0, [μ])
end
function typeIIIintermittency_rule(x, p, n)
	@inbounds mod( 1 - 2x - 1/(2π) * (1 - p[1]) * cos(2π*(x- (1/12))) , 1)
end

Ttr = 0 
T = 1000
μ = 0.01
ds = typeIIIintermittency(rand(); μ)

tr = trajectory(ds, T); t = Ttr:T

fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, t, tr)


fig = Figure()
for (i, μ) ∈ enumerate([0.01, 0.015, 0.02])
ds = typeIIIintermittency(0.1; μ)
ax = Axis(fig[i,1], title="μ = $(μ)")
tr = trajectory(ds, T, 1/3+1e-4); t = Ttr:T
tr2 = trajectory(ds, T, 0.4); t = Ttr:T
lines!(ax, t, tr, color=:black, label="ic near FP at x=1/3")
lines!(ax, t, tr2, color=(:orange, 0.8), label = "ic far from FP")
end
fig[4,1] = Legend(fig, ax, tellwidth=false)
save("typeIIIintermittency-trajectories.png", fig)