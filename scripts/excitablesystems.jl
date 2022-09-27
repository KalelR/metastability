using DrWatson
@quickactivate

using DynamicalSystems
u0 = [0.5, 0.5]
a = 3.0
b = 0.2
ε,= 0.01
I = 0.0
ds = Systems.fitzhugh_nagumo(u0; a, b, ε, I)

tr = trajectory(ds, 1000); t=0:1000

using CairoMakie
fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, t, tr[:,1])
display(fig)