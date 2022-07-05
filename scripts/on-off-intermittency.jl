using GLMakie, DynamicalSystems

# ---------------------------------------------------------------- From Ashwin, 1999: on-off and in-out --------------------------------------------------------------- #
using DynamicalSystemsBase:DDS
"""
Ashwin, 1999
"""
function intermittency_ashwin(u0 = [0.25, 0.1]; r=3.82786, s = -0.3, ν=1.82)
    return DDS(intermittency_ashwin_rule, u0, [r, s, ν])
end
function intermittency_ashwin_rule(u, p, n)
	r, s, ν = p
	x, y = u
	dx1 = r*x*(1-x) + s*x*y^2 
	dx2 = ν*exp(-x)*y - y^3
	return SVector(dx1, dx2)
end

# r = 3.82786 #on-off, λₜ = 0.0024 (case 1 in Fig 1a)
r = 3.88615 #in-out, λₜ = 0.023  (case 2 in Fig 1a)

Ttr = 25000
T = 10000
ds = intermittency_ashwin(; r)
tr = trajectory(ds, T; Ttr); t=Ttr:Ttr+T

fig = Figure()
ax1 = Axis(fig[1,1], yscale=log10)
lines!(ax1, t, tr[:,2])
ax2 = Axis(fig[2,1])
scatter!(ax2, t, tr[:,1], markersize=4)
ax3 = Axis(fig[3,1])
lines!(ax3, t, tr[:,2])

fig = Figure()
for (i, r) ∈ enumerate([3.82786, 3.88615])
	ds = intermittency_ashwin(; r)
	T = [200000, 10000][i]; Ttr=[0, 25000][i] 
	tr = trajectory(ds, T; Ttr); t=Ttr:Ttr+T
	ax1 = Axis(fig[1,i], yscale=log10, title=["on-off intermittency", "in-out intermittency"][i], ylabel="y")
	lines!(ax1, t, tr[:,2], color=:black)
	ax2 = Axis(fig[2,i], ylabel="x")
	scatter!(ax2, t, tr[:,1], markersize=1.5, color=:black)
	ax3 = Axis(fig[3,i], ylabel="y", xlabel="x")
	lines!(ax3, t, tr[:,2], color=:black)
end
save("onoff-and-inout-intermittency-ashwin1999.png", fig)


# ----------------------------------------------------------------- Coupled rossler: on-off ---------------------------------------------------------------- #
# http://www.scholarpedia.org/article/Bubbling_transition
function coupled_roessler_scholar(u0=[1, -2, 0, 0.11, 0.2, 0.1];
    K = 0.045, a=0.2, b=0.2, c=6.6)
    p = [K, a, b, c]
    return ContinuousDynamicalSystem(coupled_roessler_scholar_rule, u0, p)
end
function coupled_roessler_scholar_rule(u,p,t)
	K, a, b, c = p
	x1, y1, z1, x2, y2, z2 = u
	dx1 = -y1 - z1 + K*(x2-x1)
	dy1 = x1 + a*y1 +K*(y2 - y1)
	dz1 = b + z1*(x1 - c) + K*(z2 - z1)
	dx2 = -y2 - z2 + K*(x1-x2)
	dy2 = x2 + a*y2 +K*(y1 - y2)
	dz2 = b + z2*(x2 - c) + K*(z1 - z2)
    return SVector(dx1, dy1, dz1, dx2, dy2, dz2)
end

Ttr = 0
T = 10000
Δt = 0.01

# K = 0.045
# K = 0.02
K = 0.03

cr = coupled_roessler_scholar(; K)
tr = trajectory(cr, T; Ttr, Δt); t = Ttr:Δt:T 

x1 = tr[:,1];
x2 = tr[:,4];
δx = x2 .- x1;
lines(t, δx, color=δx .< 0.3)


fig = Figure()
ax1 = Axis(fig[1,1], ylabel="x" )
ax2 = Axis(fig[2,1], ylabel="y")
ax3 = Axis(fig[3,1], ylabel="z", xlabel="t")
for i =1:2
	x = tr[:, (i-1)*3 + 1]
	y = tr[:, (i-1)*3 + 2]
	z = tr[:, (i-1)*3 + 1]
	lines!(ax1, t, x)
	lines!(ax2, t, y)
	lines!(ax3, t, z)
end
δx = tr[:,4] .- tr[:,1];
δy = tr[:,5] .- tr[:,2];
δz = tr[:,6] .- tr[:,3];
ax4 = Axis(fig[1,2], ylabel="δx" )
lines!(ax4, t, δx, color=:black)
ax5 = Axis(fig[2,2], ylabel="δy")
lines!(ax5, t, δy, markersize=1.5, color=:black)
ax6 = Axis(fig[3,2], ylabel="δz", xlabel="t")
lines!(ax6, t, δz, color=:black)

ax7 = Axis3(fig[4:5,1])
lines!(ax7, tr[:,1], tr[:,2], tr[:,3]) 
ax8 = Axis3(fig[4:5,2])
lines!(ax8, tr[:,4], tr[:,5], tr[:,6]) 

ax9 = Axis(fig[6:7,1], xlabel="x_1", ylabel="x_2")
lines!(ax9, tr[:,1], tr[:,4], markersize=1, color=(:blue, 0.1))
lines!(ax9, minimum(tr[:,1]):0.01:maximum(tr[:,1]), minimum(tr[:,1]):0.01:maximum(tr[:,1]), color=:red, linestyle="--")
ax10 = Axis(fig[6:7,2], xlabel="y_1", ylabel="y_2")
lines!(ax10, tr[:,2], tr[:,5], markersize=1, color=(:blue, 0.1))
lines!(ax10, minimum(tr[:,2]):0.01:maximum(tr[:,2]), minimum(tr[:,2]):0.01:maximum(tr[:,2]), color=:red, linestyle="--")
# ax11 = Axis3(fig[6,3], xlabel="z_1", ylabel="z_2")
# scatter!(ax11, tr[:,3], tr[:,6], markersize=1, color=(:cyan, 0.7))
save("onoff-intermittency-coupledrosslers-K_$(K).png", fig, px_per_unit=4)
