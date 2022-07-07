using GLMakie, DynamicalSystems

# ----------------------------------------- Logistic map has attractor-merging crisis when attractors merge (duh) haha ---------------------------------------- #

# Cobweb interactive
# the second range is a convenience for intermittency example of logistic
using InteractiveDynamics 
n = 2000
Ttr = 2000
rrange = 3.4:0.0001:4.0; #interesting part of the diagram 
lo = Systems.logistic(0.4; r = rrange[1]);
# interactive_cobweb(lo, rrange, 5)
output = orbitdiagram(lo, 1, 1, rrange; n, Ttr);

L = length(rrange);
x = Vector{Float64}(undef, n*L);
y = copy(x);
for j in 1:L
    x[(1 + (j-1)*n):j*n] .= rrange[j]
    y[(1 + (j-1)*n):j*n] .= output[j]
end

fig, ax = scatter(x, y; axis = (xlabel = L"r", ylabel = L"x"),
    markersize = 0.8, color = ("black", 0.05),
)

##plot iterates of logistic map starting at 0.5
n = 20; Ttr=  0
output = orbitdiagram(lo, 1, 1, rrange; n, Ttr, u0=0.5);

L = length(rrange);
x = Vector{Float64}(undef, n*L);
y = copy(x);
for j in 1:L
    x[(1 + (j-1)*n):j*n] .= rrange[j]
    y[(1 + (j-1)*n):j*n] .= output[j]
end

scatter!(ax, x, y; markersize = 3, color = ("red", 0.3))

# ------------------------------------------------- Example on particle in double well (eg duffing). Following Ishii, 1986 (Physics Letters A) ------------------------------------------------ #
# """
# U(x) = ax^4/4 - bx 2/2 + cx (c: assymetry paramter); gamma: dissipation, f forcing omega forcing freq, n noise strength
# Alternativly, mdx2/dt2 = -νdx/dt - dV/dx + f*sen(ωt)
# """
using DynamicalSystemsBase:CDS
function duffing_assymetric(u0=[0.1, 0.1]; a=100, b=10, c=0, d=1, ω=3.5, f=0.849)
	return CDS(duffing_assymetric_rule, u0, [a,b,c,d,f,ω])
end

@inbounds function duffing_assymetric_rule(du, u, p, t)
    a,b,c,d,f, ω = p
    x, v = u
    du[1] = v
    du[2] = (-a*x^3 + b*x - c) - d*v + f*cos(ω*t)
end


# f = p = F = 
d = 1 # = ν/m = γ/m; m =1
b = 10 #-a in the paper 
a = 100 # = b in the paper 
c = 0 
ω = 3.5 # = Ω

T = 1000
Ttr = 0
Δt=0.1



fig = Figure(resolution=(1920, 1080)) 
axs = []
for (i,f) ∈ enumerate([0.849, 0.853])
	df = duffing_assymetric(; f)
	tr = trajectory(df, T; Ttr, Δt); t=Ttr:Δt:Ttr+T
	x = tr[:,1]
	ax = Axis(fig[i,1], ylabel="x", xlabel="t", title=["Before attractor-merging crisis (f = $(f))", "After attractor-merging crisis (f = $(f))"][i]); push!(axs, ax)
	lines!(t, x, color=x.>0, colormap=(:red, :black))
	ylims!(-0.6, 0.6)
end

f = 0.849
fig = Figure(resolution=(1920, 1080)) 
ax1 = Axis(fig[1,1], ylabel="x", xlabel="t", title=["Before attractor-merging crisis (f = $(f))", "After attractor-merging crisis (f = $(f))"][1]); push!(axs, ax)
df = duffing_assymetric(; f)
tr = trajectory(df, T, [0.1,0.1]; Ttr, Δt); t=Ttr:Δt:Ttr+T; x = tr[:,1];
# lines!(ax, t, x, color=[el > 0 ? :red : :black for el in x], label="initial condition on the right")
lines!(ax1, t, x, color=:black, label="initial condition on the right")
tr = trajectory(df, T, [-0.1,-0.1]; Ttr, Δt); t=Ttr:Δt:Ttr+T; x = tr[:,1];
lines!(ax1, t, x, color=:red, label="initial condition on the left")
# lines!(ax, t, x, color=[el > 0 ? :red : :black for el in x], label="initial condition on the left")
ylims!(-0.5, 0.5)
xlims!(Ttr, T)

f = 0.853
ax2 = Axis(fig[2,1], ylabel="x", xlabel="t", title=["Before attractor-merging crisis (f = $(f))", "After attractor-merging crisis (f = $(f))"][1]); push!(axs, ax)
df = duffing_assymetric(; f)
tr = trajectory(df, T, [0.1,0.1]; Ttr, Δt); t=Ttr:Δt:Ttr+T; x = tr[:,1];
# lines!(ax, t, x, color=[el > 0 ? :red : :black for el in x], label="initial condition on the right")
lines!(ax2, t, x, color=:black)
ylims!(-0.5, 0.5)
xlims!(Ttr, T)

Legend(fig[3,1], ax1, tellwidth = false, tellheight = true)
save("attractormergingcrisis-doublewellforcing.png", fig, px_per_unit=4)
