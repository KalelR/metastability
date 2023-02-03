using DrWatson
@quickactivate "metastability"
using DynamicalSystems, Statistics
# using GLMakie
using CairoMakie
include("$(srcdir())/paperplottheme.jl")
include("$(scriptsdir())/utils.jl")


# ds = Systems.logistic(r = 3.8282)
# tr = trajectory(ds, 150)
# t = 1:151
# 
# fig = Figure(; resolution) 
# ax = Axis(fig[1,1])
# lines!(ax, t, tr, color=:black)
# xmin = 50; xmax = 140
# xlims!(ax, xmin, xmax)
# xticks = range(xmin, xmax, step=20)
# ax.xticks = (xticks, string.( xticks .- xmin) )
# # ax.xticklabels = string.(0:20:60)
# hidespines!(ax, :r, :t, :b, :l)
# hideydecorations!(ax)
# hidedecorations!(ax)
# # gca().set_xticklabels(string.(0:20:100))
# # savefig("logistic-intermittency.png", dpi=300)
# save("$(plotsdir())/paper/illustration-definition-activitypattern.pdf", fig)


width = 190 #pt 
height = 80 #pt
resolution = 1 .*(width, height)
include("$(srcdir())/systems/logistic.jl")
function plot_varactivitypattern_definition()
T = 15000
Ttr = 500
r = rc-1e-7
lo = Systems.logistic(0.4; r); t = Ttr:Ttr+T
traj = trajectory(lo, T; Ttr)
colors, traj, t = logistic_get_colors_trajectory!(traj, t);

fig = Figure(; resolution)
# ax = Axis(fig[1,1])
ax = Axis(fig[1,1]; alignmode=Outside(-5))
scatterlines!(ax, t, traj, color=colors, markersize=2, linewidth=0.8)
ylims!(ax, 0, 1)
xlims!(ax, 10700, 10950)
hidespines!(ax, :r, :t, :b, :l)
hidedecorations!(ax)
save("$(plotsdir())/paper/illustration-definition-activitypattern.pdf", fig)
end


function plot_varsync_definition()
fig = Figure(; resolution) 
ax = Axis(fig[1,1])
t = 0:0.01:16pi
y = sin.(t)

y2 = deepcopy(y) 
half = floor(Int64, length(t)/2)
y2[1:half] = 1.25 .* y[1:half]
y2[half+1:end] .= 1.25 .* sin.(2 .* t[half+1:end])

y  = y .+ 0.1 .* rand(length(y))
y2  = y2 .+ 0.1 .* rand(length(y2))

lines!(ax, t,y, color=:black)
# lines!(ax, t,y2,color=:purple, linestyle=:dash)
lines!(ax, t,y2,color=:purple)
xlims!(ax, 5,42)

hidespines!(ax, :r, :t, :b, :l)
hidedecorations!(ax)
# ylabel("Phases")
# xlabel("t")
# gcf().tight_layout()
save("$(plotsdir())/paper/illustration-definition-varsyn.pdf", fig)
end