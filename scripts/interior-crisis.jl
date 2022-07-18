using GLMakie, DynamicalSystems

# ----------------------------------------- Logistic map has interior crisis at the end of periodic windows ---------------------------------------- #

# Cobweb interactive
# the second range is a convenience for intermittency example of logistic
using InteractiveDynamics 
n = 2000
Ttr = 2000
rrange = 3.4:0.0005:4.0; #interesting part of the diagram 
rrange = 3.82:0.00001:3.88; #period 3
rrange = 3.855:0.000001:3.859; #to find the boundary crisis point; critical value around 3.8568; can pick before 3.8567 and after at 3.8569
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

T = 2000
fig = Figure(resolution=(1920, 1080)) 
axs = []
for (i, r) ∈ enumerate([3.8567, 3.8569])
	lo = Systems.logistic(0.4; r);
	tr = trajectory(lo, T; Ttr)
	x = tr[:,1];  t = Ttr:Ttr+T
	ax = Axis(fig[1,i], ylabel="x_n", xlabel="n", title=["Before interior crisis (r = $(r))", "After interior crisis (r = $(r))"][i]); push!(axs, ax)
	# scatterlines!(ax, t, x, color=(:black, 1.0), markersize=2)
	scatter!(ax, t, x, color=(:black, 1.0), markersize=2)
end
save("interiorcrisis-logisticmap.png", fig, px_per_unit=4)


# ---------------------------------------------------- Ikeda map also famous for interior crisis --------------------------------------------------- #
# translation from ds to my own notation: t = χ, c=C, d=K (= a in Alligood), a=A, b=B 
a = 0.84
b = 0.9 
c = 0.4 
# d = -7.24
d = 7.1
ik = Systems.ikedamap(;a, b, c, d)

T = 1e6
Ttr = 1e5
# tr = trajectory(ik, T; Ttr)
# x = tr[:,1]; y = tr[:,2]

# fig = Figure() 
# ax = Axis(fig[1,1])
# scatter!(ax, x, y, markersize=2, color=:black)


using DataStructures
function histmeasure(tr, numdigits)
	v = [[tr[i,1], tr[i,2]] for i=1:length(tr)]
	vround = [round.(el, digits=numdigits) for el in v]
	c = counter(vround)
	measure = [c[elround] for elround in vround ]
	measure ./ length(v)
end


fig = Figure(resolution=(1920, 1080)) 
axs = []
for (i, d) ∈ enumerate([7.1, 7.3])
	ik = Systems.ikedamap(;a, b, c, d)
	tr = trajectory(ik, T; Ttr)
	x = tr[:,1]; y = tr[:,2]; t = Ttr:Ttr+T
	measure = histmeasure(tr, 2)
	ax = Axis(fig[1:3,i], ylabel="y", xlabel="x", title=["Before interior crisis (d = $(d))", "After interior crisis (d = $(d))"][i]); push!(axs, ax)
	scatter!(ax, x, y, markersize=3, color=measure)
	xlims!(-0.6, 1.7)
	ylims!(-1.8, 1.4)
	ax = Axis(fig[4, i], ylabel="x")
	lines!(t, x, color=:black)
	xlims!(T-800, T)
	ylims!(-0.6, 1.7)
	ax = Axis(fig[5, i], ylabel="y", xlabel="t")
	lines!(t, y, color=:black)
	xlims!(T-800, T)
	ylims!(-1.8, 1.4)
end
save("interiorcrisis-ikedamap.png", fig, px_per_unit=4)





##--------------------------------------------------animations --------------------a = 0.84
include("utils.jl")
plotsdir() = "../plots/"
a = 0.84
b = 0.9 
c = 0.4 
# d = 7.1
d = 7.3
ik = Systems.ikedamap(;a, b, c, d)

# T = 10e3
T = 5e3
Ttr = 1e3
tr = trajectory(ik, T; Ttr)

x = tr[:,1]; y = tr[:,2]; t = Ttr:Ttr+T
# sol = [ [tr[i,1], tr[i,2]] for i=1:length(tr)]
p = [a, b, c, d]
speeds = norm.([ik.jacobian([tr[i,1], tr[i,2]], p, 1) for i=1:length(tr)])

Δt = 1.0
Tplot= T-Ttr
framerate=200
# framerate=400
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


points = Observable(Point2f[(tr_plot[1,1], tr_plot[1,2])])
colors_ob = Observable([colors[1]])
time = Observable(t_plot[1])
fig = Figure(resolution=(800, 600))
ax = Axis(fig[1,1], title= @lift("t = $($time)"))
s=scatter!(ax, points, color=colors_ob, markersize=3)
hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
xlims!(-0.6, 1.7)
ylims!(-1.8, 1.4)
record(fig, "../plots/interiorcrisis/ikedamap-animation-a_$(a)-b_$(b)-c_$(c)-d_$(d).mp4", frames;
        framerate) do frame
    time[] = t_plot[frame]
    new_point = Point2f(tr_plot[frame,1], tr_plot[frame,2])
	points[] = push!(points[], new_point)
	colors_ob[] = push!(colors_ob[], colors[frame])
end

















ax = Axis(fig[1,1], ylabel="y", xlabel="x"); 
scatter!(ax, x, y, markersize=3, color=measure)
