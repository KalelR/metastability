using GLMakie, DynamicalSystems

# ----------------------------------------- Logistic map has boundary crisis at r = 4 ---------------------------------------- #

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
for (i, r) ∈ enumerate([4.0, 4.00001])
	lo = Systems.logistic(0.4; r);
	tr = trajectory(lo, T; Ttr)
	x = tr[:,1];  t = Ttr:Ttr+T
	ax = Axis(fig[1,i], ylabel="x_n", xlabel="n", title=["Before interior crisis (r = $(r))", "After interior crisis (r = $(r))"][i]); push!(axs, ax)
	# scatterlines!(ax, t, x, color=(:black, 1.0), markersize=2)
	scatter!(ax, t, x, color=(:black, 1.0), markersize=2)
	ylims!(-6, 1)
end
save("boundarycrisis-logisticmap.png", fig, px_per_unit=4)


# ---------------------------------------------------- Ikeda map also famous for boundary crisis --------------------------------------------------- #
# translation from ds to my own notation: t = χ, c=C, d=K (= a in Alligood), a=A(=r in Alligood), b=B 
b = 0.9 
c = 0.4 
d = 6.0
a = 1.0 # a=A=r varying
a = 1.003
# u0 = [0.875, -1.1]
u0 = [-2.0, -2.0]
ik = Systems.ikedamap(u0;a, b, c, d)

T = 1e6
Ttr = 0
tr = trajectory(ik, T; Ttr); t = Ttr:Ttr+T
x = tr[:,1]; y = tr[:,2]
measure = histmeasure(tr, 2)

fig = Figure() 
ax = Axis(fig[1:3,1])
scatter!(ax, x, y, markersize=2, color=measure)

ax = Axis(fig[4, 1], ylabel="x")
lines!(t, x, color=:black)
save("$(plotsdir())/boundarycrisis/boundarycrisis-ikedamap-statespace-timeseries-a_$(a).png", fig, px_per_unit=4)



fig = Figure(resolution=(1920, 1080)) 
axs = []
for (i, a) ∈ enumerate([0.997, 1.003])
	ik = Systems.ikedamap(;a, b, c, d)
	tr = trajectory(ik, T, u0; Ttr)
	x = tr[:,1]; y = tr[:,2]; t = Ttr:Ttr+T
	measure = histmeasure(tr, 2)
	ax = Axis(fig[1:3,i], ylabel="y", xlabel="x", title=["Before boundary crisis (a = $(a))", "After boundary crisis (a = $(a))"][i]); push!(axs, ax)
	scatter!(ax, x, y, markersize=3, color=measure)
	xlims!(-0.6, 6.2)
	ylims!(-2.8, 6)
	ax = Axis(fig[4, i], ylabel="x")
	lines!(t, x, color=:black)
	xlims!(0, 20000)
	ylims!(-0.6, 6.2)
	ax = Axis(fig[5, i], ylabel="y", xlabel="t")
	lines!(t, y, color=:black)
	xlims!(0, 20000)
	ylims!(-0.6, 6.2)
end
save("boundarycrisis-ikedamap.png", fig, px_per_unit=4)
#just time-series
for (i, a) ∈ enumerate([0.997, 1.003])
	ik = Systems.ikedamap(;a, b, c, d)
	tr = trajectory(ik, T, u0; Ttr)
	x = tr[:,1]; y = tr[:,2]; t = Ttr:Ttr+T
	# ax = Axis(fig[1:3,i], ylabel="y", xlabel="x", title=["Before boundary crisis (a = $(a))", "After boundary crisis (a = $(a))"][i]); push!(axs, ax)
	fig = Figure(resolution=(800, 250))
	ax = Axis(fig[1, 1], ylabel="x", xlabel="t", xautolimitmargin=(1.0, 1.0))
	lines!(t, x, color=:black)
	xlims!(0, 20500)
	ylims!(-0.6, 6.2)
	ax.xticks = ([0,10000, 20000], ["0", "10000", "20000"])
	# ax = Axis(fig[2, 1], ylabel="y", xlabel="t")
	# lines!(t, y, color=:black)
	# xlims!(0, 20000)
	# ylims!(-0.6, 6.2)
	save("../plots/boundarycrisis/boundarycrisis-ikedamap-timeseries-a_$(a).png", fig, px_per_unit=4)

end


##distribution of times in saddle 




#------------------------------------------- animation ----------
include("utils.jl")
plotsdir() = "../plots/"
a = 0.997
# a = 1.003 #after crisis, CS
T = 16400
Ttr = 15500
ik = Systems.ikedamap(;a, b, c, d)
tr = trajectory(ik, T, u0; Ttr)
x = tr[:,1]; y = tr[:,2]; t = Ttr:Ttr+T
sol = [ [tr[i,1], tr[i,2]] for i=1:length(tr)]
p = [a, b, c, d]
speeds = norm.([ik.jacobian([tr[i,1], tr[i,2]], p, 1) for i=1:length(tr)])

Δt = 1.0
Tplot= T-Ttr
framerate=50
# framerate=50
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
xlims!(-0.6, 6.2)
ylims!(-2.8, 6)
record(fig, "../plots/boundarycrisis/ikedmap-animation-a_$(a)-b_$(b)-c_$(c)-d_$(d).mp4", frames;
        framerate) do frame
    time[] = t_plot[frame]
    new_point = Point2f(tr_plot[frame,1], tr_plot[frame,2])
	points[] = push!(points[], new_point)
	colors_ob[] = push!(colors_ob[], colors[frame])
end



