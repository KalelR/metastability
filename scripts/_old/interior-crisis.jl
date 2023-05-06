using DrWatson
@quickactivate "metastability"
using GLMakie, DynamicalSystems
dirname = "interiorcrisis"
include("$(scriptsdir())/utils.jl")

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
save("$(plotsdir())/$(dirname)/interiorcrisis-logisticmap.png", fig, px_per_unit=4)


# ---------------------------------------------------- Ikeda map also famous for interior crisis --------------------------------------------------- #
# translation from ds to my own notation: t = χ, c=C, d=K (= a in Alligood), a=A, b=B 
a = 0.84
b = 0.9 
c = 0.4 
# d = -7.24
d = 7.1
# ik = Systems.ikedamap(;a, b, c, d)

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

function trajectoryspeed(ik, tr, t)
	Δxs = diff(Matrix(tr), dims=1)
	speeds = mapslices(x->norm(x), Δxs, dims=2)[:,1]
	tr_cut = tr[1:end-1, :]; t_cut = t[1:end-1]
	return speeds, tr_cut, t_cut
end

fps_5 = [-0.149036 -0.111942; 0.6901164331310736 -0.07534279207228838; 0.7081550024831971 0.6107258394579234; 0.08793801081244812 -0.377787118392552; 1.0755007518118098 -0.25769927892208205]
fig = Figure(resolution=(1920, 1080)) 
axs = []
for (i, d) ∈ enumerate([7.1, 7.4])
	ik = Systems.ikedamap(;a, b, c, d)
	tr = trajectory(ik, T; Ttr); t = Ttr:Ttr+T
	# measure = histmeasure(tr, 2); measure_label = "ρ(x)"
	measure, tr, t = trajectoryspeed(ik, tr, t); measure_label = "speed(x)"
	x = tr[:,1]; y = tr[:,2]; 
	ax = Axis(fig[1:3,i][1,1], ylabel="y", xlabel="x", title=["Before interior crisis (d = $(d))", "After interior crisis (d = $(d))"][i]); push!(axs, ax)
	s=scatter!(ax, x, y, markersize=1.5, color=measure)
	Colorbar(fig[1:3,i][1,2], s, label=measure_label)
	scatter!(ax, fps_5[:,1], fps_5[:,2], color=:red, markersize=4)
	xlims!(-0.6, 1.7)
	ylims!(-1.8, 1.4)
	ax = Axis(fig[4, i], ylabel="x")
	lines!(t, x, color=measure)
	for j=1:5 hlines!(ax, fps_5[j,1], color=(:red, 0.5), linestyle="--") end
	xlims!(T-1500, T)
	ylims!(-0.6, 1.7)
	ax = Axis(fig[5, i], ylabel="y", xlabel="t")
	lines!(t, y, color=measure)
	for j=1:5 hlines!(ax, fps_5[j,2], color=(:red, 0.5), linestyle="--") end
	xlims!(T-1500, T)
	ylims!(-1.8, 1.4)
end
save("$(plotsdir())/$(dirname)/interiorcrisis-ikedamap-colorasspeed.png", fig, px_per_unit=4)
# save("$(plotsdir())/$(dirname)/interiorcrisis-ikedamap-colorasdensity.png", fig, px_per_unit=4)


# color red for low measure, blue for high measure : does not work well
# fig = Figure(resolution=(1920, 1080)) 
# axs = []
# ms_threshold = 
# for (i, d) ∈ enumerate([7.22, 7.24])
# 	ik = Systems.ikedamap(;a, b, c, d)
# 	tr = trajectory(ik, T; Ttr)
# 	x = tr[:,1]; y = tr[:,2]; t = Ttr:Ttr+T
# 	measure = histmeasure(tr, 1)
# 	if i == 1 ms_threshold = minimum(measure) end
# 	ax = Axis(fig[1:3,i], ylabel="y", xlabel="x", title=["Before interior crisis (d = $(d))", "After interior crisis (d = $(d))"][i]); push!(axs, ax)
# 	scatter!(ax, x, y, markersize=3, color=[el > ms_threshold ? :blue : :red for el in measure])
# 	xlims!(-0.6, 1.7)
# 	ylims!(-1.8, 1.4)
# 	ax = Axis(fig[4, i], ylabel="x")
# 	lines!(t, x, color=:black)
# 	xlims!(T-800, T)
# 	ylims!(-0.6, 1.7)
# 	ax = Axis(fig[5, i], ylabel="y", xlabel="t")
# 	lines!(t, y, color=:black)
# 	xlims!(T-800, T)
# 	ylims!(-1.8, 1.4)
# end
# save("$(plotsdir())/$(dirname)/interiorcrisis-ikedamap-coloroutsidecoreaslowmeasure.png", fig, px_per_unit=4)

#outside core based on box containing core before crisis
# fig = Figure(resolution=(1920, 1080)) 
# axs = []
# ms_threshold = 
# for (i, d) ∈ enumerate([7.1, 7.3])
# 	ik = Systems.ikedamap(;a, b, c, d)
# 	tr = trajectory(ik, T; Ttr)
# 	x = tr[:,1]; y = tr[:,2]; t = Ttr:Ttr+T
# 	measure = histmeasure(tr, 2)
# 	if i == 1 box = [minimum(tr[:,1]), maximum()] end
# 	ax = Axis(fig[1:3,i], ylabel="y", xlabel="x", title=["Before interior crisis (d = $(d))", "After interior crisis (d = $(d))"][i]); push!(axs, ax)
# 	scatter!(ax, x, y, markersize=3, color=[el > ms_threshold ? :blue : :red for el in measure])
# 	xlims!(-0.6, 1.7)
# 	ylims!(-1.8, 1.4)
# 	ax = Axis(fig[4, i], ylabel="x")
# 	lines!(t, x, color=:black)
# 	xlims!(T-800, T)
# 	ylims!(-0.6, 1.7)
# 	ax = Axis(fig[5, i], ylabel="y", xlabel="t")
# 	lines!(t, y, color=:black)
# 	xlims!(T-800, T)
# 	ylims!(-1.8, 1.4)
# end

# d=7.22
# ik = Systems.ikedamap(;a, b, c, d)
# tr = trajectory(ik, T; Ttr)
# x = tr[:,1]; y = tr[:,2]; t = Ttr:Ttr+T
# using GeometricalPredicates

# ll = [minimum(tr[:,1]), minimum(tr[:,2])]
# lr = [maximum(tr[:,1]), minimum(tr[:,2])]
# ul = [minimum(tr[:,1]), maximum(tr[:,2])]
# ur = [maximum(tr[:,1]), maximum(tr[:,2])]
# points = [ll, ul, ur, lr]
# poly = Polygon([GeometricalPredicates.Point(pt...) for pt in points]...)
# m_points = reduce(hcat, push!(points, points[1]))'

# fig = Figure() 
# ax = Axis(fig[1,1])
# scatter!(ax, x, y, markersize=3)
# scatterlines!(ax, m_points[:,1], m_points[:,2], color=:red, markersize=10)

# d = 7.24
# ik = Systems.ikedamap(;a, b, c, d)
# tr = trajectory(ik, T; Ttr)
# x = tr[:,1]; y = tr[:,2]; t = Ttr:Ttr+T
# fig = Figure() 
# ax = Axis(fig[1,1])
# scatter!(ax, x, y, markersize=3)
# scatterlines!(ax, m_points[:,1], m_points[:,2], color=:red, markersize=10)

# inpolygon(poly, GeometricalPredicates.Point(1.1,1.1)) 

# ---------------------------- Lyapunov exponents ---------------------------- #
ds = 6.8:0.1:7.4
all_λs = []
T=5e7
Ttr=2e6
for d ∈ ds
	ik = Systems.ikedamap(;a, b, c, d)
	# λs = lyapunovspectrum(ik, Int(T); Ttr=Int(Ttr))
	λ = lyapunov(ik, T; Ttr)
	push!(all_λs, λ)
end
all_λs = reduce(hcat, all_λs)'
	
fig = Figure()
ax = Axis(fig[1,1], xlabel="d", ylabel="Lyapunov exponents")
scatterlines!(ax, ds, all_λs[:,1])
# scatterlines!(ax, ds, all_λs[:,2])
save("$(plotsdir())/$(dirname)/lyapunovexponents-interiorcrisis-ikedamap.png", fig, px_per_unit=4)

# ------------------------------ Recurrence plot ----------------------------- #
a = 0.84
b = 0.9 
c = 0.4 
d = 7.35
T = 1000; Ttr = 2e3-T
# T=20
# T = 800; Ttr = 0
ϵ = 0.1
u0 = [1.0, 1.0] #default
u0 = [-3.0, -3.0]
ik = Systems.ikedamap(;a, b, c, d)
tr = trajectory(ik, T, u0; Ttr); ts = Ttr:Ttr+T
fig = Figure(resolution=(1500, 500), fontsize=30, figure_padding=(5, 30, 5, 30))
axs = plot_RM!(fig, ts, tr, ϵ; tsmode="line")
RM = RecurrenceMatrix(tr, ϵ)
mrt = meanrecurrencetime(RM)
axs[2].title = "mean recurrence time = $(round(mrt, digits=1))"

# rqa(RM)
rs = recurrencestructures(RM, recurrencetimes=true)
recurrencetimes = rs["recurrencetimes"]
ax = Axis(fig[1,3], yscale=log10, ylabel="PDF(recurrence times)", xlabel="recurrence times")
numbins = 10; weights, bins = histogram(recurrencetimes, numbins); 
weights .+= 1e-7; ylims!(ax, 1e-6, maximum(weights))
scatterlines!(ax, bins[1:end-1], weights, color=:black)
save("$(plotsdir())/$(dirname)/interiorcrisis-ikedamap-recurrenceplot-recurrenceth_$(ϵ)-u0_$(u0)-Ttr_$(Ttr).png", fig, px_per_unit=4)
# ---------------------------- frequency analysis ---------------------------- #
include("$(scriptsdir())/utils.jl")
using FFTW

a = 0.84
b = 0.9 
c = 0.4 
d = 7.3
# d = 6.8 #periodic
# d = 5.589
T = 50000; Ttr = 1e6-800
T=50; Ttr=1e3 #periodic, d=6.8
T=1000; Ttr=1e3 #periodic, d=6.8
# T = 20
# T = 800; Ttr = 0
u0 = [1.0, 1.0] #default
u0 = rand(2)
ik = Systems.ikedamap(;a, b, c, d)
tr = trajectory(ik, T, u0; Ttr); ts = Ttr:Ttr+T


# signal = tr[:,1]
Ts = Δt = 1


using GLMakie
axs = []
fig = Figure()
for i=1:2
	signal = tr[:,i]
	variable = ["x", "y"][i]
	F, freqs, periods, maxperiod, Fmax = getspectrum2(signal, ts, Ts)
	ax1 = Axis(fig[1,i], title="Tmax = $maxperiod", ylabel="$variable", xlabel="t")
	time_domain = lines!(ax1, ts, signal, color=:black)
	# F .+= 1
	# ax2 = Axis(fig[2,i], ylabel="FFT($variable)", xlabel="Frequency", yscale=log10)
	ax2 = Axis(fig[2,i], ylabel="FFT($variable)", xlabel="Frequency")
	freq_domain = lines!(ax2, freqs, abs.(F), color=:black) 
	scatter!(ax2, [1/maxperiod], [Fmax], color=:red)
	push!(axs, [ax1, ax2]...)
end
linkxaxes!(axs[1], axs[3])
linkaxes!(axs[2], axs[4])
# ylims!(axs[2], 1e-4, 1e6)

save("$(plotsdir())/$(dirname)/interiorcrisis-spectrumanalysis-d_$(d)-all.png", fig)



##--------------------------------------------------animations --------------------a = 0.84
# include("utils.jl")
# plotsdir() = "../plots/"
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
points2 = Observable(Point2f[(t_plot[1], tr_plot[1,2])])
colors_ob = Observable([colors[1]])
tanim = Observable(t_plot[1])
fig = Figure(resolution=(700, 600))
ax = Axis(fig[1:2,1], title= @lift("t = $($tanim)"), ylabel="y", xlabel="x")
s=scatter!(ax, points, color=colors_ob, markersize=3)
ax2 = Axis(fig[3,1], ylabel="x", xlabel="t")
lines!(ax2, points2)
hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
xmin = -0.6; xmax=1.7
xlims!(ax, xmin, xmax); ylims!(ax, -1.8, 1.4)
xlims!(ax2, t_plot[1], t_plot[end]); ylims!(ax2, xmin, xmax)
record(fig, "$(plotsdir())/interiorcrisis/ikedamap-animation-withtimeseries-a_$(a)-b_$(b)-c_$(c)-d_$(d).mp4", frames;
        framerate) do frame
    tanim[] = t_plot[frame]
    new_point = Point2f(tr_plot[frame,1], tr_plot[frame,2])
    new_point2 = Point2f(t_plot[frame,1], tr_plot[frame,1])
	points[] = push!(points[], new_point)
	points2[] = push!(points2[], new_point2)
	colors_ob[] = push!(colors_ob[], colors[frame])
end



# ------------------------------ color outside edge  ----------------------------- #
T = 15e3; Ttr = 1e3
T = 3e3; Ttr = 1e3
a = 0.84; b = 0.9; c = 0.4 
d = 7.22
ik = Systems.ikedamap(;a, b, c, d)
traj = trajectory(ik, T; Ttr); t = Ttr:Ttr+T
innercore = deepcopy(Matrix(traj))
fps_5 = [-0.149036 -0.111942; 0.6901164331310736 -0.07534279207228838; 0.7081550024831971 0.6107258394579234; 0.08793801081244812 -0.377787118392552; 1.0755007518118098 -0.25769927892208205]

# measure = histmeasure(tr, 1)
# inner_threshold = minimum(measure)*10
maxdist_threshold = 0.1
maxdist_threshold = 0.04

# d = 7.22
# d = 7.248
d = 7.3
c1 = d <= 7.22 ? :black : :green
ik = Systems.ikedamap(;a, b, c, d)
traj = trajectory(ik, T; Ttr); t = Ttr:Ttr+T
# measure = histmeasure(tr, 2)
arewithin = mapslices(x->withinset(x, innercore, maxdist_threshold), Matrix(traj), dims=2)[:,1]; #a bit expensive, have to compare each point to all other points in the set

Δt = 1.0
Tplot= 1000
framerate=5

t_plot, tr_plot, frames = animationdata_tr(traj, t, Tplot, Δt, Δt)
arewithin_plot = mapslices(x->withinset(x, innercore, maxdist_threshold), Matrix(tr_plot), dims=2)[:,1]; #a bit expensive, have to compare each point to all other points in the set
# measure_plot = histmeasure(tr_plot, 1)
# colors = pointspeed_as_colors(speeds_plot);

#observables
points = Observable(Point2f[(tr_plot[1,1], tr_plot[1,2])])
points2 = Observable(Point2f[(t_plot[1], tr_plot[1,1])])
# colors_ob = Observable([colors[1]])
tanim = Observable(t_plot[1])

#plot
fig = Figure(resolution=(700, 600))
ax = Axis(fig[1:3,1][1,1], title= @lift("t = $($tanim)"), ylabel="y", xlabel="x")
s=scatter!(ax, traj[:,1], traj[:,2,], color=[el == 1  ? (c1, 0.6) : (:purple, 0.6) for el in arewithin], markersize=3)
s=scatter!(ax, points, color=:orange, markersize=12)
scatter!(ax, fps_5[:,1], fps_5[:,2], color=:red, markersize=10)

ax2 = Axis(fig[4,1], ylabel="x", xlabel="t")
lines!(ax2, t_plot, tr_plot[:,1], color=[el == 1  ? (c1,0.8) : (:purple,0.8) for el in arewithin_plot], markersize=5)
scatter!(ax2, points2, color=:orange, markersize=12)
for j=1:5 hlines!(ax2, fps_5[j,1], color=(:red, 0.5), linestyle="--") end
#decoration
hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
xmin = -0.6; xmax=1.7
xlims!(ax, xmin, xmax); ylims!(ax, -1.8, 1.4)
xlims!(ax2, t_plot[1], t_plot[end]); ylims!(ax2, xmin, xmax)
record(fig, "$(plotsdir())/interiorcrisis/ikedamap-animation-withtimeseries-coloredbasedondistance-a_$(a)-b_$(b)-c_$(c)-d_$(d)-framerate_$(framerate).mp4", frames;
        framerate) do frame
    tanim[] = t_plot[frame]
    new_point = Point2f(tr_plot[frame,1], tr_plot[frame,2])
    new_point2 = Point2f(t_plot[frame,1], tr_plot[frame,1])
	points[] = [new_point]
	points2[] = [new_point2]
	# points2[] = push!(points2[], new_point2)
	# colors_ob[] = push!(colors_ob[], colors[frame])
end







# ------------------------ finding the period-5 saddle ----------------------- #
a = 0.84
b = 0.9 
c = 0.4 
d = 7.1
# d = 7.3
u0 = rand(2)
p = [a, b, c, d]

fp = periodicorbits(ik, 5, [rand(2) for i=1:100])
scatter!(fp[:,1], fp[:,2], color=:red)

@inbounds function ikedamap_rule(u, p, n)
    a,b,c,d  = p
    t = c - d/(1 + u[1]^2 + u[2]^2)
    dx = a + b*(u[1]*cos(t) - u[2]*sin(t))
    dy = b*( u[1]*sin(t) + u[2]*cos(t) )
    return SVector{2}(dx, dy)
end

# @inbounds function ikedamap_5thorder_rule(u,p,n)
# 	u_mapped = deepcopy(u)
# 	for i=1:5 
# 		u_mapped = Systems.ikedamap_rule(u_mapped, p, n)
# 	end
#     return SVector{2}(u_mapped[1], u_mapped[2])
# end


@inbounds function ikedamap_5thorder_rule(u,p,n)
	for i=1:5 
		u = Systems.ikedamap_rule(u, p, n)
	end
    return SVector{2}(u[1], u[2])
end


#test if 5th iterate map works
# Ttr=0
# T=200
# ik = Systems.ikedamap(u0; a, b, c, d)
# ik_5order = DiscreteDynamicalSystem(ikedamap_5thorder_rule, u0, p)
# tr = trajectory(ik, T, u0; Ttr); t = Ttr:Ttr+T
# tr_5th = trajectory(ik_5order, T, u0; Ttr); t_5th = Ttr:Ttr+T/5
# tr = tr[1:5:end, :]; t = t[1:5:end]

#find fp: first look through 'periodicorbits`, which is faster. Give ics around fp by eyeballing. If it does indeed find a fp that is close, use that to give refined ics to `fixedpoints` which is then faster. 
xs = range(-0.15, -0.148, length=25)
ys = range(-0.112, -0.111, length=25)
ics = [[x,y] for x ∈ xs for y ∈ ys]
ik_5order = DiscreteDynamicalSystem(ikedamap_5thorder_rule, u0, p)
fps = periodicorbits(ik_5order, 1, ics)
# fp_3 = [0.0879, -0.3778]
tr = trajectory(ik_5order, 100, fp_leftmost) #not very good, escapes after ~2 iterates
scatter!(ax, tr[:,1], tr[:,2], color=:orange, markersize=10)
fp_leftmost = [ -0.149036, -0.111942]
# fp_leftmost =  [-0.149, -0.1119]
# fp_leftmost =  [-0.156, -0.2293]

#using fixedpoints
x = interval(-0.14904, -0.14902)
y = interval(-0.11195, -0.11194);
box = x × y
fps, eigs, stable = fixedpoints(ik_5order, box; tol=1e-20)
tr = trajectory(ik_5order, 100, fps[1]) #much better!
fps_5 = trajectory(ik, 4, fp_leftmost)




T = 3e3
Ttr = 1e3
fig = Figure(resolution=(1920, 1080)) 
axs = []
# for (i, d) ∈ enumerate([7.1, 7.2])
d = 7.1
p = [a, b, c, d]
i = 1
ik = Systems.ikedamap(;a, b, c, d)
tr = trajectory(ik, T; Ttr)
x = tr[:,1]; y = tr[:,2]; t = Ttr:Ttr+T
measure = :black
ax = Axis(fig[1,i], ylabel="y", xlabel="x", title="d=$(d)")
scatter!(ax, x, y, markersize=3, color=measure)
xlims!(-0.6, 1.7)
ylims!(-1.8, 1.4)
scatter!(ax, fps_5[:,1], fps_5[:,2], color=:red, markersize=10)
# end

save("$(plotsdir())/$(dirname)/interiorcrisis-ikedamap-period5saddle-d_$(d).png", fig)
fps_5 = [-0.149036 -0.111942; 0.6901164331310736 -0.07534279207228838; 0.7081550024831971 0.6107258394579234; 0.08793801081244812 -0.377787118392552; 1.0755007518118098 -0.25769927892208205]
 



# ax = Axis(fig[4, 1], ylabel="x")
# lines!(t, tr[:,1], color=:black)
# lines!(t_5th, tr_5th[:,1], color=:orange)
# ax = Axis(fig[5, 1], ylabel="y", xlabel="t")
# lines!(t, tr[:,2], color=:black)
# lines!(t_5th, tr_5th[:,2], color=:orange)

# grid = [[x,y] for x ∈ range(-0.5, 1.5, length=10) for y ∈ range(-2, 1, length=10)]
# grid = reduce(hcat, grid)
# z = mapslices(x->ikedamap_5thorder_rule(x, p, 1), grid, dims=1)

# heatmap(grid[1,:], grid[2,:], z)








# ax = Axis(fig[1,1], ylabel="y", xlabel="x"); 
# scatter!(ax, x, y, markersize=3, color=measure)


# --------------------- estimate durations in inner core --------------------- #
# fps_5 = [-0.149036 -0.111942; 0.6901164331310736 -0.07534279207228838; 0.7081550024831971 0.6107258394579234; 0.08793801081244812 -0.377787118392552; 1.0755007518118098 -0.25769927892208205]
using CurveFit
T = 1000e3; Ttr = 1e3
a = 0.84; b = 0.9; c = 0.4 
#get inner core
d = 7.22
ik = Systems.ikedamap(;a, b, c, d)
traj = trajectory(ik, T; Ttr); t = Ttr:Ttr+T
innercore = deepcopy(Matrix(traj))
#integrate real system
d = 7.3
ik = Systems.ikedamap(;a, b, c, d)
traj = trajectory(ik, T; Ttr); t = Ttr:Ttr+T


maxdist_threshold = 0.1
maxdist_threshold = 0.04
arewithin = mapslices(x->withinset(x, innercore, maxdist_threshold), Matrix(traj), dims=2)[:,1]; #a bit expensive, have to compare each point to all other points in the set

T, bool_innercore = length_samevalues_allowfluctuations(arewithin)
τs = T[1]
numbins = 20; weights, bins = histogram(τs, numbins); 

fig = Figure(resolution=(800, 600), fontsize=30,figure_padding=(5, 35, 5, 30))
ax = Axis(fig[1,1], yscale=log10, ylabel="P(τ)", xlabel="τ")
scatterlines!(ax, bins[1:end-1], weights, color=:black)

xfit = bins[1:end-1]; yfit = weights[1:end]; A, B = exp_fit(xfit, yfit)
l=lines!(ax, xfit, A .* exp.(B .* xfit), color=:red, label="P(τ) = $(round(A, digits=4)) exp($(round(B, digits=4)) τ)" )
# fig[2,1] = Legend(fig, ax, orientation= :horizontal, tellwidth = false, tellheight = true)
axislegend(ax)
xlims!(ax, -20, 1500)
save("$(plotsdir())/$(dirname)/interiorcrisis-ikedamap-period5saddle-d_$(d)-numbins_$(numbins).png", fig)
using DelimitedFiles
writedlm("$(plotsdir())/$(dirname)/durationininnercore-d_$(d).dat", τs)
fig = Figure(resolution=(1920, 1080)) 
ax = Axis(fig[1,1], ylabel="P(τ)", xlabel="τ", title="d=$(d)", xscale=log10)
hist!(ax, durationinneredge, 50)
ax = Axis(fig[1,1], ylabel="y", xlabel="x", title="d=$(d)")
scatter!(ax, traj[:,1], traj[:,2], markersize=3, color=measure)
ax = Axis(fig[2,1], ylabel="y", xlabel="x", title="d=$(d)")
lines!(ax, traj[:,1])
