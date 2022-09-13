using DrWatson
@quickactivate "metastability"
using GLMakie,  DifferentialEquations

include("$(scriptsdir())/utils.jl")

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
	traj = trajectory(lo, T; Ttr)
	x = traj[:,1];  t = Ttr:Ttr+T
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
traj = trajectory(ik, T; Ttr); t = Ttr:Ttr+T
x = traj[:,1]; y = traj[:,2]
measure = histmeasure(traj, 2)

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
	traj = trajectory(ik, T, u0; Ttr)
	x = traj[:,1]; y = traj[:,2]; t = Ttr:Ttr+T
	# measure = histmeasure(traj, 2)
	measure = :black
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
#just tanim-series
for (i, a) ∈ enumerate([0.997, 1.003])
	ik = Systems.ikedamap(;a, b, c, d)
	traj = trajectory(ik, T, u0; Ttr)
	x = traj[:,1]; y = traj[:,2]; t = Ttr:Ttr+T
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
	# save("../plots/boundarycrisis/boundarycrisis-ikedamap-timeseries-a_$(a).png", fig, px_per_unit=4)

end



# ---------------------- distribution of times in saddle --------------------- #
include("$(scriptsdir())/utils.jl")
T = 2e5
Ttr = 0
b = 0.9 
c = 0.4 
d = 6.0
# a = 1.0 # a=A=r varying
a = 1.003
u0 = [-2.0, -2.0]
ik = Systems.ikedamap(u0;a, b, c, d)
# traj = trajectory(ik, T, u0; Ttr);  ts = Ttr:Δt:Ttr+T
using DynamicalSystems, ProgressMeter
tinteg = tangent_integrator(ik, 2)
λs_vec, traj_vec, t= lyapunovspectrum_convergence_discrete(tinteg, Int(T), 1, Ttr, true)
λs = reduce(hcat, λs_vec)'
traj = reduce(hcat, traj_vec)'

fig = Figure()
ax = Axis(fig[1,1])
scatter!(ax, traj[:,1], traj[:,2])
ax = Axis(fig[2,1])
lines!(ax, t, traj[:,1])
ax = Axis(fig[3,1])
lines!(ax, t, λs[:,1])
lines!(ax, t, λs[:,2])

## just by trajectories entering neighborhood of FP; its NOT VERY ELEGANT! haha but a quick way to do it
fp = [ 2.9715737943410305, 4.153134755539537 ]
iswithinneighborhood(point, fp, threshold) = evaluate(Euclidean(), point, fp) ≤ threshold ? true : false  
iswithinneighborhood([3.0, 3.0], fp, 0.1)

using DynamicalSystemsBase:obtain_access, get_state
function trajectory_discrete(integ, t, u0 = nothing;
        Δt::Int = 1, Ttr = 0, a = nothing, diffeq = nothing, fp, threshold
    )
    !isnothing(u0) && reinit!(integ, u0)
    Δt = round(Int, Δt)
    T = eltype(get_state(integ))
    t0 = current_time(integ)
    tvec = (t0+Ttr):Δt:(t0+t+Ttr)
    L = length(tvec)
    T = eltype(get_state(integ))
    X = isnothing(a) ? dimension(integ) : length(a)
    data = Vector{SVector{X, T}}(undef, L)
    Ttr ≠ 0 && step!(integ, Ttr)
    data[1] = obtain_access(get_state(integ), a)
	t_convergence = 0
    for i in 2:L
        step!(integ, Δt)
		if iswithinneighborhood(integ.u, fp, threshold) for j=i:L data[j] = SVector{X, T}(fp) end; t_convergence=i; break end
        data[i] = SVector{X, T}(obtain_access(get_state(integ), a))
    end
    return Dataset(data), t_convergence
end

function time_to_converge(integ, t, u0 = nothing;
	Δt::Int = 1, Ttr = 0, a = nothing, diffeq = nothing, fp, threshold
)
	!isnothing(u0) && reinit!(integ, u0)
	Δt = round(Int, Δt)
	t0 = current_time(integ)
	tvec = (t0+Ttr):Δt:(t0+t+Ttr)
	L = length(tvec)
	Ttr ≠ 0 && step!(integ, Ttr)
	t_convergence = 0
	for i in 2:L
		step!(integ, Δt)
		if iswithinneighborhood(integ.u, fp, threshold)  return i; break end
	end
	return L
end

T=1e6; Ttr=1e2
threshold = 1
xs = range(-0.0, 1.0, length=100)
ys = range(-2.0, 0, length=100)
u0s = [[x, y] for x ∈ xs for y ∈ ys]
τs = zeros(Int, length(u0s))
for (i, u0) ∈ enumerate(u0s)
	integ = integrator(ik, u0) 
	t_conv = time_to_converge(integ, T, u0; Ttr, fp, threshold)
	τs[i] = t_conv
end

filter!(x->x>=10, τs)
filter!(x->x<T, τs)
# numbins = 10; weights, bins = histogram(τs, numbins); 
numbins = 15;mode = :pdf
bins = collect(range(minimum(τs), maximum(τs), length=numbins))
 weights, bins = histogram(τs, bins; mode=mode); 

fig = Figure(resolution=(800, 600), fontsize=30,figure_padding=(5, 35, 5, 30))
ax = Axis(fig[1,1], yscale=log10, ylabel="P(τ)", xlabel="τ")
scatterlines!(ax, bins[1:end-1], weights, color=:black)

xfit = bins[1:end-1]; yfit = weights[1:end]; A, B = exp_fit(xfit, yfit)
l=lines!(ax, xfit, A .* exp.(B .* xfit), color=:red, label="P(τ) = $(round(A, digits=6)) exp($(round(B, digits=6)) τ)" )
axislegend(ax)
# xlims!(ax, -20, 1500)
save("$(plotsdir())/$(dirname)/boundarycrisis-ikedamap-distributiontimeschaoticsaddle-a_$(a)-b_$(b)-c_$(c)-d_$(d)-numbins_$(numbins).png", fig)

using DelimitedFiles
writedlm("$(plotsdir())/$(dirname)/durationinchaoticsaddle-d_$(d).dat", τs)

#------------------------------------------- animation ----------
include("$(scriptsdir())/utils.jl")
# a = 0.997
a = 1.003 #after crisis, CS
T = 101500
Ttr = 100500

ik = Systems.ikedamap(;a, b, c, d)
traj = trajectory(ik, T, u0; Ttr)
x = traj[:,1]; y = traj[:,2]; t = Ttr:Ttr+T
sol = [ [traj[i,1], traj[i,2]] for i=1:length(traj)]
p = [a, b, c, d]
speeds = norm.([ik.jacobian([traj[i,1], traj[i,2]], p, 1) for i=1:length(traj)])

Δt = 1.0
Tplot= T-Ttr
framerate=25
# framerate=50
tplot = Int64(Tplot/Δt);
Δt_plot = round(Int64, Δt/Δt);
t_plot = t[1:Δt_plot:tplot];
tr_plot = traj[1:Δt_plot:tplot, :];
speeds_plot = log10.(1 ./ speeds[1:Δt_plot:tplot]);
frames = 2:length(t_plot)
duration = length(frames) / framerate

#transform the vector with the info for the colors onto a Int vector going from 1 to 264; this is used to index the colormap (wihch has 264 colors); basically transforming it into an vector of indices
v = (speeds_plot .- minimum(speeds_plot)) ./ (maximum(speeds_plot) .- minimum(speeds_plot)) .* 255 .+ 1;
v = round.(Int, v );
colors = ColorSchemes.viridis[v];


# points = Observable(Point2f[(tr_plot[1,1], tr_plot[1,2])])
# colors_ob = Observable([colors[1]])
# tanim = Observable(t_plot[1])
# fig = Figure(resolution=(800, 600))
# ax = Axis(fig[1,1], title= @lift("t = $($tanim)"))
# s=scatter!(ax, points, color=colors_ob, markersize=3)
# hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
# xlims!(-0.6, 6.2)
# ylims!(-2.8, 6)

# record(fig, "../plots/boundarycrisis/ikedmap-animation-a_$(a)-b_$(b)-c_$(c)-d_$(d).mp4", frames;
#         framerate) do frame
#     tanim[] = t_plot[frame]
#     new_point = Point2f(tr_plot[frame,1], tr_plot[frame,2])
# 	points[] = push!(points[], new_point)
# 	colors_ob[] = push!(colors_ob[], colors[frame])
# end



#without

c1 = "#440154FF"
c2 = "#FDE725FF"
cball = "#FF1400"
points = Observable(Point2f[(tr_plot[1,1], tr_plot[1,2])])
points2 = Observable(Point2f[(t_plot[1], tr_plot[1,1])])
tanim = Observable(t_plot[1])
fig = Figure(resolution=(800, 600))
ax = Axis(fig[1,1], title= @lift("t = $($tanim)"))
s=scatter!(ax, tr_plot[:,1], tr_plot[:,2], color=[t_plot_el < 101330 ? :green : :purple for t_plot_el in t_plot], markersize=2)
s=scatter!(ax, points, color=:red, markersize=16)
hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
xlims!(-0.6, 6.2)
ylims!(-2.8, 6)

ax = Axis(fig[2, 1], ylabel="x", xlabel="t")
# lines!(ax, t_plot, tr_plot[:,1], color=[el > 0 ? c1 : c2 for el in tr_plot[:,1]])
lines!(ax, t_plot, tr_plot[:,1], color=[t_plot_el < 101330 ? :green : :purple for t_plot_el in t_plot])
# lines!(ax, t_plot, tr_plot[:,1], color=:black)
scatter!(ax, points2, color=cball, markersize=16)

record(fig, "$(plotsdir())/boundarycrisis/ikedamap-animation-timeseries-a_$(a)-b_$(b)-c_$(c)-d_$(d).mp4", frames;
        framerate) do frame
    tanim[] = t_plot[frame]
    new_point = Point2f(tr_plot[frame,1], tr_plot[frame,2])
    new_point2 = Point2f(t_plot[frame,1], tr_plot[frame,1])
	points[] = [new_point]
	points2[] = [new_point2]
end



