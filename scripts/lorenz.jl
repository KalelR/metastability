using DrWatson
@quickactivate "metastability"

using GLMakie, DynamicalSystems

T = 1000
Ttr=100
Δt = 0.005

ρ = 28 #normal CA 
# ρ = 166 #period (LC)
# ρ = 166.06 #LC
# ρ = 166.068 #intermittency (CA)
# ρ = 166.3 #intermittency (CA)
# ρ = 166.3 #intermittency (CA)
lo = Systems.lorenz(; σ=10, β=8/3, ρ) #lorenz 63; σ==a; β==b; ρ==r
tr = trajectory(lo, T; Δt, Ttr); t = 0:Δt:T

fig = Figure()
ax1 = Axis(fig[1, 1])
lines!(ax1, t, tr[:, 3], linewidth=1)
# save("lorenz-zvariable-r_$(ρ).png", fig)



fig = Figure(resolution=(1000,600), fontsize = 20)
ax1 = Axis3(fig[1, 1], aspect=(1,1,0.5))
scatter!(ax1, tr[:,1], tr[:,2], tr[:,3], color=t, colormap=:thermal, markersize=1000)
save("lorenzattractor-r_$(ρ).png", fig)
# save("$(plotsdir())/lorenzattractor-r_$(ρ).png", fig)


#plot everthign
fig = Figure(resolution=(1000,600), fontsize = 20)
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])
ax3 = Axis(fig[1, 3])
ax4 = Axis3(fig[2, :])
lines!(ax1, t, tr[:, 1], linewidth=1, color=:black)
lines!(ax2, t, tr[:, 2], linewidth=1, color=:black)
lines!(ax3, t, tr[:, 3], linewidth=1, color=:black)
xlims!(ax1, 900,1000)
xlims!(ax2, 900,1000)
xlims!(ax3, 900,1000)
scatter!(ax4, tr[:,1], tr[:,2], tr[:,3], color=t, colormap=:thermal, markersize=1000)
save("lorenzdynamics-r_$(ρ).png", fig)


#points remain, so they accumualte
fig = Figure(resolution=(1000,600), fontsize = 20)
ax = Axis3(fig[1,1], aspect = (1,1,0.5), azimuth = -0.3π, elevation = π/9)
lines!(ax, tr[:,1], tr[:,2], tr[:,3], color=(:black, 0.3))
points = Observable(Point3f[])
color = Observable(Float64[])
# attr = (color=color, transparency=true, colormap=:plasma)
attr = (color=:red, transparency=true, colormap=:plasma)
# pltobj = lines!(ax, points; attr...)
scatter!(ax, points; markersize = 0.02, markerspace=SceneSpace, attr...)
# cbar = Colorbar(fig[1,2], pltobj, label = "t", width = 15, ticksize=15, tickalign = 1, height = Relative(0.5))
path = "lorenzattractor-r_$(ρ).mp4"
record(fig, path, enumerate(t), framerate = 24*16) do (i, t)
    xyz = Point3f(tr[i, :])
    push!(points[], xyz)
    # push!(color[], t)
    notify(points)
    # notify(color)
    # autolimits!(ax)
end

#only one points
fig = Figure(resolution=(1000,600), fontsize = 20)
ax = Axis3(fig[1,1], aspect = (1,1,0.5), azimuth = -0.3π, elevation = π/9)
lines!(ax, tr[:,1], tr[:,2], tr[:,3], color=(:black, 0.05), priority=0.1)
points = Observable(Point3f[])
attr = (color=:red, transparency=true, colormap=:plasma)
scatter!(ax, points; markersize = 0.04, markerspace=SceneSpace, attr..., priority=1)
path = "lorenzattractor-r_$(ρ).mp4"
record(fig, path, enumerate(t), framerate = 24*32) do (i, t)
    xyz = Point3f(tr[i, :])
	points[] = [xyz]
    # notify(points)
end



#estimating duration on each state 

## i'll just count the time in x<0 (one wing) and x>0 (another wing); consider the transition to be fast enough that it doesnt mattter much 

function count_condition(vec, value)
	sum = 0
	durations = Int64[]
	for x ∈ vec
		if x > value
			sum += 1
		else
			if sum > 0
				push!(durations, sum)
				sum = 0
			end
		end 
	end
	return durations
end
	
	
T = 1000000
Ttr=500
Δt = 0.1

ρ = 28 #normal CA 
lo = Systems.lorenz(; σ=10, β=8/3, ρ) #lorenz 63; σ==a; β==b; ρ==r
tr = trajectory(lo, T; Δt, Ttr); t = 0:Δt:T
durations = count_condition(tr[:,1], 0)
t_durations = durations .* Δt

using StatsBase
edges = 0:1.0:15
a= fit(Histogram, t_durations, edges )
weights = a.weights
y_offset  = 1
weights .+= y_offset

fig = Figure()
ax = Axis(fig[1,1], yscale=log10)
scatter!(ax, edges[1:end-1], weights)



using CurveFit 
fit_res = linear_fit(edges[1:end-1], log10.(weights))
offset = fit_res[1]; grad = fit_res[2]
y_fit = offset .+ grad .* edges[1:end-1]
x = collect(edges[1:end-1])
x_fit = x[y_fit .> 0]
y_fit = y_fit[y_fit .> 0]
lines!(ax, x_fit, 10 .^ (y_fit), color=:black, label = "y = 10^$(grad) + $(offset)")
Legend()

# should seriously improve this before showing anyone; but point is that it's exponential



# -------------------------------------------------------------------------------------------------------------------------------------------------- #
#                                                                  chaotic transient                                                                 #
# -------------------------------------------------------------------------------------------------------------------------------------------------- #
T = 200
Ttr=0
Δt = 0.01

# ρ = 13 #only stable foci
ρ = 22 #chaotic saddle + stable node 
# ρ = 23 #chaotic saddle + stable node 
u0 = [5.0,5,5]
lo = Systems.lorenz(u0; σ=10, β=8/3, ρ) #lorenz 63; σ==a; β==b; ρ==r
tr = trajectory(lo, T; Δt, Ttr); t = Ttr:Δt:T
fig = Figure(resolution=(4000, 2000), fontsize = 20)
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])
ax3 = Axis(fig[1, 3])
ax4 = Axis3(fig[2, :])
lines!(ax1, t, tr[:, 1], linewidth=1, color=:black)
lines!(ax2, t, tr[:, 2], linewidth=1, color=:black)
lines!(ax3, t, tr[:, 3], linewidth=1, color=:black)
xlims!(ax1, Ttr, T)
xlims!(ax2, Ttr, T)
xlims!(ax3, Ttr, T)
# scatter!(ax4, tr[:,1], tr[:,2], tr[:,3], color=t, colormap=:thermal, markersize=1000)
lines!(ax4, tr[:,1], tr[:,2], tr[:,3], color=t, colormap=:thermal, markersize=1000)
save("lorenzdynamics-r_$(ρ).png", fig)
