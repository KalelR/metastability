using DrWatson
@quickactivate "metastability"
using GLMakie, OrdinaryDiffEq
include("$(scriptsdir())/utils.jl")

# @inbounds @inline function modvanderpol!(du, u, p, t)
# 	ϵ = p[1]
# 	x, y = u
#     du[1] = y 
#     du[2] = -x + ϵ * ((y^3/3) -y)
# end



# ϵ = 1
# T = 2e3; Ttr= 500
# Δt = 0.1
# u0s = [[0.1, 0.1], [2.0, 0.0]]
# fig = Figure(resolution=(800, 800), fontsize = 20)
# ax1 = Axis(fig[1, 1], xlabel="t", ylabel="x")
# ax2 = Axis(fig[2, 1], xlabel="t", ylabel="xdot")
# for u0 ∈ u0s
#     p = [ϵ]
#     solver = Tsit5()
#     tspan = (0, T)
#     prob = ODEProblem(modvanderpol!, u0, tspan, p)
#     sol = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, maxiters=1e9)
#     t = sol.t ; 
#     tr = sol[:,:]'
#     lines!(ax1, t, tr[:, 1])
#     lines!(ax2, t, tr[:, 2])
# end
# linkxaxes!(ax1, ax2)



# """ prey = N; predator = P """
@inbounds @inline function predator_prey!(du, u, p, t)
	α, γ, ϵ, ν, h, K, m = p
	N, P = u
    du[1] = α*N*(1 - N/K) - γ*N*P / (N+h)
    du[2] = ϵ * (  ν*γ*N*P/(N+h) - m*P    )
	# return SVector{2}(du1, du2)
end


γ = 2.5
h = 1 
ν = 0.5 
m = 0.4
ϵ = 1.0
T = 250
Ttr= 100
Δt = 0.01
α = 1.5 
K = 2.2 
p = [α, γ, ϵ, ν, h, K, m]
solver = Tsit5()
tspan = (0, T)


fig = Figure(resolution=(1200, 800), fontsize = 20)
ax1 = Axis(fig[1, 1], xlabel="t", ylabel="x")
ax2 = Axis(fig[2, 1], xlabel="t", ylabel="y")
ax3 = Axis(fig[1:2, 2:3], xlabel="x", ylabel="y")
u0 = [0.1, 0.1]
prob = ODEProblem(predator_prey!, u0, tspan, p)
for perturbation ∈ [0, 2]
    condition(u,t,integrator) = t==126
    affect!(integrator) = integrator.u[1] += perturbation
    cb = DiscreteCallback(condition,affect!)

    sol = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, callback=cb, tstops=[126.0])
    t = sol.t; tr = sol[:,:]'
    lines!(ax1, t, tr[:,1])
    lines!(ax2, t, tr[:,2])
    lines!(ax3, tr[:,1], tr[:,2])
end

sol_att = solve(prob, solver, saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8)
t_att = sol_att; tr_att = sol_att[:,:]'

u0 = [1.2, 0.8]
prob = ODEProblem(predator_prey!, u0, tspan, p)
sol = solve(prob, solver, saveat=0:Δt:T, abstol=1e-8, reltol=1e-8)
t = sol; tr = sol[:,:]'

fig = Figure(resolution=(1200, 800), fontsize = 20)
ax = Axis(fig[1, 1], xlabel="x", ylabel="y")
lines!(ax, tr_att[:,1], tr_att[:,2], color=(:gray, 0.2))
lines!(ax, tr[:,1], tr[:,2], color=:orange)


# ----------------------------------- MOVIE ---------------------------------- #


Tplot = 70
framerate=400
Δtanimation = 0.01
# Δtanimation = 1
t_plot, tr_plot, speeds_plot, frames = animationdata(sol, Tplot, Δt)
colors = pointspeed_as_colors(speeds_plot);

points = Observable(Point2f[(tr_plot[1,1], tr_plot[1,2])])
colors_ob = Observable([colors[1]])
time = Observable(t_plot[1])
fig = Figure(resolution=(800, 1000))
ax = Axis(fig[1,1]);
lines!(ax, tr_att[:,1], tr_att[:,2], color=(:black, 1.0))
scatter!(ax, points, color=(:orange, 0.6), markersize=6)
# hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
hidedecorations!(ax)
hidespines!(ax)
limits!(ax, 0, 1.35, 0.25, 1.2)
hidespines!(ax, :t, :r) 
record(fig, "$(plotsdir())/stability/stabilityexample.mp4", frames;
        framerate) do frame
    new_point = Point2f(tr_plot[frame,1], tr_plot[frame,2])
    points[] = push!(points[], new_point)
	colors_ob[] = push!(colors_ob[], colors[frame])
    time[] = t_plot[frame]
end


# ---------------------------- SCALING EXPONENTIAL --------------------------- #
f(x, A, α) = A*exp(-α*x) 
xs = collect(range(0, 1e3, length=100))
ys = f.(xs, 1000, 1e-2)
fig=Figure()
ax=Axis(fig[1,1], yscale=log10)
# lines!(ax, xs, ys.+1e-5, color=:black)
scatter!(ax, xs, ys.+1e-5, color=:black)