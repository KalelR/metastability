using DrWatson
@quickactivate "metastability"
using GLMakie, OrdinaryDiffEq
include("$(scriptsdir())/utils.jl")


@inbounds function heteroclinic_cycle_gh!(du, u, p, t)
	μ, a, b, c = p
    du[1] = μ*u[1] + u[1]*(a*u[1]^2 + b*u[2]^2 + c*u[3]^2)
    du[2] = μ*u[2] + u[2]*(a*u[2]^2 + b*u[3]^2 + c*u[1]^2)
    du[3] = μ*u[3] + u[3]*(a*u[3]^2 + b*u[1]^2 + c*u[2]^2)
    nothing
end

@inbounds function jac!(J, u, p, t)
    x,y,z = u
    l, a, b, c = p
    J[1,1] = l + 3*a*x^2 + b*y^2 + c*z^2
    J[1,2] = 2*b*x*y
    J[1,3] = 2*c*x*z
    J[2,1] = 2*c*y*x
    J[2,2] = l+3*a*y^2 + b*z^2 + c*x^2
    J[2,3] = 2*b*y*z
    J[3,1] = 2*b*z*x
    J[3,2] = 2*c*z*y
    J[3,3] = l+3*a*z^2 + b*x^2 +c*y^2
    nothing
end



# T = 5000
T = 10000
μ = 1.0 

##GH - stable HC!!!!!
a = -0.2
b = -0.1
# a = -0.3 
# b = -0.2
# a = -0.5
# b = -0.01
c = -1 - a - b

#inpisred by ashwin (with fps being (1,0,0)) but does not work
# a = -1.0
# b = -1.0
# c = -1.0 

# ##or 
# μ = 0.2 
# a = b = c = -μ

## stabilitizes and becomes node apparently
# T = 1500
# a = -1
# b = -1.5
# c = (2a - b) - 0.1 #unstable ## 
# c = (2a - b) + 0.1#not even a HC 


## parameters from scholarpedia (adapted so c>a>b); unstable
# c = -0.55 
# a = -1.0
# b = -1.5

## directly from scholarpedia; unstable
# a = -1.0 
# b = -0.55 
# c = -1.5

@warn c > a > b #Krupa exists if true 
@warn 2a > b + c #Krupa stable if true
@warn a + b + c == -1 #GH
@warn -1/3 < a < 0 #GH 
@warn c < a < b < 0 #GH


p = [μ, a, b, c]


xpcoor(μ,p) = sqrt(-μ/p)
fps_x = [[xpcoor(μ, a), 0, 0], [xpcoor(μ, b), 0, 0], [xpcoor(μ, c), 0, 0]]; fps_x = [fps_x; -1 .* fps_x]
fps_y = [[0, xpcoor(μ, a), 0], [0, xpcoor(μ, b), 0], [0, xpcoor(μ, c), 0]]; fps_y = [fps_y; -1 .* fps_y]
fps_z = [[0, 0, xpcoor(μ, a)], [0, 0, xpcoor(μ, b)], [0, 0, xpcoor(μ, c)]]; fps_z = [fps_z; -1 .* fps_z]
fps = hcat([fps_x; fps_y; fps_z]...)'

# u0 = fps_y[3] .+ [0.1, -0.5, 0.1] #still converges to node, this fp is repeller
# u0 = fps_y[1] .+ [0.1, +0.5, 0.1]
u0 = [0.1, 0.2, 0.3] #off the coordinate plane or axis should be attracted towards HC
# u0 = 5 .* rand(3)
# u0 = [5.3, 5.2, 5.4]
# u0 = [-0.1, -0.2, -0.3] 
# u0 = [10,10,10] #goes to a node
# u0 = [10,-10,10] #goes to a node

tspan = (0, T)
# hcgh = ODEProblem(heteroclinic_cycle_gh!, u0, tspan, p)
ff = ODEFunction(heteroclinic_cycle_gh!; jac=jac!)
hcgh = ODEProblem(ff, u0, tspan, p)
# sol = solve(hcgh, saveat=0:0.01:T, abstol=1e-16, reltol=1e-16); t = sol.t;
# sol = solve(hcgh, AutoTsit5(Rosenbrock23()), saveat=0:0.01:T, abstol=1e-8, reltol=1e-8); t = sol.t;
# sol = solve(hcgh, Rodas5(), saveat=0:0.01:T, abstol=1e-10, reltol=1e-10, maxiters=1e9); t = sol.t;
# sol = solve(hcgh, VCABM(), saveat=0:0.01:T, abstol=1e-10, reltol=1e-10, maxiters=1e9); t = sol.t;
sol = solve(hcgh, RK4(), dt=1/1000); t = sol.t;
# sol = solve(hcgh, saveat=0:0.01:T); t = sol.t;
# sol = solve(hcgh, RadauIIA3(), saveat=0:0.01:T, abstol=1e-10, reltol=1e-10, maxiters=1e9); t = sol.t;


fig = Figure()
ax1 = Axis(fig[1, 1])
lines!(ax1, t, sol[1,:])
lines!(ax1, t, sol[2,:])
lines!(ax1, t, sol[3,:])
save("../plots/heterocliniccycle/guckenheimerholmes.png", fig)

# fig2 = Figure()
ax2 = Axis3(fig[1,2])
# scatter(sol[1,:], sol[2,:], sol[3,:])
lines!(ax2, sol[1,:], sol[2,:], sol[3,:], color=t)
scatter!(ax2, fps[:,1], fps[:,2], fps[:,3], color=:red, markersize=10)
scatter!(ax2, [0], [0], [0], color=:blue, markersize=13)


#------------------------movie ---------------------------
fps_x = [[xpcoor(μ, a), 0, 0]];
fps_y = [[0, xpcoor(μ, a), 0]];
fps_z = [[0, 0, xpcoor(μ, a)]];
fps = hcat([fps_x; fps_y; fps_z]...)'


include("utils.jl")
plotsdir() = "../plots/"
Ttr = 240
T = 420
Ttr = 0.0; T = 240; 
Δt= 0.1
tspan = (0, T)
ff = ODEFunction(heteroclinic_cycle_gh!; jac=jac!)
hcgh = ODEProblem(ff, u0, tspan, p)
sol = solve(hcgh, RK4(), dt=1/1000, saveat=Ttr:Δt:T); t = sol.t; tr = sol[:, :]';
Tplot= T-Ttr
framerate = 150 #166.3
tplot = Int64(Tplot/Δt);
Δt_plot = Int64(Δt/Δt);
t_plot = t[1:Δt_plot:tplot];
tr_plot = tr[1:Δt_plot:tplot, :];
speeds = norm.(timederivative(sol))
speeds_plot = log10.(1 ./ speeds[1:Δt_plot:tplot]);
frames = 2:length(t_plot);

az =0.7875222039230623 
el =0.3190658503988664 

#transform the vector with the info for the colors onto a Int vector going from 1 to 264; this is used to index the colormap (wihch has 264 colors); basically transforming it into an vector of indices
# v = (speeds_plot .- minimum(speeds_plot)) ./ (maximum(speeds_plot) .- minimum(speeds_plot)) .* 255 .+ 1;
v = (t .- minimum(t)) ./ (maximum(t) .- minimum(t)) .* 255 .+ 1;
v = round.(Int, v );
colors = ColorSchemes.viridis[v];

points = Observable(Point3f[(tr[1,1], tr[1,2], tr[1,3])])
colors_ob = Observable([colors[1]])
fig = Figure(resolution=(800, 600))
time_p = Observable(t_plot[1])
ax = Axis3(fig[1,1], azimuth = az, elevation = el, title= @lift("t = $($time_p)"))
# scatter!(ax, points, color=colors_ob)
lines!(ax, points, color=colors_ob)
scatter!(ax, fps[:,1], fps[:,2], fps[:,3], color=:red, markersize=10)
# scatter!(ax, [0], [0], [0], color=:blue, markersize=13)
hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
limits!(ax, 0.0, 2.4, 0, 2.4, 0, 2.4)
# hidespines!(ax, :t, :r) 
record(fig, "../plots/heterocliniccycle/guckenheimerholmes-Ttr_$(Ttr).mp4", frames;
        framerate) do frame
    time_p[] = t_plot[frame]
    new_point = Point3f(tr_plot[frame,1], tr_plot[frame,2], tr_plot[frame,3])
    points[] = push!(points[], new_point)
	colors_ob[] = push!(colors_ob[], colors[frame])
end


#------------------------movie gray background ---------------------------
fps_x = [[xpcoor(μ, a), 0, 0]];
fps_y = [[0, xpcoor(μ, a), 0]];
fps_z = [[0, 0, xpcoor(μ, a)]];
fps = hcat([fps_x; fps_y; fps_z]...)'


Ttr = 240
T = 420
Ttr = 0.0; T = 240; 
Δt= 0.1
tspan = (0, T)
ff = ODEFunction(heteroclinic_cycle_gh!; jac=jac!)
hcgh = ODEProblem(ff, u0, tspan, p)
sol = solve(hcgh, RK4(), dt=1/1000, saveat=Ttr:Δt:T); t = sol.t; tr = sol[:, :]';
Tplot= T-Ttr
framerate = 150 #166.3


t_plot, tr_plot, speeds_plot, frames = animationdata(sol, Tplot, Δt, Δt)
colors = pointspeed_as_colors(speeds_plot);

az =0.7875222039230623 
el =0.3190658503988664 

points = Observable(Point3f[(tr_plot[1,1], tr_plot[1,2], tr_plot[1,3])])
colors_ob = Observable([colors[1]])
tanim = Observable(t_plot[1])

fig = Figure(resolution=(800, 600))
ax = Axis3(fig[1:2,1], azimuth = az, elevation = el, title= @lift("t = $($tanim)"))
lines!(ax, tr_plot[:,1], tr_plot[:,2], tr_plot[:,3], color=colors)
scatter!(ax, points, color=:orange)
scatter!(ax, fps[:,1], fps[:,2], fps[:,3], color=:red, markersize=10)
limits!(ax, 0.0, 2.4, 0, 2.4, 0, 2.4)
hidedecorations!(ax, ticks=false, label=false, ticklabels=false)

ax = Axis(fig[3,1])
lines!(ax, t_plot, tr_plot[:,1], color=(:gray, 0.4))
scatter!(ax, points2, color=(:orange, 1.0))
hidespines!(ax, :t, :r) 


record(fig, "$(plotsdir())/heterocliniccycle/guckenheimerholmes-Ttr_$(Ttr)-graybackground.mp4", frames;
        framerate) do frame
    tanim[] = t_plot[frame]
    new_point = Point3f(tr_plot[frame,1], tr_plot[frame,2], tr_plot[frame,3])
    new_point2 = Point2f(t_plot[frame], tr_plot[frame,1])
    points[] = [new_point]
    points2[] = [new_point2]
	# colors_ob[] = push!(colors_ob[], colors[frame])
endz



## Stability analysis
using LinearAlgebra
xp = sol[end]
xp = fps_y[2]
# xp = fps_y[2]
# xp = fps_y[3]
# xp = [sqrt(-1/a), 0, 0]
J = gh_jac(xp, p)
evals, evecs = eigen(J)


# With noise
@inbounds function heteroclinic_cycle_gh_noise(du, u, p, t)
	μ, a, b, c, n = p
    du[1] = μ*u[1] + u[1]*(a*u[1]^2 + b*u[2]^2 + c*u[3]^2)
    du[2] = μ*u[2] + u[2]*(a*u[2]^2 + b*u[3]^2 + c*u[1]^2)
    du[3] = μ*u[3] + u[3]*(a*u[3]^2 + b*u[1]^2 + c*u[2]^2)
end
function noise_hc_gh(du,u,p,t)
	μ, a, b, c, n = p
    du[1] = abs(n) 
    du[2] = abs(n)
    du[3] = abs(n)
end
Δt = 0.01
n = 1e-8
p = [μ, a, b, c, n]
prob_hc = SDEProblem(heteroclinic_cycle_gh_noise, noise_hc_gh, u0, tspan, p; seed=0)
sol = solve(prob_hc, SOSRA(), saveat=0:Δt:T); t = sol.t #nonstiff, didnt test
# sol = solve(prob_hc, SKenCarp(), saveat=0:Δt:T); t = sol.t #stiff, worked well but slow

fig = Figure()
ax1 = Axis(fig[1, 1])
lines!(ax1, t, sol[1,:])
lines!(ax1, t, sol[2,:])
lines!(ax1, t, sol[3,:])

# fig2 = Figure()
ax2 = Axis3(fig[1,2])
# scatter(sol[1,:], sol[2,:], sol[3,:])
lines!(ax2, sol[1,:], sol[2,:], sol[3,:], color=t)

scatter!(ax2, fps[:,1], fps[:,2], fps[:,3], color=:red, markersize=5000)
scatter!(ax2, [0], [0], [0], color=:blue, markersize=8000)





#--------------------------------------------3 COUPLED HODGKIN-HUXLEY NEURONS----------------------------------------------------------------
using GLMakie,  DifferentialEquations

@inbounds function syncoup(gi, S, V, Vrev)
    Isyn = 0.0
    for j = 1:length(gi)
        Isyn += gi[j] * S * (V[j] - Vrev) 
    end
    return Isyn
end


@inbounds function hodgkinhuxley_coeffs(V)
        αn = 0.032 * (-50 - V)/(exp(-(V+50)/5) - 1)
        αm = 0.32 * (-52 -V)/(exp(-(V+52)/4) - 1)
        αh = 0.128 * exp(-(V+48)/18.0 )
        βn = 0.5 * exp(-(V+55)/40.0)
        βm = 0.28 * (25+V) / (exp( (V+25)/5.0 ) -1 )
        βh = 4.0/(1 + exp(-(V+25)/5))
        return αn, αm, αh, βn, βm, βh
end

heaviside(x) = x > 0 ? 1 : 0
mutable struct params
    gs :: Matrix{Float64}
    I :: Float64
    Vna :: Float64
    Vk :: Float64
    Vl :: Float64
    gna :: Float64
    gk :: Float64
    gl :: Float64
    C :: Float64
    Vrev :: Float64
    κ :: Float64
    Smax :: Float64
    τ :: Float64
    Vth :: Float64
end

@inbounds function hodgkinhuxley_rule!(du, u, p, t)
        I = p.I; Vna = p.Vna; Vk=p.Vk; Vl=p.Vl; gna=p.gna; gk=p.gk; gl=p.gl; C=p.C; Vrev=p.Vrev; κ=p.κ; Smax=p.Smax; τ=p.τ; Vth=p.Vth; gs = p.gs;
        N = size(gs, 1)    
        for i=1:N
            V, n, m, h, S, R = @view u[(i-1)*6 + 1 : i*6]
            Vs = @view u[((1:N) .-1) .*6 .+ 1] 
            αn, αm, αh, βn, βm, βh = hodgkinhuxley_coeffs(V)
            gi = @view gs[i,:]
            Isyn = syncoup(gi, S, Vs, Vrev)
            # if i == 1 println("$(Isyn), $(S), $(R), $(V), $(n^4*gk*(V-Vk))") end
            # Isyn = 0.0
            du[(i-1)*6 + 1] = (-I -n^4*gk*(V-Vk) - m^3*h*gna*(V-Vna) - gl*(V-Vl) -Isyn)/C
            du[(i-1)*6 + 2] = αn*(1-n) - βn*n
            du[(i-1)*6 + 3] = αm*(1-m) - βm*m
            du[(i-1)*6 + 4] = αh*(1-h) - βh*h
            du[(i-1)*6 + 5] = ( (R - κ*S)*(Smax - S)/Smax ) / τ
            du[(i-1)*6 + 6] = ( heaviside(V - Vth) - R ) / τ
        end
    return nothing
end

# u = du = rand(7*3)
# t = 0.0
# @btime hodgkinhuxley_rule!(du, u, p, t)
# @code_warntype hodgkinhuxley_rule!(du, u, p, t)
# @btime syncoup(gs[1,:], 1, rand(3), 1)
# @code_warntype syncoup(gs[1,:], 1, rand(3), 1)
# @btime hodgkinhuxley_coeffs(-50.0) 

C = 0.143
gl = 0.02672
Vl = -63.563
gna = 7.15 
Vna = 50.
gk = 1.43 
Vk = -95.
Vrev = -80.
κ = 1/2
Smax = 0.045
τ = 50.
Vth = -20.
N = 3
# g = 50.0
# g = 10.0
g = 0.5
gs = [0 g g; g 0 g; g g 0]
I = 0.08
# I = -5.0 #nice spikes

p = params(gs, I, Vna, Vk, Vl, gna, gk, gl, C, Vrev, κ, Smax, τ, Vth)
T = 100
Ttr = 0
u0 = rand(3*6);

tspan = (0, T);
hcgh = ODEProblem(hodgkinhuxley_rule!, u0, tspan, p);
# sol = solve(hcgh, Rosenbrock23()); t = sol.t;
sol = solve(hcgh, Rodas5()); t = sol.t;
# sol = solve(hcgh, Tsit5()); t = sol.t;
# sol = solve(hcgh, AutoTsit5(Rosenbrock23()), saveat=0:0.01:T, abstol=1e-8, reltol=1e-8); t = sol.t;

t = sol.t; tr=sol[:,:];
Vs = @view tr[((1:N) .-1) .*6 .+ 1, :] 

fig = Figure()
ax = Axis(fig[1,1])
for i=1:N lines!(ax, t, Vs[i,:]) end 