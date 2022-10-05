using DrWatson
@quickactivate "metastability"
using GLMakie, OrdinaryDiffEq
include("$(scriptsdir())/utils.jl")
include("$(srcdir())/paperplottheme.jl")



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



T = 10000; Ttr = 100; Δt=1.0
μ = 1.0; a = -0.2; b = -0.1; c = -1 - a - b; p = [μ, a, b, c]
u0 = [0.1, 0.2, 0.3] #off the coordinate plane or axis should be attracted towards HC

xpcoor(μ,p) = sqrt(-μ/p)
fps_x = [xpcoor(μ, a), 0, 0];
fps_y = [0, xpcoor(μ, a), 0];
fps_z = [0, 0, xpcoor(μ, a)];
fps = reduce(hcat, [fps_x, fps_y, fps_z]...)

ff = ODEFunction(heteroclinic_cycle_gh!; jac=jac!)
hcgh = ODEProblem(ff, u0, (0, T), p)
sol = solve(hcgh, RK4(), dt=1/1000, saveat=Ttr:Δt:T); ts = sol.t;

tplot = 500;
sol_plot = sol[:, 1:tplot]; t_plot = sol.t[1:tplot];
fp_colors = [:purple, :green, :orange]
azim = 6.665530633326982
elev = 0.19269908169872405

neigh_th = 0.05
traj_state_idxs = zeros(Int, size(sol, 2))
for (i, pt) in enumerate(eachcol(sol) )
    state = 0 #transition
    for j=1:3
        if iswithinneighborhood(pt, fps[:, j], neigh_th)  state = j end
    end
    traj_state_idxs[i] = state
end
dwelltimes, _ = length_samevalues_allowfluctuations(traj_state_idxs)
traj_colors = replace(traj_state_idxs[1:tplot], (0:3 .=> [:blue; fp_colors])...)

fig = Figure(resolution=(columnsize_pt, 1.0*width_pt))
ax1 = Axis(fig[1, 1], ylabel="x,y,z", xlabel="t")
lines!(ax1, t_plot, sol_plot[1,:], label="x", color=traj_colors)
lines!(ax1, t_plot, sol_plot[2,:], label="y", color=traj_colors)
lines!(ax1, t_plot, sol_plot[3,:], label="z", color=traj_colors)
# axislegend(ax1, position=:lb; framevisible=false, labelsize=8,orientation = :vertical, margin=(-8,-8,-8,-8))
ax2 = Axis3(fig[2, 1], azimuth=azim, elevation=elev, xlabeloffset=20, ylabeloffset=20, zlabeloffset=30)
lines!(ax2, sol_plot[1,:], sol_plot[2,:], sol_plot[3,:], color=traj_colors)
scatter!(ax2, fps[1,:], fps[2,:], fps[3,:], color=fp_colors, markersize=10)
ax3 = Axis(fig[3, 1], ylabel="τ_occurrence", xlabel="occurrence")
for i=1:3 scatter!(dwelltimes[i], color=fp_colors[i]) end
# rowsize!(fig.layout, 2, Relative(1/2))
rowgap!(fig.layout, 5)
save("$(plotsdir())/paper/heterocliniccycle_ghrk4.png", fig, px_per_unit=4)