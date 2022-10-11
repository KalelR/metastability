using DrWatson
@quickactivate "metastability"
using GLMakie, OrdinaryDiffEq, CurveFit
include("$(scriptsdir())/utils.jl")
include("$(srcdir())/paperplottheme.jl")

# -------------------------------- Rate model -------------------------------- #
function F(x, ϵ, α)
    if x ≤ 0 return 0.0 end
    return exp(-ϵ/x) * x^α
end

@inbounds function ratecoup_Z(gi, Zs, Smax)
    Isyn = 0.0
    for j = 1:length(Zs)
        Isyn += gi[j] * (Smax - exp(Zs[j]))
    end
    return Isyn
end

@inbounds function ratemodel_rule_Z!(du, u, p, t)
    Smax = p.Smax; x₀ = p.x₀; I = p.I; gs = p.gs; τ = p.τ; ϵ = p.ϵ; α = p.α
    N = size(gs, 1)
    for i=1:N
        Z, r = @view u[(i-1)*2 + 1 : i*2]
        Zs = @view u[((1:N) .-1) .*2 .+ 1]
        gi = @view gs[i,:]
        du[(i-1)*2 + 1] = (1/(τ*Smax)) * ( ((Smax - exp(Z))/2) - r)
        du[(i-1)*2 + 2] = x₀ * F(I - ratecoup_Z(gi, Zs, Smax), ϵ, α) - r/τ
    end
return nothing
end

z_to_s(z, Smax) = Smax - exp(z);

mutable struct rateparams
    gs :: Matrix{Float64}
    I :: Float64
    ϵ :: Float64
    Smax :: Float64
    τ :: Float64
    x₀ :: Float64
    α :: Float64
end

""" for each point in the trajectory, assign to closest fixed point if within neigh_th; otherwise classify as transition point
fps[:, j] needs to be the j-th fixed point
traj[:,j] needs to be the j-th point in the trajectory (at time t[j])
"""
function classify_points_into_fps(traj, fps; neigh_th=0.5)
    traj_state_idxs = zeros(Int, size(traj, 2))
    for (i, pt) in enumerate(eachcol(traj) )
        state = 0 #transition
        for j=1:3
            if iswithinneighborhood(pt, fps[:, j], neigh_th)  state = j end
        end
        traj_state_idxs[i] = state
    end
    return traj_state_idxs
end

"""utility: create matrix traj to be used in classify points"""
function classify_points_into_fps_from_sol(sol, fps, tplot=0; neigh_th)
    if tplot == 0 tplot = size(sol, 2) end
    sol_a = sol[:, 1:tplot];
    zs_a = @view sol_a[((1:N) .-1) .*2 .+ 1, :]; ss_a = map(x->z_to_s(x, Smax), zs_a); rs_a = @view sol_a[((1:N) .-1) .*2 .+ 2, :];
    traj = zeros(Float64, (size(sol_a))); traj[1:2:end, :] .= ss_a; traj[2:2:end, :] .= rs_a;
    traj_state_idxs = classify_points_into_fps(traj, fps; neigh_th = 0.01);
end


#Model parameters
N = 3; τ = 50.; ϵ = 1e-3; I = 0.145; Smax = 0.045; g₁ = 3.0; g₂ = 0.7; x₀ = 2.57e-3; α = 0.564;
gs = zeros(Float64, (3,3));
gs[2,1] = gs[3,2] = gs[1,3] = g₁;
gs[1,2] = gs[2,3] = gs[3,1] = g₂;
gs[1,1] = gs[2,2] = gs[3,3] = 0.;
p = rateparams(gs, I, ϵ, Smax, τ, x₀, α);
#Sim params
T = 1e6; Ttr = 1000; Δt = 1.0;
# u0 = rand(3*2);
u0 = [0.5, 0.2, 0.4, 0.9, 0.5, 0.6];
tspan = (0, T);
hcgh = ODEProblem(ratemodel_rule_Z!, u0, tspan, p);
sol = solve(hcgh, Vern9(), maxiters=1e9, saveat=Ttr:Δt:T); ts = sol.t;
zs = @view sol[((1:N) .-1) .*2 .+ 1, :]; ss = map(x->z_to_s(x, Smax), zs);
tplot = 35000;
ss_plot = ss[:, 1:tplot]; t_plot = sol.t[1:tplot];


#system's fps
r1 = 0.00866; r2 = 0.03733; r3 = 0.0;
fp1 = [0 r3 Smax r2 Smax r1];
fp2 = [Smax r2 Smax r1 0 r3];
fp3 = [Smax r1 0 r3 Smax r2];

fps = [fp1; fp2; fp3]';
azim = 6.995530633326985; elev = 0.3926990816987241;
fp_colors = [:purple, :green, :orange];

traj_state_idxs_all = classify_points_into_fps_from_sol(sol, fps, 0; neigh_th = 0.01);
dwelltimes, _ = length_samevalues_allowfluctuations(traj_state_idxs_all); #dwell times in each saddle
traj_colors = replace(traj_state_idxs_all[1:tplot], (0:3 .=> [:blue; fp_colors])...); #get colors for trajectory


fig = Figure(resolution=(columnsize_pt, 1.0*width_pt))
ax1 = Axis(fig[1, 1], ylabel="x", xlabel="t")
for i=1:1 lines!(ax1, t_plot, ss_plot[i,:], label=["x", "y", "z"][i], color=traj_colors) end
ax2 = Axis3(fig[2, 1], azimuth=azim, elevation=elev, xlabeloffset=20, ylabeloffset=20, zlabeloffset=40, protrusions=0, xticklabelpad=2, viewmode=:stretch)
lines!(ax2, ss_plot[1,:], ss_plot[2,:], ss_plot[3,:], color=traj_colors)
scatter!(ax2, fps[1:2:end, :][:,1], fps[1:2:end, :][:, 2], fps[1:2:end, :][:, 3], color=fp_colors, markersize=10)
hidedecorations!(ax2, label=false, grid=false);
# hidedecorations!(ax2);
# hidespines!(ax2);
ax3 = Axis(fig[3, 1], ylabel="τ_occurrence", xlabel="occurrence", xscale=log10, yscale=log10)
for i=1:3
     scatter!(ax3, dwelltimes[i], color=fp_colors[i], markersize=4)
    # xfit = 1:length(dwelltimes[i]); yfit = dwelltimes[i]; A, B = CurveFit.linear_fit(xfit, yfit)
    # lines!(ax3, 10^A .+ B.*xfit, color=fp_colors[i])
    # xfit = 1:length(dwelltimes[i]); yfit = dwelltimes[i]; A, B = CurveFit.power_fit(xfit, yfit)
    # A = trunc(A, sigdigits=2); B = trunc(B, sigdigits=2);
    # lines!(ax3, A .* xfit .^ B, color=fp_colors[i], label="y = $A * x^$B")
    end
# axislegend(ax3, position=:rb; framevisible=false, labelsize=8,orientation = :vertical, margin=(-8,-8,-8,-8), rowgap=-10)
rowsize!(fig.layout, 2, Relative(0.6))
# rowgap!(fig.layout, 0)
ax1.xticks = 0:20000:tplot;
ax2.xticks = 0:0.02:0.04; ax2.yticks = 0:0.02:0.04; ax2.zticks = 0:0.02:0.04;
# save("$(plotsdir())/paper/heterocliniccycle.png", fig, px_per_unit=4)
save("$(plotsdir())/paper/heterocliniccycle.pdf", fig)

# ------------------------- Guckenheimer holmes model ------------------------ #
# (not working nicely, prob the same transofrmation as in the rate model is needed)
#=
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
=#