using DrWatson
@quickactivate
using OrdinaryDiffEq, GLMakie, Peaks
include("$(srcdir())/paperplottheme.jl")
include("$(scriptsdir())/utils.jl")

@inbounds function hindmarshrose_rule!(du, u, p, t)
    @unpack a,b,c,d,r,s, xr, I = p
    du[1] = u[2] - a*u[1]^3 + b*u[1]^2 -u[3] + I
    du[2] = c - d*u[1]^2 - u[2]
    du[3] = r*(s*(u[1] - xr) - u[3])
end

"""
Given the variables z from Hindmarsh-Rose, find the bursts
as the intervals between z's maxima and minima
"""
function classify_points_in_burst(zs)
    count_min = count_max = 1
    state = 0 #silence
    states = zeros(Bool, length(zs));
    idxs_min = Peaks.findminima(zs)[1];
    idxs_max = Peaks.findmaxima(zs)[1];
    for idx in eachindex(zs)
        if count_min ≤ length(idxs_min)
            if idx == idxs_min[count_min]
            state = 1
            count_min += 1;
            end
        end
        if count_max ≤ length(idxs_max)
            if idx == idxs_max[count_max]
                state = 0
                count_max += 1;
            end
        end
        states[idx] = state
    end
return states
end

function construct_delta_dist(x; offset, numpts)
    xs = range(x-offset, x+offset, length=numpts); dist = zeros(length(xs)); dist[11] = 1.0;
    return xs, dist
end

function integrate_and_plot_bursting()
    #integrate
    u0 = [-1.0, 0, 0];
    a=1; b=3; c=1; d=5; xr=-8/5; s=4; r=0.001; I=2.0;
    p = @strdict a b c d xr s r I;
    T = 3000; Ttr=1000; Δt = 0.1;
    prob = ODEProblem(hindmarshrose_rule!, u0, (0.0, T), p);
    sol = solve(prob, AutoTsit5(Rosenbrock23()); saveat=Ttr:Δt:T, maxiters=1e9);
    c1 = :green; c2 = :purple;

    #identify states
    states = classify_points_in_burst(sol[3,:]);
    colors = [el == 1 ? c1 : c2 for el in states];
    dwelltimes, _ = length_samevalues_allowfluctuations(states);
    dwelltimes_lc = dwelltimes[1] .*  Δt
    dwelltimes_node = dwelltimes[0] .*  Δt
    period_lc = dwelltimes_lc[1]; period_node = dwelltimes_node[1];
    #construct dist
    xs_lc, dist_lc = construct_delta_dist(period_lc, offset=100, numpts=21)
    xs_node, dist_node = construct_delta_dist(period_node, offset=100, numpts=21)

    fig = Figure(resolution=(columnsize_pt, 1.0*width_pt))
    ax1 = Axis(fig[1,1], ylabel="V", xlabel="t")
    lines!(ax1, sol.t,  sol[1,:], color=colors);
    xlims!(ax1, high=2900)
    ax1.xticks = [1000, 2000]

    azi = 7.635530633326996; ele = 0.3526990816987241;
    ax2 = Axis3(fig[2,1], xlabel="V", ylabel="y", zlabel="z", xlabeloffset=20, ylabeloffset=20, zlabeloffset=30,
                azimuth=azi, elevation=ele, protrusions=0, viewmode=:stretch)
    l=lines!(ax2, sol[1,:], sol[2,:], sol[3, :], color=colors, colormap=:viridis);
    ax2.zticks = [1.8, 2.1]

    ax3 = Axis(fig[3,1], ylabel="PDF(τ)", xlabel="τ")
    lines!(ax3, xs_lc, dist_lc, color=c1)
    lines!(ax3, xs_node, dist_node, color=c2)
    ax3.xticks = [150, 250]
    xlims!(ax3, 150, 300)
    rowsize!(fig.layout, 2, Relative(0.6))
    return fig, [ax1, ax2, ax3]
end

fig, axs = integrate_and_plot_bursting()
fig
safesave("$(plotsdir())/mechanisms/paper-version/bursting-hindmarshrose.png", fig)