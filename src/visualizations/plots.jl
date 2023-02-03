
include("$(srcdir())/paperplottheme.jl")
include("$(srcdir())/visualizations/plot_tools.jl")

function color_trajectory_hc(traj_state_idxs_all; c1=:purple, c2=:green, c3=:orange, fpalpha=1.0, trajcolor=:blue)
    fp_colors = [(c1, fpalpha), (c2, fpalpha), (c3, fpalpha)];
    traj_colors = replace(traj_state_idxs_all, (0:3 .=> [trajcolor; fp_colors])...); #get colors for trajectory
    return fp_colors, traj_colors
end

function plot_heterocliniccycle_paper(ts, xs, ys, zs, p, dwelltimes, fps; azim=7, elev=0.39, c1 = :green, c2 = :purple, c3=:red, fig=nothing, idxrow=1, idxcol=1, markersize=12, neigh_th=1e-3)
    
    traj_state_idxs_all = classify_points_into_fps([xs ys zs]', fps[1:2:end, :]; neigh_th);
    fp_colors, traj_colors = color_trajectory_hc(traj_state_idxs_all; c1, c2, c3, trajcolor=:black)
    
    
    if isnothing(fig) fig = Figure(resolution=(columnsize_pt, 1.0*width_pt)) end
    ax1 = Axis(fig[idxrow, idxcol], ylabel="x", xlabel="t")
    for i=1:1 lines!(ax1, ts, xs, label=["x", "y", "z"][i], color=traj_colors, linewidth=4) end

    ax2 = Axis3(fig[idxrow+1, idxcol], azimuth=azim, elevation=elev)
    lines!(ax2, xs, ys, zs, color=traj_colors)
    scatter!(ax2, fps[1:2:end, :][:,1], fps[1:2:end, :][:, 2], fps[1:2:end, :][:, 3], color=[c1, c2, c3], markersize=markersize)

    ax3 = Axis(fig[idxrow+2, idxcol], ylabel="τ_occurrence", xlabel="occurrence", xscale=log10, yscale=log10)
    for i=1:3
        scatter!(ax3, dwelltimes[i], color=fp_colors[i], markersize=6)
    end

    return fig, [ax1, ax2, ax3]
end


function prepare_flattened_arrays(varying_param_vals, values_allparams)
	xs_flat = Float64[]; ys_flat = Float64[];
	for (i, (param, vals)) in enumerate(zip(varying_param_vals, values_allparams))
		xs = [param for j in vals]
		ys = vals
		push!(xs_flat, xs...)
		push!(ys_flat, ys...)
	end
    return xs_flat, ys_flat
end

function plot_noisy_bistable(ts, xs, xdots, distribution_info; idxrow = 1, idxcol = 1, fig=nothing, c1 = :green, c2 = :purple, markersize)
    
   # xs_sol = sol[1,:];  #big_extreme = maximum(extrema(xs_sol)); #xs = range(-big_extreme, +big_extreme, length=100); #Us = U.(xs, a, b, c)
    # tplot = Int(5000/Δt);  sol_plot = sol[:, 1:tplot]; t_plot = sol.t[1:tplot];
    colors = [el > 0 ? c1 : c2 for el in xs]

    if isnothing(fig) fig = Figure(resolution=(columnsize_pt, 1.0*width_pt)) end
    ax1 = Axis(fig[idxrow, idxcol], ylabel="x", xlabel="t")
    lines!(ax1, ts, xs, color=colors)


    ax2 = Axis(fig[idxrow+1, idxcol], ylabel="dx/dt", xlabel="x")
    scatter!(ax2, xs, xdots;  color=[el > 0 ? c1 : c2 for el in xs], markersize)

    @unpack bins_l, weights_l, xfit_l, yfit_l, A_l, B_l = distribution_info
    ax3 = Axis(fig[idxrow+2, idxcol], ylabel="PDF(τ)", xlabel="τ", yscale=log10)
    scatter!(ax3, bins_l, weights_l, color=:green)
    lines!(ax3, xfit_l, yfit_l, color=:green, label="$(trunc(A_l, sigdigits=2)) exp($(trunc(B_l, sigdigits=3)) τ)")
    ylims!(ax3, 1e-8, 1e-1)
    ax3.yticks = ([1e-8, 1e-5, 1e-2], ["1e-8", "1e-5", "1e-2"])

    # rowsize!(fig.layout, 2, Relative(0.6))
    # axislegend(ax3, position=:rt; framevisible=false, labelsize=10, orientation = :horizontal, margin=(-25,-20,-20,-8))
    return fig, [ax1, ax2, ax3]
end


function plot_attractormerging_crisis(ts, xs, ys, dwelltimes; fig=nothing, idxrow=1, idxcol=1, cleft=:black, cright=:red, numbins=10)
    if isnothing(fig) fig = Figure(resolution=(columnsize_pt, 1.0*width_pt)) end
    
    idx_states_time_plot = [x >= 0 ? 1 : 0 for x in xs]
    colors_ax1 = replace(idx_states_time_plot, 1=>cright, 0=>cleft)
    alpha =0.8
    colors_ax2 = replace(idx_states_time_plot, 1=>(cright, alpha), 0=>(cleft, alpha))

    ax1 = Axis(fig[idxrow, idxcol], ylabel="x", xlabel="t")
    lines!(ax1, ts, xs, color=colors_ax1, linewidth=0.7)
    ax1.xticks = [10000, 10500, 11000]

    ax2 = Axis(fig[idxrow+1, idxcol], ylabel="dx/dt", xlabel="x")
    # lines!(ax2, xs, ys, color=colors, linewidth=0.5)
    scatter!(ax2, xs, ys, color=colors_ax2, markersize=1.0)
    
    ax3 = Axis(fig[idxrow+2, idxcol], ylabel="PDF(τ)", xlabel="τ", yscale=log10)
    # ax3.xticks = 0:2000:4000
    
    for (i, dwelltimes_one) in enumerate(dwelltimes)
        weights, bins = histogram(dwelltimes_one, collect(range(0, maximum(dwelltimes_one); length=numbins)))
        # weights .+=1e-8
        scatterlines!(ax3, bins, weights, color=[cleft, cright][i])
    end
    
    return fig, [ax1, ax2, ax3]
end

function plot_typeI_intermittency(ts, xs, ys, zs, dwelltimes, colors; fig=nothing, idxrow=1, idxcol=1, c1, c2, azimuth=5.47, elevation=0.05, numbins=10)
    
    ax1 = Axis(fig[idxrow, idxcol], ylabel="x", xlabel="t")
    lines!(ax1, ts, xs, color=colors)
    xlims!(ax1, 790, 848)
    ax1.xticks = 800:20:840

    ax2 = Axis3(fig[idxrow+1, idxcol]; azimuth, elevation)#, xlabeloffset=5, ylabeloffset=5, zlabeloffset=5, protrusions=5, viewmode=:stretch)
    # hidedecorations!(ax2, label=false, grid=false);
    scatter!(ax2, xs, ys, zs; markersize=2, color=colors) 

    ax3 = Axis(fig[idxrow+2, idxcol], ylabel="PDF(τ)", xlabel="τ", yscale=log10)
    for (i, dwelltimes_one) in enumerate(values(dwelltimes))
        bins = collect(range(0, 1500, length=numbins))
        weights, _bins = histogram(dwelltimes_one, bins)
        @show weights
        scatterlines!(ax3, _bins, weights .+ 1e-6, color=[(c1, 0.8), (c2, 0.8)][i])
    end
    # rowsize!(fig.layout, 2, Relative(0.6))
    ax3.yticks = [1e-3, 1e-2]

    return fig, [ax1, ax2, ax3]
end

function plot_chaotic_saddle(ts, xs, ys, dwelltimes, colors; fig=nothing, idxrow=1, idxcol=1)

    ax1 = Axis(fig[idxrow, idxcol], ylabel="x", xlabel="t")
    lines!(ax1, t, xs, color=colors);
    ax1.xticks = 0:100000:200000;
    ax1.yticks = 0:3:6;
    

    ax2 = Axis(fig[idxrow+1, idxcol], ylabel="y", xlabel="x")
    scatter!(ax2, xs, ys; markersize=1, color=colors);

    ax3 = Axis(fig[idxrow+2, idxcol], ylabel="PDF(τ)", xlabel="τ", ytickformat="{1:1.0e}")
    numbins = 15; 
    weights, bins = distribution_times_chaotic_saddle(dwelltimes, numbins; T);
    scatterlines!(ax3, bins, weights, color=:green, markersize=4);
    # l=lines!(ax3, xfit, A .* exp.(B .* xfit), color=c1)
    # ax3.xticks = 0:4e5:8e5;
    # ax3.yticks = [1e-5, 1e-4]
    # axislegend(ax3, position=:lb; framevisible=false, labelsize=8,orientation = :horizontal, margin=(-8,-8,-8,-8))

    return fig, [ax1, ax2, ax3]
end