
include("$(srcdir())/visualizations/plot_tools.jl")

function color_trajectory_hc(traj_state_idxs_all; c1=:purple, c2=:green, c3=:orange, fpalpha=1.0, trajcolor=:blue)
    fp_colors = [(c1, fpalpha), (c2, fpalpha), (c3, fpalpha)];
    traj_colors = replace(traj_state_idxs_all, (0:3 .=> [trajcolor; fp_colors])...); #get colors for trajectory
    return fp_colors, traj_colors
end

function prepare_ticks_twotick_plot(xs, ys)
    return [floor(Int64, xs[1]), ceil(Int64, xs[end])], [floor(minimum(ys)), ceil(maximum(ys))]
end

function plot_heterocliniccycle_paper(ts, xs, ys, zs, p, dwelltimes, fps; azim=7, elev=0.39, c1 = :green, c2 = :purple, c3=:red, fig=nothing, idxrow=1, idxcol=1, markersize=12, neigh_th=1e-3, viewmode=:stretch)
    
    traj_state_idxs_all = classify_points_into_fps([xs ys zs]', fps[1:2:end, :]; neigh_th);
    fp_colors, traj_colors = color_trajectory_hc(traj_state_idxs_all; c1, c2, c3, trajcolor=:black)
    
    
    if isnothing(fig) fig = Figure(resolution=(columnsize_pt, 1.0*width_pt)) end
    ax1 = Axis(fig[idxrow, idxcol], ylabel="x", xlabel="t")
    for i=1:1 lines!(ax1, ts, xs, label=["x", "y", "z"][i], color=traj_colors, linewidth=1.4) end
    # xticks, yticks = prepare_ticks_twotick_plot(ts, xs)
    # ax1.xticks = xticks
    # xticks!(ax1, xticks)
    # ax1.yticks = yticks
    # ylims!(ax1, yticks...)

    ax2 = Axis3(fig[idxrow+1, idxcol]; azimuth=azim, elevation=elev, viewmode)
    lines!(ax2, xs, ys, zs, color=traj_colors, linewidth=0.5)
    scatter!(ax2, fps[1:2:end, :][:,1], fps[1:2:end, :][:, 2], fps[1:2:end, :][:, 3], color=[c1, c2, c3], markersize=6)

    ax3 = Axis(fig[idxrow+2, idxcol], ylabel="τ_i", xlabel="i", yscale=log10)
    for i=1:3
        scatter!(ax3, dwelltimes[i], color=fp_colors[i], markersize=3)
    end
    yticks = get_ticks_in_powers(dwelltimes[1])
    @show yticks
    ax3.yticks = yticks
    ylims!(ax3, yticks[1]...)

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
function plot_noisy_bistable(ts, xs, xdots, distribution_info; idxrow = 1, idxcol = 1, fig=nothing, c1 = :green, c2 = :purple)
    colors = [el > 0 ? c1 : c2 for el in xs] #right side (above) receives c1

    if isnothing(fig) fig = Figure(resolution=(columnsize_pt, 1.0*width_pt)) end
    ax1 = Axis(fig[idxrow, idxcol], ylabel="x", xlabel="t")
    lines!(ax1, ts, xs, color=colors, linewidth=0.7)
    # xticks, yticks = prepare_ticks_twotick_plot(ts, xs)
    # ax1.xticks = xticks
    # ax1.yticks = yticks

    ax2 = Axis(fig[idxrow+1, idxcol], ylabel="dx/dt", xlabel="x", yticks=WilkinsonTicks(3))
    scatter!(ax2, xs, xdots;  color=[el > 0 ? c1 : c2 for el in xs], markersize=1.0)

    ax3 = Axis(fig[idxrow+2, idxcol], ylabel="PDF(τ)", xlabel="τ", yscale=log10)
    for (idx, (side, info)) in enumerate(distribution_info)
        @unpack bins, weights, xfit, yfit, A, B = info
        @show weights, yfit, xfit
        scatter!(ax3, bins, weights, color=color=[c2, c1][idx], markersize=4)
        lines!(ax3, xfit, yfit, color=[c2, c1][idx], label="$(trunc(A, sigdigits=2)) exp($(trunc(B, sigdigits=3)) τ)")
        yticks = get_ticks_in_powers(weights)
        @show yticks
        ax3.yticks = yticks
        ylims!(ax3, yticks[1]...)
    end
    
    # axislegend(ax3, position=:rt; framevisible=false, labelsize=10, orientation = :horizontal, margin=(-25,-20,-20,-8))
    return fig, [ax1, ax2, ax3]
end


function plot_attractormerging_crisis(ts, xs, ys, dwelltimes_info; fig=nothing, idxrow=1, idxcol=1, cleft=:black, cright=:red, numbins=10)
    if isnothing(fig) fig = Figure(resolution=(columnsize_pt, 1.0*width_pt)) end
    
    idx_states_time_plot = [x >= 0 ? 1 : 0 for x in xs]
    colors_ax1 = replace(idx_states_time_plot, 1=>cright, 0=>cleft)
    alpha =0.8
    colors_ax2 = replace(idx_states_time_plot, 1=>(cright, alpha), 0=>(cleft, alpha))

    ax1 = Axis(fig[idxrow, idxcol], ylabel="x", xlabel="t")
    lines!(ax1, ts, xs, color=colors_ax1, linewidth=0.7)
    ax1.xticks = [10000, 11000]

    ax2 = Axis(fig[idxrow+1, idxcol], ylabel="dx/dt", xlabel="x")
    scatter!(ax2, xs, ys, color=colors_ax2, markersize=1.0)
    
    ax3 = Axis(fig[idxrow+2, idxcol], ylabel="PDF(τ)", xlabel="τ", yscale=log10)
    for (i, (key, dwelltimes)) in enumerate(dwelltimes_info)
        weights, bins = histogram(dwelltimes, collect(range(0, maximum(dwelltimes); length=numbins)))
        if minimum(weights) <= 0 weights .+=1e-8 end
        scatterlines!(ax3, bins, weights, color=[cleft, cright][i], markersize=4, linewidth=0.7)
        ax3.yticks = get_ticks_in_powers(weights)
    end
    
    return fig, [ax1, ax2, ax3]
end

function plot_typeI_intermittency(ts, xs, ys, zs, dwelltimes, colors; fig=nothing, idxrow=1, idxcol=1, c1, c2, azimuth=5.47, elevation=0.05, numbins=10, viewmode=:stretch)
    
    ax1 = Axis(fig[idxrow, idxcol], ylabel="x", xlabel="t")
    lines!(ax1, ts, xs, color=colors, linewidth=0.6)
    xlims!(ax1, 790, 848)
    ax1.xticks = 800:20:840

    ax2 = Axis3(fig[idxrow+1, idxcol]; azimuth, elevation, viewmode)#, xlabeloffset=5, ylabeloffset=5, zlabeloffset=5, protrusions=5, viewmode=:stretch)
    # hidedecorations!(ax2, label=false, grid=false);
    scatter!(ax2, xs, ys, zs; markersize=1.0, color=colors) 

    ax3 = Axis(fig[idxrow+2, idxcol], ylabel="PDF(τ)", xlabel="τ", yscale=log10)
    for (i, (withinlc, dwelltimes_one)) in enumerate(dwelltimes)
        bins = collect(range(0, 1500, length=numbins))
        weights, _bins = histogram(dwelltimes_one, bins)
        weights  .+= 1e-6
        scatterlines!(ax3, _bins, weights, color=[(c2, 0.8), (c1, 0.8)][i], markersize=4, linewidth=0.7)
        ax3.yticks = get_ticks_in_powers(weights)
    end
    # rowsize!(fig.layout, 2, Relative(0.6))
    # ax3.yticks = [1e-3, 1e-2]

    return fig, [ax1, ax2, ax3]
end

function plot_chaotic_saddle(ts, xs, ys, dwelltimes, colors; fig=nothing, idxrow=1, idxcol=1, color_pdf=:green)

    ax1 = Axis(fig[idxrow, idxcol], ylabel="x", xlabel="t")
    lines!(ax1, ts, xs, color=colors, linewidth=0.7);
    ax1.xticks = 0:100000:200000;
    ax1.yticks = 0:3:6;
    

    ax2 = Axis(fig[idxrow+1, idxcol], ylabel="y", xlabel="x")
    scatter!(ax2, xs, ys; markersize=0.7, color=colors);

    ax3 = Axis(fig[idxrow+2, idxcol], ylabel="PDF(τ)", xlabel="τ", yscale=log10)
    numbins = 15; 
    weights, bins = distribution_times_chaotic_saddle(dwelltimes, numbins);
    scatterlines!(ax3, bins, weights, color=color_pdf, markersize=4, linewidth=0.7);
    ax3.yticks = get_ticks_in_powers(weights)
    # l=lines!(ax3, xfit, A .* exp.(B .* xfit), color=c1)
    # ax3.xticks = 0:4e5:8e5;
    # ax3.yticks = [1e-5, 1e-4]
    # axislegend(ax3, position=:lb; framevisible=false, labelsize=8,orientation = :horizontal, margin=(-8,-8,-8,-8))

    return fig, [ax1, ax2, ax3]
end