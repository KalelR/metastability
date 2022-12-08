
function plot_heterocliniccycle_paper(t_plot, ss_plot, traj_colors, dwelltimes, fp_colors; azim=7, elev=0.39)
    fig = Figure(resolution=(columnsize_pt, 1.0*width_pt))
    ax1 = Axis(fig[1, 1], ylabel="x", xlabel="t")
    for i=1:1 lines!(ax1, t_plot, ss_plot[i,:], label=["x", "y", "z"][i], color=traj_colors) end

    ax2 = Axis3(fig[2, 1], azimuth=azim, elevation=elev, xlabeloffset=30, ylabeloffset=30, zlabeloffset=50, protrusions=0,  viewmode=:stretch)
    lines!(ax2, ss_plot[1,:], ss_plot[2,:], ss_plot[3,:], color=traj_colors)
    scatter!(ax2, fps[1:2:end, :][:,1], fps[1:2:end, :][:, 2], fps[1:2:end, :][:, 3], color=fp_colors, markersize=10)
    # hidedecorations!(ax2, label=false, grid=false, ticks=false);

    ax3 = Axis(fig[3, 1], ylabel="τ_occurrence", xlabel="occurrence", xscale=log10, yscale=log10)
    for i=1:3
        scatter!(ax3, dwelltimes[i], color=fp_colors[i], markersize=4)
    end

    rowsize!(fig.layout, 2, Relative(0.6))
    return fig, [ax1, ax2, ax3]
end
