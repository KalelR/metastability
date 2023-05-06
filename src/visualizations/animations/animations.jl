
# const _FONTSIZE = 16; const _SPINEWIDTH = 1.8; const _LABELSIZE=28; const _TITLESIZE = 20
_FONTSIZE = 16;_SPINEWIDTH = 1.8; _LABELSIZE=28; _TITLESIZE = 20
const _RESOLUTION = (1000, 600)
animationtheme = Theme(
    fontsize = _FONTSIZE,
    Axis=(
        spinewidth=_SPINEWIDTH, rightspinevisible=false, topspinevisible=false, xtickalign=1, ytickalign=1, ygridvisible=false, xgridvisible=false,
        xlabelsize = _LABELSIZE, ylabelsize = _LABELSIZE, titlesize = _TITLESIZE
        ),
    Axis3=(
        xspinewidth=_SPINEWIDTH-0.5, yspinewidth=_SPINEWIDTH-0.5, zspinewidth=_SPINEWIDTH-0.5,
        xtickalign=1, ytickalign=1,
        xlabelsize = _LABELSIZE, ylabelsize = _LABELSIZE, zlabelsize = _LABELSIZE, titlesize = _TITLESIZE
        ),
    Lines = (linewidth=2, ),
    # Figure = (resolution=(1000, 600), ),
)


function animate_heterolinic_cycle(filename, ts, xs, ys, zs, fps; framerate=100, numframes=100, alphalines = 1.0, mksize=12, az=0.78, el = 0.31, ball_color=:orange, c1 = :green, c2 = :purple, c3=:red, fulltrajcolor=:black)
    xyz_points = Observable(Point3f[(xs[1], ys[1], zs[1])])
    tx_points = Observable(Point2f[(ts[1], xs[1])])
    ty_points = Observable(Point2f[(ts[1], ys[1])])
    tz_points = Observable(Point2f[(ts[1], zs[1])])
    tanim = Observable(ts[1])

    traj_state_idxs_all, dwelltimes = dwelltimes_heteroclinicycle(xs, ys, zs, fps; neigh_th=0.001)
    fp_colors,  traj_colors = color_trajectory_hc(traj_state_idxs_all; c1, c2, c3, trajcolor=fulltrajcolor)

    fig = Figure(; resolution); axs = [];
    ax = Axis3(fig[1:2,1:3], azimuth = az, elevation = el, title= @lift("t = $($tanim)")); push!(axs, ax);

    lines!(ax, xs, ys, zs, color=fulltrajcolor) #full trajectory, black line
    scatter!(ax, fps[1:2:end, :][:,1], fps[1:2:end, :][:, 2], fps[1:2:end, :][:, 3], color=[c1, c2, c3], markersize=mksize+4)
    scatter!(ax, xyz_points, color=:orange, markersize=mksize, overdraw=true)
    hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
    ax.xticks = [0.0, 0.02, 0.04]
    ax.yticks = [0.0, 0.02, 0.04]
    ax.zticks = [0.0, 0.02, 0.04]


    timeseries_all = [xs, ys, zs]; ylabel_all = ["x", "y", "z"]; point_position_all = [tx_points, ty_points, tz_points]
    for (idx, idxcol) in enumerate(1:3)
        ylabel = ylabel_all[idx] 
        ax = Axis(fig[3, idxcol]; xlabel="t", ylabel); push!(axs, ax);
        timeseries = timeseries_all[idx]
        lines!(ax, ts, timeseries, color=traj_colors)
        point_position = point_position_all[idx]
        scatter!(ax, point_position; color=(ball_color, 1.0), markersize=mksize)
        hidespines!(ax, :t, :r)
        ax.xticks = [minimum(ts), (maximum(ts)+minimum(ts))/2, maximum(ts)]
        ax.yticks = [0.0, 0.02, 0.04]
    end
    linkxaxes!(axs[2:end]...)
    linkyaxes!(axs[2:end]...)


    mkpath(dirname(filename))
    ts_iterator = floor.(Int, range(1, length(ts), length=numframes))
    record(fig, filename, ts_iterator; framerate) do idx_t
        tanim[] = ts[idx_t]
        new_xyz_point = Point3f(xs[idx_t], ys[idx_t], zs[idx_t])
        new_tx_point = Point2f(ts[idx_t], xs[idx_t])
        new_ty_point = Point2f(ts[idx_t], ys[idx_t])
        new_tz_point = Point2f(ts[idx_t], zs[idx_t])
        xyz_points[] = [new_xyz_point]
        tx_points[] = [new_tx_point]
        ty_points[] = [new_ty_point]
        tz_points[] = [new_tz_point]
    end
end


function animate_double_well(filename, xs, ts, xs_well, double_well::Function; c1 = "#440154FF", c2 = "#FDE725FF", ball_color = "#FF1400", framerate=75, numframes=100, ball_size=20, resolution=(800, 600))
    time = Observable(ts[1])
    colors_double_well = [el > 0 ? c1 : c2 for el in xs_well]
    points_double_well = Observable(Point2f[(xs[1], double_well.(xs[1]))])
    Us_well = double_well.(xs_well)
    colors_time_series = [el > 0 ? c1 : c2 for el in xs]
    points_time_series = Observable(Point2f[(ts[1], xs[1])])
    
    fig = Figure(; resolution); axs = [];
    ax = Axis(fig[1,1]; ylabel = "Potential energy", xlabel = "position", title= @lift("t = $($time)")); push!(axs, ax)
    lines!(ax, xs_well, Us_well; linewidth=4, color=colors_double_well)
    s=scatter!(ax, points_double_well, color=ball_color, markersize=ball_size)


    ax = Axis(fig[2, 1], ylabel="position", xlabel="time"); push!(axs, ax)
    lines!(ax, ts, xs; color=colors_time_series)
    scatter!(ax, points_time_series; color=ball_color, markersize=ball_size)

    ts_iterator = floor.(Int, range(1, length(ts), length=numframes))
    mkpath(dirname(filename))
    record(fig, filename, ts_iterator; framerate) do idx_t
        time[] = ts[idx_t]
        new_dw_point = Point2f(xs[idx_t], double_well(xs[idx_t]))
        new_ts_point = Point2f(ts[idx_t], xs[idx_t])
        points_double_well[] = [new_dw_point]
        points_time_series[] = [new_ts_point]
    end
end

"""
differentiate between xs animatino and xs fixed structure because i just want the animatino to have points close to the transition from saddle to focus, while for the structure in state space ofthe entire saddle i want more points
"""
function animate_chaotic_saddle(filename, xs_animation, ys_animation, ts_animation, xs_fixedstructures, ys_fixedstructures, colors_state_space, markers_state_space, colors_animation; c1 = "#440154FF", c2 = "#FDE725FF", ball_color = "#FF1400", framerate=75, numframes=100, ball_size=20, resolution=(800, 600))
    
    points_state_space = Observable(Point2f[(xs_animation[1], ys_animation[1])])
    points_time_series = Observable(Point2f[(ts_animation[1], xs_animation[1])])
    tanim = Observable(ts_animation[1])

    fig = Figure(;resolution)
    ax = Axis(fig[1,1]; title= @lift("t = $($tanim)"), ylabel="y", xlabel="x")

    s=scatter!(ax, xs_fixedstructures, ys_fixedstructures; color=colors_state_space, markersize=markers_state_space)
    s=scatter!(ax, points_state_space; color=ball_color, markersize=ball_size)
    hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
    xlims!(-0.6, 6.2)
    ylims!(-2.8, 6)

    ax = Axis(fig[2, 1]; ylabel="x", xlabel="t")
    scatter!(ax, ts_animation, xs_animation, color=colors_animation, markersize=4)
    scatter!(ax, points_time_series; color=ball_color, markersize=ball_size)

    ts_iterator = floor.(Int, range(1, length(ts_animation), length=numframes))
    mkpath(dirname(filename))
    record(fig, filename, ts_iterator; framerate) do idx_t
        tanim[] = ts_animation[idx_t]
        new_ss_point = Point2f(xs_animation[idx_t], ys_animation[idx_t])
        new_ts_point = Point2f(ts_animation[idx_t], xs_animation[idx_t])
        points_state_space[] = [new_ss_point]
        points_time_series[] = [new_ts_point]
    end

end

function animate_type_I_intermittency(filename, xs, ys, zs, ts, colors_state_space, xs_fixedstructures, ys_fixedstructures, zs_fixedstructures; framerate=100, numframes=100, azimuth = -1.14, elevation = 0.09, resolution=(800, 600), c1=:green, c2=:purple, alpha=0.8, ball_color=:orange, ball_size=15)
    points_state_space = Observable(Point3f[(xs[1], ys[1], zs[1])])
    points_time_series = Observable(Point2f[(ts[1], xs[1])])
    tanim = Observable(ts[1])

    fig = Figure(; resolution)
    ax = Axis3(fig[1,1]; azimuth, elevation, title= @lift("t = $((round($tanim; digits=1)))"))
    lines!(ax, xs_fixedstructures, ys_fixedstructures, zs_fixedstructures, color=colors_state_space)
    scatter!(ax, points_state_space, color=ball_color, markersize=ball_size)
    hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
    limits!(ax, -50, 50, -100, 100, 50, 250)

    ax = Axis(fig[2,1], ylabel="x", xlabel="t")
    lines!(ax, ts, xs, color=colors_state_space)
    scatter!(ax, points_time_series, color=ball_color, markersize=ball_size)
    hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
    hidespines!(ax, :t, :r)
    rowsize!(fig.layout, 1, Relative(0.7))

    ts_iterator = floor.(Int, range(1, length(ts), length=numframes))
    mkpath(dirname(filename))
    record(fig, filename, ts_iterator; framerate) do idx_t
        tanim[] = ts[idx_t]
        new_ss_point = Point3f(xs[idx_t], ys[idx_t], zs[idx_t])
        new_ts_point = Point2f(ts[idx_t], xs[idx_t])
        points_state_space[] = [new_ss_point]
        points_time_series[] = [new_ts_point ]
    end
end


function animate_interior_crisis(filename, xs, ys, ts, colors_state_space_anim, xs_fixedstructures, ys_fixedstructures, colors_state_space, fps_5; framerate=100, numframes=100, azimuth = -1.14, elevation = 0.09, resolution=(800, 600), c1=:green, c2=:purple, alpha=0.8, ball_color=:orange, ball_size=15)

    points_state_space = Observable(Point2f[(xs[1], ys[1])])
    points_time_series = Observable(Point2f[(ts[1], xs[1])])
    tanim = Observable(ts[1])

    fig = Figure(; resolution)
    ax = Axis(fig[1,1], title= @lift("t = $($tanim)"), ylabel="y", xlabel="x")
    s=scatter!(ax, xs_fixedstructures, ys_fixedstructures, color=colors_state_space, markersize=3)
    s=scatter!(ax, points_state_space; color=ball_color, markersize=ball_size)
    scatter!(ax, fps_5[:,1], fps_5[:,2], color=:red, markersize=10)
    hidedecorations!(ax, ticks=false, label=false, ticklabels=false)
    xmin = -0.6; xmax=1.7
    xlims!(ax, xmin, xmax); ylims!(ax, -1.8, 1.4)

    ax2 = Axis(fig[2,1], ylabel="x", xlabel="t")
    lines!(ax2, ts, xs; color=colors_state_space_anim, markersize=5)
    scatter!(ax2, points_time_series; color=ball_color, markersize=ball_size)
    for j=1:5 hlines!(ax2, fps_5[j,1], color=(:red, 0.5), linestyle="--") end
    rowsize!(fig.layout, 1, Relative(0.7))
    
    ts_iterator = floor.(Int, range(1, length(ts), length=numframes))
    mkpath(dirname(filename))
    record(fig, filename, ts_iterator; framerate) do idx_t
        tanim[] = ts[idx_t]
        new_ss_point = Point2f(xs[idx_t], ys[idx_t])
        new_ts_point = Point2f(ts[idx_t], xs[idx_t])
        points_state_space[] = [new_ss_point]
        points_time_series[] = [new_ts_point ]

    end



end