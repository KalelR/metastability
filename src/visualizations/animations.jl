function dwelltimes_heteroclinicycle(xs, ys, zs, fps; neigh_th = 0.001)
    traj_state_idxs_all = classify_points_into_fps([xs ys zs]', fps[1:2:end, :]; neigh_th);
    dwelltimes, _ = length_samevalues_allowfluctuations(traj_state_idxs_all); #dwell times in each saddle
    return traj_state_idxs_all, dwelltimes
end

function color_trajectory_hc(traj_state_idxs_all; c1=:purple, c2=:green, c3=:orange, fpalpha=1.0, trajcolor=:blue)
    fp_colors = [(c1, fpalpha), (c2, fpalpha), (c3, fpalpha)];
    traj_colors = replace(traj_state_idxs_all, (0:3 .=> [trajcolor; fp_colors])...); #get colors for trajectory
    return fp_colors, traj_colors
end

_fontsize = 16; _spinewidth = 1.8; _labelsize=20
animationtheme = Theme(
    fontsize = _fontsize,
    Axis=(
        spinewidth=_spinewidth, rightspinevisible=false, topspinevisible=false, xtickalign=1, ytickalign=1, ygridvisible=false, xgridvisible=false,
        xlabelsize = _labelsize, ylabelsize = _labelsize
        ),
    Axis3=(
        xspinewidth=_spinewidth-0.5, yspinewidth=_spinewidth-0.5, zspinewidth=_spinewidth-0.5,
        xtickalign=1, ytickalign=1,
        xlabelsize = _labelsize, ylabelsize = _labelsize, zlabelsize = _labelsize
        ),
    Lines = (linewidth=2, ),
)


function animate_heterolinic_cycle(filename, ts, xs, ys, zs, fps; framerate=100, numframes=100, alphalines = 1.0, mksize=12,
    az=0.78, el = 0.31, c1 = :green, c2 = :purple, c3=:red, fulltrajcolor=:black)
    xyz_points = Observable(Point3f[(xs[1], ys[1], zs[1])])
    tx_points = Observable(Point2f[(ts[1], xs[1])])
    ty_points = Observable(Point2f[(ts[1], ys[1])])
    tz_points = Observable(Point2f[(ts[1], zs[1])])
    tanim = Observable(ts[1])

    traj_state_idxs_all, dwelltimes = dwelltimes_heteroclinicycle(xs, ys, zs, fps; neigh_th=0.001)
    fp_colors,  traj_colors = color_trajectory_hc(traj_state_idxs_all; c1, c2, c3, trajcolor=fulltrajcolor)
# traj_colors = :black

    fig = Figure(resolution=(800, 600)); axs = [];
    ax = Axis3(fig[1:2,1:3], azimuth = az, elevation = el, title= @lift("t = $($tanim)")); push!(axs, ax);

    lines!(ax, xs, ys, zs, color=fulltrajcolor) #full trajectory, black line
    scatter!(ax, fps[1:2:end, :][:,1], fps[1:2:end, :][:, 2], fps[1:2:end, :][:, 3], color=[c1, c2, c3], markersize=mksize+4)
    scatter!(ax, xyz_points, color=:orange, markersize=mksize, overdraw=true)
    hidedecorations!(ax, ticks=false, label=false, ticklabels=false)


    ax = Axis(fig[3,1], xlabel="t", ylabel="x", xticks = WilkinsonTicks(3)); push!(axs, ax);
    lines!(ax, ts, xs, color=traj_colors)
    scatter!(ax, tx_points, color=(:orange, 1.0), markersize=mksize)
    hidespines!(ax, :t, :r)

    ax = Axis(fig[3,2], xlabel="t", ylabel="y"); push!(axs, ax);
    lines!(ax, ts, ys, color=traj_colors)
    scatter!(ax, ty_points, color=(:orange, 1.0), markersize=mksize)
    hidespines!(ax, :t, :r)

    ax = Axis(fig[3,3], xlabel="t", ylabel="z"); push!(axs, ax);
    lines!(ax, ts, zs, color=traj_colors)
    scatter!(ax, tz_points, color=(:orange, 1.0), markersize=mksize)
    hidespines!(ax, :t, :r)
    linkxaxes!(axs[2:end]...)
    linkyaxes!(axs[2:end]...)
    for i=2:4 axs[i].xticks = [minimum(ts), (maximum(ts)+minimum(ts))/2, maximum(ts)] end


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