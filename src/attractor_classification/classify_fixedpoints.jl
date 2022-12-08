""" for each point in the trajectory, assign to closest fixed point if within neigh_th; otherwise classify as transition point
fps[:, j] needs to be the j-th fixed point
traj[:,j] needs to be the j-th point in the trajectory (at time t[j]) ie cols contains the state space points
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
