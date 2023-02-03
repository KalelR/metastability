reduced_r(r, rc) = abs(r-rc)/rc
reduced_to_normal(μ, rc) = rc*(1-μ)
rc = 1+√8
# savedir="typeI"

function logistic_get_colors_trajectory!(traj, t)
    Ts, _= logistic_laminarperiods(traj, 3; atol=0.0, rtol=0.02);
    _, corr_lam_periods_bool= logistic_laminarperiods(traj, 3; atol=0.0, rtol=0.02); #for colors
    colors = [el > 0 ? :green : :purple for el in corr_lam_periods_bool];
    t = t[1:length(corr_lam_periods_bool)]; 
    traj = traj[1:length(corr_lam_periods_bool)];
    return colors, traj, t 
end

