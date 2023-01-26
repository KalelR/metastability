
include("$(srcdir())/paperplottheme.jl")
include("$(srcdir())/visualizations/plot_tools.jl")
function plot_heterocliniccycle_paper(t, ss_plot, traj_colors, dwelltimes, fp_colors; azim=7, elev=0.39)
    zs = @view sol[((1:N) .-1) .*2 .+ 1, :]; ss = map(x->z_to_s(x, Smax), zs);
    tplot = 35000; ss_plot = ss[:, 1:tplot]; t_plot = sol.t[1:tplot];
    fps = fixedpoints_ratemodel(p)
    traj_state_idxs_all, dwelltimes = dwelltimes_heteroclinicycle(ss[1,:], ss[2,:], ss[3,:], fps; neigh_th=0.001)
    fp_colors, traj_colors = color_trajectory_hc(traj_state_idxs_all[1:tplot])
    
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

function plot_noisybistable(sol; T=1e4, obtaindwelltime="readfromfile", numbins = 15)
    if obtaindwelltime == "readfromfile"
        τs_l = readdlm("$(datadir())/noisybistable-dwelltimes-left-n_0.18-T_1.0e7.dat")[:,1];
        # τs_r = readdlm("$(datadir())/noisybistable-dwelltimes-right-n_0.18-T_1e.0e8.dat")[:,1];
    else
        τs_l, _  = length_samevalues_allowfluctuations(bistable_laminarperiods(sol[1,:]), 0)
        τs_l = τs_l[1];
        @show τs_l
        τs_r, _  = length_samevalues_allowfluctuations(bistable_laminarperiods(sol[1,:]), 1)
        τs_r = τs_r[1];
        # writedlm("$(datadir())/noisybistable-dwelltimes-left-n_$n-T_$T.dat", τs_l)
        # writedlm("$(datadir())/noisybistable-dwelltimes-right-n_$n-T_$T.dat", τs_r)
    end

    weights_l, bins_l =  histogram(τs_l, numbins);
    xfit_l = bins_l; yfit_l = weights_l; A_l, B_l = CurveFit.exp_fit(xfit_l, yfit_l .+ 1e-7);
    yfit_l = A_l .* exp.(B_l .* xfit_l)
    yfit_l .+= 1e-8; #small value so log10 doesnt break
    weights_l .+= 1e-8; #small value so log10 doesnt break
    println("Exponents of exp fit are $A_l and $B_l")

    c1 = :green; c2 = :purple

    xs_sol = sol[1,:];  big_extreme = maximum(extrema(xs_sol)); xs = range(-big_extreme, +big_extreme, length=100); Us = U.(xs, a, b, c)
    tplot = Int(5000/Δt);  sol_plot = sol[:, 1:tplot]; t_plot = sol.t[1:tplot];

    fig = Figure(resolution=(columnsize_pt, 1.0*width_pt))
    ax1 = Axis(fig[1, 1], ylabel="x", xlabel="t")
    lines!(ax1, t_plot, sol_plot[1,:], color=[el > 0 ? c1 : c2 for el in sol_plot[1,:]])

    ax2 = Axis(fig[2, 1], ylabel="dx/dt", xlabel="x")
    scatter!(ax2, sol_plot[1,:], sol_plot[2,:],  color=[el > 0 ? c1 : c2 for el in sol_plot[1,:]])

    ax3 = Axis(fig[3, 1], ylabel="PDF(τ)", xlabel="τ", yscale=log10)
    scatter!(ax3, bins_l, weights_l, color=:green)
    lines!(ax3, xfit_l, yfit_l, color=:green, label="$(trunc(A_l, sigdigits=2)) exp($(trunc(B_l, sigdigits=3)) τ)")
    ylims!(ax3, 1e-8, 1e-1)
    ax3.yticks = ([1e-8, 1e-5, 1e-2], ["1e-8", "1e-5", "1e-2"])

    rowsize!(fig.layout, 2, Relative(0.6))
    axislegend(ax3, position=:rt; framevisible=false, labelsize=10, orientation = :horizontal, margin=(-25,-20,-20,-8))
    return fig, [ax1, ax2, ax3]
end