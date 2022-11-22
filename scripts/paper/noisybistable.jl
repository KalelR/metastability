using DrWatson
@quickactivate "metastability"
using GLMakie
using DifferentialEquations, CurveFit, DelimitedFiles
include("$(srcdir())/paperplottheme.jl")
include("$(scriptsdir())/utils.jl")

# """
# U(x) = ax^4/4 - bx 2/2 + cx (c: assymetry paramter); gamma: dissipation, f forcing omega forcing freq, n noise strength
# """
@inbounds function duffing_assymetric_rule(du, u, p, t)
    a,b,c,d,f, ω, n₁, n₂ = p
    x, v = u
    du[1] = v
    du[2] = (-a*x^3 + b*x - c) - d*v + f*cos(ω*t)
end

function noise_duffing(du,u,p,t)
    a,b,c,d,f, ω, n₁, n₂ = p
    du[1] = n₁
    du[2] = n₂
end

function distribution_times_noisywell(xs, numbins, side=0) #0 is left, 1 is right
    τs, _  = length_samevalues_allowfluctuations(bistable_laminarperiods(xs_sol), side)
    ws, bins = histogram(τs[1], numbins);
end

U(x, a, b, c) = (a/4)*x^4 - (b/2) * x^2 + c*x #double well

function  integrate_and_plot_noisybistable(; T=1e4, obtaindwelltime="readfromfile")
    a=0.5; b=8.0; c=0.0; d = 0.2; f = 0.0; ω = 1.0; n₁ = n₂ = n = 0.18;
    Ttr = 0; Δt = 0.5
    u0 = [0.0, 0];
    p =  [a, b, c, d, f, ω, n₁, n₂];
    prob_duffing = SDEProblem(duffing_assymetric_rule, noise_duffing, u0, (0, T), p; seed=0)
    sol = solve(prob_duffing, SKenCarp(), saveat=Ttr:Δt:T, maxiters=1e9, progress=true);

    numbins = 15
    if obtaindwelltime == "readfromfile"
        τs_l = readdlm("$(datadir())/noisybistable-dwelltimes-left-n_0.18-T_10000.0.dat")[:,1];
        # τs_r = readdlm("$(datadir())/noisybistable-dwelltimes-right-n_0.18.dat")[:,1];
    else
        τs_l, _  = length_samevalues_allowfluctuations(bistable_laminarperiods(sol[1,:]), 0)
        τs_l = τs_l[1];
        @show τs_l
        τs_r, _  = length_samevalues_allowfluctuations(bistable_laminarperiods(sol[1,:]), 1)
        τs_r = τs_r[1];
        writedlm("$(datadir())/noisybistable-dwelltimes-left-n_$n-T_$T.dat", τs_l)
        writedlm("$(datadir())/noisybistable-dwelltimes-right-n_$n-T_$T.dat", τs_r)
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
# scatter!(ax3, weights_r, bins_r, color=:purple, label="right")
# lines!(ax3, xfit_r, yfit_r, color=:purple)

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

fig, axs = integrate_and_plot_noisybistable(T=1e7; obtaindwelltime="generate"); #from the paper, takes a while to integrate this long
# fig, axs = integrate_and_plot_noisybistable(T=1e4; obtaindwelltime="readfromfile"); #quick version: integrates for a short time and reads dwell times from long integration made previously
# safesave("$(plotsdir())/mechanisms/paper-version/noisybistable.png", fig)
fig