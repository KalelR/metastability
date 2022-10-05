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

a=0.5; b=8.0; c=0.0; d = 0.2; f = 0.0; ω = 1.0; n₁ = n₂ = n = 0.18;
T = 1e4; Ttr = 0; Δt = 0.5
u0 = [0.0, 0];
p =  [a, b, c, d, f, ω, n₁, n₂];
prob_duffing = SDEProblem(duffing_assymetric_rule, noise_duffing, u0, (0, T), p; seed=0)
sol = solve(prob_duffing, SKenCarp(), saveat=Ttr:Δt:T, maxiters=1e9, progress=true);

#plot double well
U(x, a, b, c) = (a/4)*x^4 - (b/2) * x^2 + c*x
xs_sol = sol[1,:]
big_extreme = maximum(extrema(xs_sol)); xs = range(-big_extreme, +big_extreme, length=100); Us = U.(xs, a, b, c)
tplot = Int(5000/Δt)

sol_plot = sol[:, 1:tplot]; t_plot = sol.t[1:tplot];

#distribution of times
numbins=15
function distribution_times_noisywell(xs, numbins, side=0) #0 is left, 1 is right
    τs, _  = length_samevalues_allowfluctuations(bistable_laminarperiods(xs_sol), side)
    ws, bins = histogram(τs[1], numbins);
end
τs_l = readdlm("$(datadir())/noisybistable-dwelltimes-left-n_0.18.dat")[:,1]
τs_r = readdlm("$(datadir())/noisybistable-dwelltimes-right-n_0.18.dat")[:,1]
weights_l, bins_l =  histogram(τs_l, numbins);
weights_r, bins_r = histogram(τs_r, numbins);
xfit_l = bins_l; yfit_l = weights_l; A_l, B_l = CurveFit.exp_fit(xfit_l, yfit_l .+ 1e-7)
xfit_r = bins_r; yfit_r = weights_r; A_r, B_r = CurveFit.exp_fit(xfit_r, yfit_r .+ 1e-7)

c1 = :green; c2 = :purple

fig = Figure(resolution=(columnsize_pt, 1.0*width_pt))
ax1 = Axis(fig[1, 1], ylabel="x", xlabel="t")
lines!(ax1, t_plot, sol_plot[1,:], color=[el > 0 ? c1 : c2 for el in sol_plot[1,:]])
ax2 = Axis(fig[2, 1], ylabel="dx/dt", xlabel="x")
scatter!(ax2, sol_plot[1,:], sol_plot[2,:],  color=[el > 0 ? c1 : c2 for el in sol_plot[1,:]])
ax3 = Axis(fig[3, 1], ylabel="PDF(τ)", xlabel="τ", yscale=log10)
scatter!(ax3, bins_l, weights_l, color=:green)
lines!(ax3, xfit_l, yfit_l, color=:green, label="P(τ) = $(trunc(A_l, sigdigits=2)) exp($(trunc(B_l, sigdigits=3)) τ)")
# scatter!(ax3, weights_r, bins_r, color=:purple, label="right")
# lines!(ax3, xfit_r, yfit_r, color=:purple)
rowsize!(fig.layout, 2, Relative(0.6))
axislegend(ax3, position=:rt; framevisible=false, labelsize=8,orientation = :horizontal, margin=(-8,-8,-8,-8))
# tit = Label(fig[0,1], "Noisy bistable"; textsize=14, tellwidth=false)
save("$(plotsdir())/paper/noisybistable.png", fig, px_per_unit=4)



T = 1e7
sol = solve(prob_duffing, SKenCarp(), saveat=Ttr:Δt:T, maxiters=1e9, progress=true);
τs, _  = length_samevalues_allowfluctuations(bistable_laminarperiods(sol[1,:]), 0)
writedlm("$(datadir())/noisybistable-dwelltimes-left-n_$n.dat", τs[1])
τs, _  = length_samevalues_allowfluctuations(bistable_laminarperiods(sol[1,:]), 1)
writedlm("$(datadir())/noisybistable-dwelltimes-right-n_$n.dat", τs[1])