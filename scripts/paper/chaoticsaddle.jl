using DrWatson
@quickactivate "metastability"
using DynamicalSystems, Statistics, CurveFit, DelimitedFiles
using GLMakie
include("$(scriptsdir())/utils.jl")
include("$(srcdir())/paperplottheme.jl")

b = 0.9; c = 0.4; d = 6.0; a = 1.003;
T = 2e5; Ttr = 100;
u0 = [-2.0, -2.0]
ik = Systems.ikedamap(u0;a, b, c, d)
traj = trajectory(ik, T; Ttr); t = Ttr:Ttr+T;


numbins = 15; dist_th = 1; icsperdim=50;
# fp = [ 2.9715737943410305, 4.153134755539537 ];
fp = [2.97, 4.15];
#getting distribution of times in chaotic saddle for various ics
filename_dwelltimes = "$(datadir())/durationinchaoticsaddle-d_0.2.dat";
weights, bins = distribution_times_chaotic_saddle(filename_dwelltimes, numbins; T=1e6);
xfit = bins; yfit = weights; A, B = CurveFit.exp_fit(xfit, yfit);


#getting colors that distinguish chaotic saddle from fixed point
bool_within_fp = map(x->iswithinneighborhood(x, fp, dist_th), traj);
colors = [el == true ? :purple : :green for el in bool_within_fp];

fig = Figure(resolution=(columnsize_pt, 1.0*width_pt))
ax1 = Axis(fig[1,1], xlabel="t", ylabel="x");
lines!(ax1, t, traj[:,1], color=colors);
ax1.xticks = 0:100000:200000;
ax1.yticks = 0:3:6;
ax2 = Axis(fig[2:4,1], xlabel="x", ylabel="y");
scatter!(ax2, traj[1:end,1], traj[1:end,2], markersize=1, color=colors);
ax3 = Axis(fig[5, 1], ylabel="PDF(τ)", xlabel="τ", yscale=log10);
scatterlines!(ax3, bins, weights, color=:green, markersize=4);
l=lines!(ax3, xfit, A .* exp.(B .* xfit), color=:red, label="P(τ) = $(trunc(A, sigdigits=2)) exp($(trunc(B, sigdigits=3)) τ)" );
ax3.xticks = 0:4e5:8e5;
axislegend(ax3, position=:lb; framevisible=false, labelsize=8,orientation = :horizontal, margin=(-8,-8,-8,-8))
# save("$(plotsdir())/paper/chaoticsaddle.png", fig, px_per_unit=4)
save("$(plotsdir())/paper/chaoticsaddle.pdf", fig)



# Generating dwell times