using DrWatson
@quickactivate "metastability"
using CairoMakie, OrdinaryDiffEq
include("$(scriptsdir())/utils.jl")

# ----------------------------------------- Logistic map has attractor-merging crisis when attractors merge (duh) haha ---------------------------------------- #

# Cobweb interactive
# the second range is a convenience for intermittency example of logistic
using InteractiveDynamics

n = 2000
Ttr = 2000
rrange = 3.4:0.0001:4.0; #interesting part of the diagram
lo = Systems.logistic(0.4; r = rrange[1]);
# interactive_cobweb(lo, rrange, 5)
output = orbitdiagram(lo, 1, 1, rrange; n, Ttr);

L = length(rrange);
x = Vector{Float64}(undef, n*L);
y = copy(x);
for j in 1:L
    x[(1 + (j-1)*n):j*n] .= rrange[j]
    y[(1 + (j-1)*n):j*n] .= output[j]
end

fig, ax = scatter(x, y; axis = (xlabel = L"r", ylabel = L"x"),
    markersize = 0.8, color = ("black", 0.05),
)


##plot iterates of logistic map starting at 0.5
n = 20; Ttr=  0
output = orbitdiagram(lo, 1, 1, rrange; n, Ttr, u0=0.5);

L = length(rrange);
x = Vector{Float64}(undef, n*L);
y = copy(x);
for j in 1:L
    x[(1 + (j-1)*n):j*n] .= rrange[j]
    y[(1 + (j-1)*n):j*n] .= output[j]
end

scatter!(ax, x, y; markersize = 3, color = ("red", 0.3))

# ------------------------------------------------- Example on particle in double well (eg duffing). Following Ishii, 1986 (Physics Letters A) ------------------------------------------------ #
include("$(srcdir())/systems/duffing.jl")
include("$(srcdir())/paperplottheme.jl")
include("$(srcdir())/visualizations/plots.jl")
function plot_time_series()
T = 1000
Ttr = 0
Δt=0.1
forcing_vals = [0.849, 0.853]

parameters_ishii1986 = Dict(
    :a => 100, #b in their paper
    :b => 10, #-a in the paper
    :c => 0.0, #not present there (symmetric wells)
    :d => 1.0, # γ in their paper
    :f => 0.0, #F
    :ω => 3.5, #Ω
)

# Time series
diffeq = (alg = Tsit5(), abstol=1e-8, reltol=1e-8)
u0_left = [-0.11, -0.11]
u0_right = [0.11, 0.11]
labels = [["left", "right"], ["left"]]
u0s = replace.(labels, "left"=>u0_left, "right"=>u0_right)
c1 = :purple
c2 = :green
rowtitles = ["Before attractor-merging", "After attractor-merging"]
fig, axs = subplotgrid(2, 1; resolution=(2.0*columnsize_pt, 1.0*columnsize_pt), rowtitles, xlabels="t", ylabels="x", sharex=true, sharey=true)
for (i,f) ∈ enumerate(forcing_vals)
    ax = axs[i]
    df = duffing_assymetric(parameters_ishii1986)
    for (j,u0) in enumerate(u0s[i])
        tr = trajectory(df, T, u0; Ttr, Δt, diffeq); t=Ttr:Δt:Ttr+T
        x = tr[:,1]
        colors = [el > 0 ? c1 : c2 for el in x]
        lines!(ax, t, x, color=colors, label=labels[i][j])
    end
	ylims!(ax, -0.6, 0.6)
end
    # Legend(fig[3,1], axs[1]; tellwidth = false, tellheight = true, framevisible=false, labelsize=12, titlesize=14, orientation=:horizontal, linewidth=3)
    # rowsize!(fig.layout, 3, 0.2)
save("$(plotsdir())/mechanisms/attractormergingcrisis/attractormergingcrisis-doublewellforcing.png", fig, px_per_unit=4)
# save("$(plotsdir())/mechanisms/attractormergingcrisis/attractormergingcrisis-doublewellforcing.pdf", fig)
fig
end

# plot_time_series()


# function bifurcation_diagram(ds, T, varying_param_vals; varying_param_idx = 1, orbit_var_idx = 1, u0 = nothing, varying_param_label = "param", diffeq = NamedTuple(), Ttr = 0.0, Δt=1.0, variable_name = "orbit")
#     orbit_vals = [Float64[] for _ in varying_param_vals]
#     for (i, param_val) in enumerate(varying_param_vals)
#         ds.p[varying_param_idx] = param_val
#         @show ds.p
#         if isnothing(u0) u0 = get_state(ds) end
#         tr = trajectory(ds, T, u0; Ttr, Δt, diffeq);
#         orbit_vals[i] = unique(tr[:, orbit_var_idx])
#     end
#     xs, ys = prepare_flattened_arrays(varying_param_vals, orbit_vals)
#     # bifurcation_diagram!(xs, ys, )
#     fig = Figure()
#     ax = Axis(fig[1,1], ylabel=variable_name, xlabel=varying_param_label)
#     scatter!(ax, xs, ys, markersize=0.8, color=:black)
#     return fig, ax
# end

function bifurcation_diagram(ds, varying_param_vals, mapper; varying_param_idx = 1, orbit_var_idx = 1, varying_param_label = "param", diffeq = NamedTuple(), variable_name = "orbit")
    info_extraction = A -> A[:, orbit_var_idx]
    continuation = RecurrencesSeedingContinuation(mapper; info_extraction);
    grid = mapper.grid
    sampler, = statespace_sampler(Random.MersenneTwister(1234);
        min_bounds = minimum.(grid), max_bounds = maximum.(grid)
    );
    fractions_curves, attractors_info = basins_fractions_continuation(
        continuation, varying_param_vals, varying_param_idx, sampler;
        show_progress = true, samples_per_parameter = 5
    );
    ukeys = unique_keys(attractors_info)
    fig = Figure(resolution=(2.0*columnsize_pt, 1.0*columnsize_pt))
    ax = Axis(fig[1,1], ylabel=variable_name, xlabel=varying_param_label)
    for key in ukeys
        # vals = attractors_info[key]
        vals = [get(attractor_info, key, Float64[]) for attractor_info in attractors_info]
        xs, ys = prepare_flattened_arrays(varying_param_vals, vals)
        scatter!(ax, xs, ys, markersize=0.5, color=(:black, 0.9))
    end

    # bifurcation_diagram!(xs, ys, )
    return fig, ax
end


using Attractors, Random
let
    parameters_ishii1986 = Dict(
    :a => 100, #b in their paper
    :b => 10, #-a in the paper
    :c => 0.0, #not present there (symmetric wells)
    :d => 1.0, # γ in their paper
    :f => 0.0, #F
    :ω => 3.5, #Ω
    )
    df = duffing_assymetric(parameters_ishii1986)
    forcing_vals = range(0, 1.0, length=200)
    varying_param_vals = forcing_vals

    xg = range(-1.1, 1.1, length = 150); grid = ntuple(x->xg, 2);
    mapper = Attractors.AttractorsViaRecurrences(df, grid;
            Δt= 1.0,
            mx_chk_safety = 1e6,
            mx_chk_fnd_att = 200,
            mx_chk_hit_bas = 200,
            mx_chk_loc_att = 10000,
            # mx_chk_att = 10,
            stop_at_Δt = true,
            # diffeq,
        );
    fig, ax = bifurcation_diagram(df, varying_param_vals, mapper; varying_param_idx = 5, varying_param_label = "forcing amplitude")
    ax.xticks = 0:0.1:1.0
    save("$(plotsdir())/mechanisms/attractormergingcrisis/orbitdiagram-doublewellforcing.png", fig, px_per_unit=4)
    fig
end
