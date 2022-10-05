using DrWatson
@quickactivate
using DifferentialEquations
# using CairoMakie
using GLMakie
include("$(srcdir())/paperplottheme.jl")
using Peaks


include("$(scriptsdir())/utils.jl")
@inbounds function hindmarshrose_rule!(du, u, p, t)
    @unpack a,b,c,d,r,s, xr, I = p
    du[1] = u[2] - a*u[1]^3 + b*u[1]^2 -u[3] + I
    du[2] = c - d*u[1]^2 - u[2]
    du[3] = r*(s*(u[1] - xr) - u[3])
end

function classify_points_in_burst(zs)
    count_min = count_max = 1
    state = 0 #silence
    states = zeros(Bool, length(zs));
    idxs_min = Peaks.minima(zs);
    idxs_max = Peaks.maxima(zs);
    for idx in eachindex(zs)
        if count_min ≤ length(idxs_min)
            if idx == idxs_min[count_min]
            state = 1
            count_min += 1;
            # if count_min > length(idxs_min) break end
            end
        end
        if count_max ≤ length(idxs_max)
            if idx == idxs_max[count_max]
                state = 0
                count_max += 1;
            # if count_max > length(idxs_max) states[idx:end] .= 0; break end
            end
        end
        states[idx] = state
    end
return states
end

u0 = [-1.0, 0, 0];
a=1; b=3; c=1; d=5; xr=-8/5; s=4; r=0.001; I=2.0;
p = @strdict a b c d xr s r I
T = 3000; Ttr=1000; Δt = 0.1
prob = ODEProblem(hindmarshrose_rule!, u0, (0.0, T), p)
sol = solve(prob, AutoTsit5(Rosenbrock23()); saveat=Ttr:Δt:T, maxiters=1e9)



azi = 4.735530633327
ele = 0.41269908169872416

# fig = Figure(resolution=(columnsize_pt, 0.50*width_pt), figure_padding=20)
# fig = Figure(resolution=(columnsize_pt, 0.50*width_pt))
# gab = fig[1:3,1] = GridLayout()
# ga = gab[1,1] = GridLayout()
# gb = gab[2:3, 1] = GridLayout()

# ax1 = Axis(ga[1,1], ylabel="V", xlabel="t", title="Bursting")
# lines!(ax1, sol.t,  sol[1,:], color=speeds);
# xlims!(ax1, high=2900)
# ax2 = Axis3(gb[1,1:2], xlabel="V", ylabel="y", zlabel="z", xlabeloffset=20, ylabeloffset=20, zlabeloffset=30, azimuth=azi, elevation=ele)
# l=lines!(ax2, sol[1,:], sol[2,:], sol[3, :], color=speeds, colormap=:viridis);
# cb=Colorbar(gb[1,3], l, label="normalized speed", tellheight=true)
# colgap!(gb, 0)
# cb.alignmode = Mixed(right = 0)
# cb.width=10
# rowgap!(gab, -20)
# fig
# resize_to_layout!(fig)
# save("$(plotsdir())/paper/bursting-hindmarshrose.png", fig, px_per_unit=3)


states = classify_points_in_burst(sol[3,:])
colors = [el == 1 ? :green : :purple for el in states]

fig = Figure(resolution=(columnsize_pt, 1.0*width_pt))
ax1 = Axis(fig[1,1], ylabel="V", xlabel="t", title="Bursting")
lines!(ax1, sol.t,  sol[1,:], color=colors);
xlims!(ax1, high=2900)
ax2 = Axis3(fig[2,1], xlabel="V", ylabel="y", zlabel="z", xlabeloffset=20, ylabeloffset=20, zlabeloffset=30, azimuth=azi, elevation=ele, width=0.9*width_pt, tellwidth=false)
hidedecorations!(ax2); hidespines!(ax2);
l=lines!(ax2, sol[1,:], sol[2,:], sol[3, :], color=colors, colormap=:viridis);
# cb.height=10
# cb.valign=:top
# cb.alignmode = Mixed(right=0)
rowsize!(fig.layout, 2, Relative(0.6))
rowgap!(fig.layout, -40)
resize_to_layout!(fig)
save("$(plotsdir())/paper/bursting-hindmarshrose.png", fig, px_per_unit=3)




# TESTING BURST START AND END ALG

fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, sol.t, sol[1,:])
idxs_min = Peaks.minima(sol[3,:]); scatter!(sol.t[idxs_min], sol[1,idxs_min], color=:orange)
idxs_max = Peaks.maxima(sol[3,:]); scatter!(sol.t[idxs_max], sol[1,idxs_max], color=:red)
states = classify_points_in_burst(sol[3,:])
colors = [el == 1 ? :green : :purple for el in states]

speeds = :black
fig = Figure(resolution=(columnsize_pt, 0.50*width_pt))
ax1 = Axis(fig[1, 1], ylabel="V", xlabel="t", title="Bursting")
lines!(ax1, sol.t,  sol[1,:], color=speeds)
ax2 = Axis3(fig[2:3, 1][1,1], xlabel="V", ylabel="y", zlabel="z")
l=lines!(ax2, sol[1,:], sol[2,:], sol[3, :], color=speeds, colormap=:viridis);
scatter!(ax2, [-1.618 for i=1:100], [-12.09 for i=1:100], range(1.7, 2.15, length=100))
# Colorbar(fig[2:3,1][1,2], l, label="normalized speed", tellheight=true)
fig
