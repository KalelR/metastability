
f_circle(x, ω) = x - ω*cos(2π*x) + ω 

function indifferent_circlemap(u, p, t)
    ω = p[1]
    return mod(f_circle(u, ω), 1)
end

using DynamicalSystems
u0 = 0.7
ω = 0.5
ds = DiscreteDynamicalSystem(indifferent_circlemap, u0, [ω])
tr = trajectory(ds, 50.0); ts = 0.0:50.0;
using CairoMakie 
fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, ts, tr)
display(fig)
ax = Axis(fig[1,2])
xs = 0:0.01:0.99
lines!(ax, xs, indifferent_circlemap.(xs, [ω], 0))
lines!(ax, xs, xs, linestyle=:dash, color=(:black, 0.5))
display(fig)



# ------------------------------- Coupled maps ------------------------------- #
function coupling_circlemap(u, i)
    iprev = (i == 1) ? length(u) : i-1
    inext = (i == length(u)) ? 1 : i+1
    return sin(2π*u[iprev]) + sin(2π*u[inext]) - 2*sin(2π*u[i])
end 
function coupled_indifferent_circlemap!(du, u, p, t)
    ω, ϵ = p; 
    for i=1:length(u)
        # du[i] = mod(f_circle(u[i], ω) + ϵ*coupling_circlemap(u, i),1)
        du[i] = mod(f_circle(u[i], ω),1) + ϵ*coupling_circlemap(u, i)
    end
    return du
end

N = 5
u0 = rand(Float64, N)
ω = 0.5
# ϵ = 0.088858; Ttr = 0; Ttr_plot = 160000; T = 5000; #Fig 4; can replicate
ϵ = 0.08885; Ttr=40000; Ttr_plot = 90000; T = 80; #Fig 3; cannot
ds = DiscreteDynamicalSystem(coupled_indifferent_circlemap!, u0, [ω, ϵ])
tr = trajectory(ds, T+Ttr_plot; Ttr); ts = Ttr:T+Ttr+Ttr_plot;
using CairoMakie 
fig = Figure()
ax = Axis(fig[1,1:2])
# lines!(ax, ts[Ttr:end], tr[Ttr:end,1])
for i=1:N lines!(ax, ts[Ttr_plot:end], tr[Ttr_plot:end,i]) end
display(fig)
ax2 = Axis(fig[2,1])
scatter!(ax2, tr[:,1], tr[:,2], markersize=1)
display(fig)
# save("chaoticitinerancy-tsudaspaper-N_$N-ω_$ω-ϵ_$ϵ-u0_$u0.png", fig)

ax = Axis3(fig[2, 2])
scatter!(ax, tr[:,1], tr[:,2], tr[:,3], color=ts)
# lines!(ax, tr[:,1], tr[:,2], tr[:,3])
display(fig)