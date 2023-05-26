using DrWatson 
@quickactivate
using GLMakie, DifferentialEquations, RandomNumbers

"""
Eq. 18 in Ashwin, 2011
"""
@inbounds function ratemodel_rule_Z!(du, u, p, t)
    Smax = p.Smax; x₀ = p.x₀; I = p.I; gs = p.gs; τ = p.τ; ϵ = p.ϵ; α = p.α
    N = size(gs, 1)
    for i=1:N
        Z, r = @view u[(i-1)*2 + 1 : i*2]
        Zs = @view u[((1:N) .-1) .*2 .+ 1]
        gi = @view gs[i,:]
        du[(i-1)*2 + 1] = (1/(τ*Smax)) * ( ((Smax - exp(Z))/2) - r)
        du[(i-1)*2 + 2] = x₀ * F(I - ratecoup_Z(gi, Zs, Smax), ϵ, α) - r/τ
    end
    return nothing
end

function F(x, ϵ, α)
    if x ≤ 0 return 0.0 end
    return exp(-ϵ/x) * x^α
end

@inbounds function ratecoup_Z(gi, Zs, Smax)
    Isyn = 0.0
    for j = 1:length(Zs)
        Isyn += gi[j] * (Smax - exp(Zs[j]))
    end
    return Isyn
end

@inbounds function ratemodel_noise_rule!(du,u,p,t)
    @unpack η = p
    du .= η
    # du[[1,3,5]] .= η #only in first variables
end


z_to_s(z, Smax) = Smax - exp(z);

mutable struct rateparams
    gs :: Matrix{Float64}
    I :: Float64
    ϵ :: Float64
    Smax :: Float64
    τ :: Float64
    x₀ :: Float64
    α :: Float64
    η :: Float64
end

function rateparams()
    N = 3
    τ = 50.
    ϵ = 1e-3
    I = 0.145
    Smax = 0.045
    g₁ = 3.0
    g₂ = 0.7
    x₀ = 2.57e-3
    α = 0.564
    gs = zeros(Float64, (3,3))
    gs[2,1] = gs[3,2] = gs[1,3] = g₁
    gs[1,2] = gs[2,3] = gs[3,1] = g₂
    gs[1,1] = gs[2,2] = gs[3,3] = 0.
    η = 0.0
    return rateparams(gs, I, ϵ, Smax, τ, x₀, α, η)
end

"""
Creates the ODEProblem to solve the heteroclinic cycle model
"""
function heteroclinic_cycle(; η=0.0, u0 = [0.5, 0.2, 0.4, 0.9, 0.5, 0.6], T = 1e6, rngseed=1)
    p = rateparams()
    if η != 0 
        p.η = η
        rng = Xorshifts.Xoroshiro128Plus(rngseed)
        W = WienerProcess(0.0, 0.0, 0.0; rng)
        hcgh = SDEProblem(ratemodel_rule_Z!, ratemodel_noise_rule!, u0, (0, T), p; noise=W, seed=rngseed)
    else
        hcgh = ODEProblem(ratemodel_rule_Z!, u0, (0, T), p);
    end
    return hcgh, p
end

function run_heteroclinic_cycle(; noise_strength=0.0, Ttr=50000, T=100000, Δt=1.0)
    u0 = [-5, 0.031, -5, 0.03, -5, 0.03]
    hcgh, params = heteroclinic_cycle(; η=noise_strength, T, rngseed=1, u0)
    solver = noise_strength == 0.0 ? Vern9() : SKenCarp(); #vern9 is for deterministic, Skencarp is a good one for stochastic
    sol = solve(hcgh, solver, maxiters=1e9, saveat=Ttr:Δt:T, abstol=1e-9, reltol=1e-9);
    N = 3
    logvars = @view sol[((1:N) .- 1) .* 2 .+ 1, :];
    ogvars = map(x->z_to_s(x, params.Smax), logvars);
    return sol.t, ogvars
end

function plot_heteroclinic_cycle(ts, xs, ys, zs; title="", azimuth=7, elevation=0.4)
    fig = Figure(); axs = [];
    
    ax = Axis3(fig[1:3, 1]; azimuth, elevation); push!(axs, ax);
    lines!(ax, xs, ys, zs, color=:black) 
    
    for (idx, variable) in enumerate([xs, ys, zs])
        ax = Axis(fig[idx, 2]); push!(axs, ax)
        lines!(ax, ts, variable; color=:black)
    end
    
    supertitle(fig, title)
    return fig, axs 
end




fig, axs = run_heteroclinic_cycle(noise_strength=1e-4; Ttr=0.0, T=1e5)  
xs, ys, zs = eachrow(vars)
fig, axs = plot_heteroclinic_cycle(ts, xs, ys, zs; title="noise strength = $(noise_strength)")