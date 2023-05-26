function F(x, ϵ, α)
    if x ≤ 0 return 0.0 end
    return exp(-ϵ/x) * x^α
end

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

# --------------------------- Naive implementation --------------------------- #
@inbounds function ratemodel_rule!(du, u, p, t)
    Smax = p.Smax; x₀ = p.x₀; I = p.I; gs = p.gs; τ = p.τ; ϵ = p.ϵ; α = p.α
    N = size(gs, 1)
    for i=1:N
        s, r = @view u[(i-1)*2 + 1 : i*2]
        ss = @view u[((1:N) .-1) .*2 .+ 1]
        gi = @view gs[i,:]
        du[(i-1)*2 + 1] = (1/τ) * (r - s/2) * (Smax - s)/Smax
        du[(i-1)*2 + 2] = x₀ * F(I - ratecoup(gi, ss), ϵ, α) - r/τ
    end
return nothing
end

@inbounds function ratecoup(gi, ss)
    Isyn = 0.0
    for j = 1:length(ss)
        Isyn += gi[j] * ss[j]
    end
    return Isyn
end




# --------------------------- Z = log(x) version --------------------------- #
z_to_s(z, Smax) = Smax - exp(z);
@inbounds function ratecoup_Z(gi, Zs, Smax)
    Isyn = 0.0
    for j = 1:length(Zs)
        Isyn += gi[j] * (Smax - exp(Zs[j]))
    end
    return Isyn
end
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

function fixedpoints_ratemodel(p)
    Smax = p.Smax
    r1 = 0.00866; r2 = 0.03733; r3 = 0.0;
    fp1 = [0 r3 Smax r2 Smax r1];
    fp2 = [Smax r2 Smax r1 0 r3];
    fp3 = [Smax r1 0 r3 Smax r2];
    fps = [fp1; fp2; fp3]';
end


@inbounds function ratemodel_noise_rule!(du,u,p,t)
    @unpack η = p
    # du .= η
    du[[1,3,5]] .= η
end


using RandomNumbers
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

include("$(srcdir())/attractor_classification/classify_fixedpoints.jl")
function dwelltimes_heteroclinicycle(xs, ys, zs, fps; neigh_th = 0.001)
    traj_state_idxs_all = classify_points_into_fps([xs ys zs]', fps[1:2:end, :]; neigh_th);
    dwelltimes, _ = length_samevalues_allowfluctuations(traj_state_idxs_all); #dwell times in each saddle
    return traj_state_idxs_all, dwelltimes
end


