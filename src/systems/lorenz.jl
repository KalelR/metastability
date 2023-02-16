
@inbounds function lorenz!(du, u, p, t)
    σ, ρ, β = p
    du[1] = σ*(u[2]-u[1])
    du[2] = u[1]*(ρ-u[3]) - u[2]
    du[3] = u[1]*u[2] - β*u[3]
end

@inbounds function lorenz(u, p, t)
    σ, ρ, β = p
    du1 = σ*(u[2]-u[1])
    du2 = u[1]*(ρ-u[3]) - u[2]
    du3 = u[1]*u[2] - β*u[3]
    return SVector{3}(du1, du2, du3)
end


function dwelltimes_info_laminar_phase(d::Dict)
    @unpack T, p_chaos, system_func, u0, p_lc, Δt_chaos, threshold = d
    return dwelltimes_info_laminar_phase(p_chaos, p_lc, T, Δt_chaos, threshold; system_func, u0)
end

function dwelltimes_info_laminar_phase(p_chaos, p_lc, T = 1e5, Δt_chaos=1.0, threshold=10; system_func=lorenz!, u0=[0.1, 0.1, 0.1])
    #get dwell times and colors
    Ttr = 1e3; 
    prob = ODEProblem(system_func, u0, (0, T), p_chaos)
    sol = solve(prob, Tsit5(), saveat=Ttr:Δt_chaos:T, abstol=1e-8, reltol=1e-8, maxiters=1e9);

    #get the limit cycle (before the bif) and save
    T = Ttr + 3; Δt = 0.0025;
    prob = ODEProblem(system_func, u0, (0, T), p_lc); 
    sol_lc = solve(prob, Tsit5(), saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, maxiters=1e9);

    arewithin = mapslices(x->withinset(x, sol_lc', threshold), sol[:, :], dims=1)[1, :];
    dwelltimes_int = length_samevalues_allowfluctuations(arewithin, 3)[1] #1 corresponds to time in  the limit cycle (so 1=laminar)
    dwelltimes = multiply_dict_number(dwelltimes_int, Δt_chaos)
    distribution_info = Dict("dwelltimes"=>dwelltimes, "sol_lc"=>sol_lc)
end

