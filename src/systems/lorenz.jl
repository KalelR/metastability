
@inbounds @inline function lorenz!(du, u, p, t)
    σ, ρ, β = p
    du[1] = σ*(u[2]-u[1])
    du[2] = u[1]*(ρ-u[3]) - u[2]
    du[3] = u[1]*u[2] - β*u[3]
end

function dwelltimes_info_laminar_phase(d::Dict)
    @unpack T, p_chaos, system_func, u0, p_lc = d
    return dwelltimes_info_laminar_phase(p_chaos, p_lc, T; system_func, u0)
end

function dwelltimes_info_laminar_phase(p_chaos, p_lc, T = 1e5; system_func=lorenz!, u0=[0.1, 0.1, 0.1])
    #get dwell times and colors
    Ttr = 1e3; Δt = 1.0;
    prob = ODEProblem(system_func, u0, (0, T), p_chaos)
    sol = solve(prob, Tsit5(), saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, maxiters=1e9);

    #get the limit cycle (before the bif) and save
    T = Ttr + 50; Δt = 0.05;
    prob = ODEProblem(system_func, u0, (0, T), p_lc); 
    sol_lc = solve(prob, Tsit5(), saveat=Ttr:Δt:T, abstol=1e-8, reltol=1e-8, maxiters=1e9);

    arewithin = mapslices(x->withinset(x, sol_lc', 2), sol[:, :], dims=1)[1, :];
    dwelltimes = length_samevalues_allowfluctuations(arewithin, 3)[1] #1 corresponds to time in  the limit cycle (so 1=laminar)
    distribution_info = Dict("dwelltimes"=>dwelltimes, "sol_lc"=>sol_lc)
end