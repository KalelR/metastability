using InteractiveDynamics, GLMakie, DynamicalSystems

# the second range is a convenience for intermittency example of logistic
# rrange = 1:0.001:4.0
# rrange = (rc = 1 + sqrt(8); [rc, rc - 1e-5, rc - 1e-3])
# rrange = 3.8248:0.00001:3.8249
# rrange = 3.7:0.0001:3.845
rrange = 3.8:0.00001:3.84
lo = Systems.logistic(0.4; r = rrange[1])
interactive_cobweb(lo, rrange, 5)



# Distribution of laminar period
lo = Systems.logistic(0.4; r = 3.8284)
Ttr = 500
T = 500000
# λt, t = ChaosTools.lyapunovspectrum_convergence(lo, T; Ttr)

win_size = 4
std_th = 0.001;

u0s = sort(rand(1000))
Ts = zeros(Float64, length(u0s))
for (idx,u0) ∈ enumerate(u0s) 
    tr = trajectory(lo, T, u0; Ttr); t_tr = Ttr:1:Ttr+T
    Ts[idx] = estimate_laminarperiod_duration(tr; win_size, std_th)
end


tr = trajectory(lo, T, rand(); Ttr); t_tr = Ttr:1:Ttr+T
Ts = estimate_laminarperiod_duration(tr; win_size, std_th)

fig = Figure()
ax = Axis(fig[1,1])
hist!(ax, Ts; bins=40)

function estimate_laminarperiod_duration(tr; win_size=4, std_th)
    std_tr3 = moving_std(tr, 3*win_size);
    laminar_period_bool = std_tr3 .< std_th #0s are laminar, 1s are nonlaminar
    T = length_samevalues_allowfluctuations(laminar_period_bool, 3)[0] #laminar period only for 0
end

function length_samevalues(v)
    unique_vals = unique(v)
    lens_values = Dict(unique_vals .=> [Int64[] for i=1:length(unique_vals)])
    current_value = v[1]
    duration = 1
    for i = 2:length(v)
        if v[i] == current_value 
            duration += 1
        else 
            push!(lens_values[current_value], duration)
            current_value = v[i]; duration=1
        end
    end
    return lens_values 
end
"""
fluctuation_th = 0 allows no repetition
"""
function length_samevalues_allowfluctuations(v, num_fluctuations=0)
    unique_vals = unique(v)
    lens_values = Dict(unique_vals .=> [Int64[] for i=1:length(unique_vals)])
    duration = 0
    for i = 1:length(v)-(1+num_fluctuations)
        if any(v[i] .== v[(i+1):(i+1)+num_fluctuations]) 
            duration += 1
        else 
            push!(lens_values[v[i-1]], duration)
            duration=1
        end
    end
    return lens_values 
end


# length_samevalues([1,1,1,0,0,0,1,1,1,0,1,0,1,1,0,0,0,0])
length_samevalues_allowfluctuations([1,1,1,0,0,0,1,1,1,0,1,0,1,1,0,0,0,0], 1)


## Below: tried a few different ways: through finite-time lyapunovs (tried mean and min, its hard because they fluctuate a lot); and finally with std of f^3. This worked best.


fig = Figure()
ax = Axis(fig[1,1])
lines!(t, λt, color=λt .> 0, colormap=:bluesreds)
ax2 = Axis(fig[2,1])
lines!(t_tr, tr)

using DynamicalSystemsBase:DDS
function lyapunovspectrum_instantaneous(ds::DDS{false, T, 1}, N; Ttr = 0) where {T}
    x = get_state(ds); f = ds.f
    p = ds.p; t0 = ds.t0
    if Ttr > 0
        for i in t0:(Ttr+t0)
            x = f(x, p, i)
        end
    end
    t = (t0+Ttr):(t0+Ttr+N)
    λs = zeros(T, length(t));
    for idx = 1:length(t)
        i = t[idx]
        x = f(x, p, i)
        @inbounds λs[idx] = log(abs(ds.jacobian(x, p, i)))
    end
    t = collect(t)
    return λs, t
end


"""
moving average in window of size ws, walking through each element (so offset is 1)
"""
function moving_average(v, ws)
    mv_average = zeros(length(v)-(ws-1))
    for i=1:length(mv_average)
        mv_average[i] = mean(v[i:i+ws-1])
    end
    return mv_average 
end


function moving_maximum(v, ws)
    mv_average = zeros(length(v)-(ws-1))
    for i=1:length(mv_average)
        mv_average[i] = maximum(v[i:i+ws-1])
    end
    return mv_average 
end

tr = trajectory(lo, T; Ttr); t_tr = Ttr:1:Ttr+T
λt, t = lyapunovspectrum_instantaneous(lo, T; Ttr)
win_size = 12
λt_fs = moving_average(λt, win_size); 
# λt_fs = moving_maximum(λt, win_size); 
t_fs = t[1:end-(win_size-1)]
t_tr = t_tr[1:end-(win_size-1)]
tr = tr[1:end-(win_size-1)]

# moving_average([1,1,1,2,2,2,3,1,1,1,1], 3)
λ_th = 0.18

fig = Figure()
ax = Axis(fig[1,1])
# lines!(t, λt, color=λt .> 0, colormap=:bluesreds)
lines!(ax, t_fs, λt_fs, color=λt_fs .> λ_th, colormap=:bluesreds)
hlines!(ax, λ_th, color=:black, linestyle="--")
ax2 = Axis(fig[2,1])
# scatterlines!(t_tr, tr, color=λt_fs .> 0, colormap=:bluesreds)
# scatterlines!(t_tr, tr, color=λt_fs .> 0, colormap=[:red, :blue])
lines!(ax2, t_tr, tr, color=λt_fs .> λ_th, colormap=[:blue, :red])
linkxaxes!(ax, ax2)

## trying with f3 
function moving_std(v, ws)
    mv_average = zeros(length(v)-(ws-1))
    for i=1:length(mv_average)
        mv_average[i] = std(v[i:3:i+ws-1])
    end
    return mv_average 
end
# tr3 = tr[1:3:end]
# t_tr3 = t[1:3:end]
win_size = 4
tr = trajectory(lo, T; Ttr); t_tr = Ttr:1:Ttr+T;
std_tr3 = moving_std(tr, 3*win_size)
t_stdtr3 = t[1:end-3*win_size+1];
t_tr = t_tr[1:end-(3*win_size-1)];
tr = tr[1:end-(3*win_size-1)];

# std_th = 0.0006; winsize = 2
std_th = 0.001;
fig = Figure()
ax = Axis(fig[1,1])
lines!(ax, t_stdtr3, std_tr3, color=std_tr3 .< std_th, colormap=[:red, :blue])
hlines!(ax, std_th, color=:black, linestyle="--")

ax2 = Axis(fig[2,1])
lines!(ax2, t_tr, tr, color=std_tr3 .< std_th, colormap=[:red, :blue])
linkxaxes!(ax, ax2)




# scaling of laminar period
r = 3.83
x = 0:0.01:1
logistic_derivative = lo.jacobian.(x, r, 1)
lines(x, logistic_derivative)


#for lorenz maybe: ChaosTools.lyapunovspectrum_convergence(tinteg, T, 1, Ttr); tinteg = tangent_integrator(lo, 1)