using DataStructures
using ColorSchemes

logrange(x1, x2; length) = (10^y for y in range(log10(x1), log10(x2), length=length))

"""
for each point in trajectory, round it to number digits (sort of equivalent of getting the eps-neighborhood of each point, eps being 10^-numdigits), count number of occurrences of rounded point and divided by total amount of points to give the measure of each point.
"""
function histmeasure(tr, numdigits)
	v = [[tr[i,1], tr[i,2]] for i=1:length(tr)]
	vround = [round.(el, digits=numdigits) for el in v]
	c = counter(vround)
	measure = [c[elround] for elround in vround ]
	measure ./ length(v)
end
function densitypointsintraj_inbins(tr::Vector, numdigits, statespace_boxes)
	v = tr
	vround = [round.(el, digits=numdigits) for el in v]
	c = counter(vround)
	ρ_bins = [c[x] for x ∈ statespace_boxes];
	return ρ_bins ./ sum(ρ_bins)
end

function densitypointsintraj_inpoints(tr::Vector, numdigits)
	v = tr
	vround = [round.(el, digits=numdigits) for el in v]
	c = counter(vround)
	measure = [c[elround] for elround in vround ]
	return measure ./ length(v)
end


using StatsBase, LinearAlgebra
function histogram(v, numbins)
	bins = range(minimum(v), maximum(v), length=numbins+1)
    a = StatsBase.fit(Histogram, v, bins, closed=:right)
    b = normalize(a, mode=:pdf)
    return b.weights, bins[1:end-1]
end

function histogram(v, bins::Vector; mode=:pdf)
    a = StatsBase.fit(Histogram, v, bins, closed=:right)
    b = normalize(a, mode=mode)
    return b.weights, bins[1:end-1]
end

timederivative(sol) = [sol(t, Val{1}) for t ∈ sol.t]
norm(v) = sum(v.^2)

"""
Return the index i of the "box" to which the value x belongs. The boxes are (-inf, box[1]], (box[2], box[3]], (box[3], box[4]], ...
"""
findbox(x, boxes) = searchsortedfirst.(Ref(boxes), x) #index of value inside (box[i], box[i+1]]

"""
# count the number of adjacent (subsequent) permanences of each point in xs inside a box in boxes
Receive the time-series xs and `boxes`, a vector that defines the bins. Returns, for each
bin, the durations that `xs` has consecutively spent within each bin. Note that boxes is
not a very good name. The bins are actually (-inf, box[1]], (box[2], box[3]], .... Check
the findbox documentation for more.  Boxes_duration is thus a vector of vectors. Each
vector i inside contains the durations that xs spent inside bin i. Also returns a vector
with the boxes appended by another bin. I think this actually does not make much sense.
"""
function durationinboxes(xs, boxes)
    boxes_durations = [Float64[] for i=1:length(boxes)+1]
    #start counting once a new box is reached; so first need to find the first transition
    initial_box_id = findbox(xs[1], boxes)
    secondbox_time_id = 0; box_id = 1; box_dur = 1; prev_box_id = 1;
    for i=2:length(xs)-1
        box_id = findbox(xs[i], boxes)
        if box_id != initial_box_id
            secondbox_time_id = i
            push!(boxes_durations[prev_box_id], box_dur)
            prev_box_id = box_id
            box_dur = 1
            break
        else 
            box_dur +=1 
        end
    end
    if secondbox_time_id == 0 boxes_durations[initial_box_id] = [length(x)]; return boxes_durations end#stayed whole series in initial box;
    for (i, x) in enumerate(xs[secondbox_time_id+1:end])
        box_id = findbox(x, boxes)
        if box_id == prev_box_id
            box_dur += 1 ##assuming time step is fixed (which we can assure with interpolation of trajs; future we extend)
        else #left the prev box
            push!(boxes_durations[prev_box_id], box_dur)
            box_dur = 1
            prev_box_id = box_id
        end
    end
    push!(boxes_durations[prev_box_id], box_dur)
    #create bins w same size as boxes_durs
    boxes_extended = deepcopy(boxes)
    if length(boxes_extended) > 1 push!(boxes_extended, boxes_extended[end] + (boxes_extended[end]-boxes_extended[end-1])) end
    return boxes_durations, boxes_extended
end
# box_durs, boxes_ext = durationinboxes([0, 1.1, 1.1, 2, 5, 1.1, 5, 5, 5, 2.1, 2.1, 1.1, 1.1, 1.1, 1.1, 0,0 ], [0,1,2,3])

function vector_to_colors(v, colormap=ColorSchemes.viridis[v])
	v = (v .- minimum(v)) ./ (maximum(v) .- minimum(v)) .* 255 .+ 1;
	v = round.(Int, v );
	colors = colormap[v];
end

"""
fluctuation_th = 0 allows no repetition
"""
function length_samevalues_allowfluctuations(v, num_fluctuations=0)
    unique_vals = unique(v)
    lens_values = Dict(unique_vals .=> [Int64[] for i=1:length(unique_vals)])
    duration = 1
    curr_val = v[1]
    corrected_v = zeros(length(v)); #corrected_v[1] = curr_val
    for i = 1:length(v)-(num_fluctuations+1)
        corrected_v[i] = curr_val
        if any(curr_val .== v[(i+1):(i+1)+num_fluctuations])
            duration += 1
        else
            push!(lens_values[curr_val], duration)
            duration=1; curr_val = v[i+1]
        end
    end
    corrected_v[end] = curr_val
    if duration > 1 push!(lens_values[curr_val], duration-1) end
    return lens_values, corrected_v
end

#second if is to avoid constant values (eg a period-1 fp instead of period-order fp)
function repetition_every_order(v, order; rtol=0.05, atol=0.0)
    vbool = zeros(Bool, length(v)-order)
    for i = 1:(length(v)-order)
        if (isapprox(v[i], v[i+order]; rtol, atol) && !isapprox(v[i], v[i+1]; rtol, atol) )
            vbool[i] = 1
         end
    end
    return vbool
end

"""
Windows std taken over window size ws and every order points (for normal series, order=1).
Used in saddle-node intermittency for finding the laminar period.
"""
function moving_std(v, ws, order=1)
    mv_average = zeros(length(v)-(ws-1))
    for i=1:length(mv_average)
        mv_average[i] = std(v[i:order:i+ws-1])
    end
    return mv_average
end




function peakperiod(F, freqs)
    _, idxmax = findmax(abs.(F))
    return 1 / abs(freqs[idxmax]), F[idxmax]
end

function getspectrum(signal, ts, Ts)
    newsize = Int(nextprod((2,), length(signal))/2)
    signal = signal[1:newsize];    ts = ts[1:newsize]
    signal .-= mean(signal)
    F = fft(signal) |> fftshift
    freqs = fftfreq(length(ts), 1.0/Ts) |> fftshift
    periods = 1 ./ freqs
    maxperiod = peakfrequency(F, freqs)
    return F, freqs, periods, maxperiod
end

function getspectrum2(signal, ts, Ts)
    newsize = Int(nextprod((2,), length(signal))/2)
    signal = signal[1:newsize];    ts = ts[1:newsize]
    signal .-= mean(signal)
    F = abs2.(rfft(signal))
    freqs = rfftfreq(length(signal))./Ts
    periods = 1 ./ freqs
    maxperiod, Fmax = peakperiod(F, freqs)
    return F, freqs, periods, maxperiod, Fmax
end

# estimate_period(signal.-mean(signal), :mt)


using Distances
"""
Check if 'point` is within `set` by seeing if its minimum distance to the set is below an absolute threshold.
"""
function withinset(point, set, threshold)
    mindist = minimum(mapslices(x->evaluate(Euclidean(), point, x), set, dims=2))
    iswithin = mindist <= threshold ? true : false
end

function withinset(point, set, threshold)
    for setpoint ∈ eachrow(set)
        dist = evaluate(Euclidean(), point, setpoint)
        if dist <= threshold return true end
    end
    return false
end

# withinset([1.0, 1.0], [0.9 0.9; 1.1 1.1; 1.3 1.3], 0.1)

# using ChaosTools:get_deviations, stateeltype, _buffered_qr
# function lyapunovspectrum_convergence(integ, N, Δt::Real, Ttr::Real = 0.0, show_progress=false)
#     if show_progress
#         progress = ProgressMeter.Progress(N; desc = "Lyapunov Spectrum: ", dt = 1.0)
#     end
#     B = copy(get_deviations(integ)) # for use in buffer
#     if Ttr > 0
#         t0 = integ.t
#         while integ.t ≤ t0 + Ttr
#             step!(integ, Δt, true)
#             Q, R = _buffered_qr(B, get_deviations(integ))
#             set_deviations!(integ, Q)
#         end
#     end

#     k = size(get_deviations(integ))[2]
#     T = stateeltype(integ)
#     t0 = integ.t; t = zeros(T, N); t[1] = t0
#     λs = [zeros(T, k) for i in 1:N];
#     tr = [zeros(T, k) for i in 1:N];

#     for i in 2:N
#         step!(integ, Δt, true)
#         Q, R = _buffered_qr(B, get_deviations(integ))
#         for j in 1:k
#             @inbounds λs[i][j] = log(abs(R[j,j]))
#         end
#         t[i] = integ.t
#         tr[i] = integ.u[:,1]
#         set_deviations!(integ, Q)
#         show_progress && ProgressMeter.update!(progress, i)
#     end
#     popfirst!(λs); popfirst!(t); popfirst!(tr)
#     return λs, tr, t
# end

# function lyapunovspectrum_convergence_discrete(integ, N, Δt::Real, Ttr::Real = 0.0, show_progress=false)
#     if show_progress
#         progress = ProgressMeter.Progress(N; desc = "Lyapunov Spectrum: ", dt = 1.0)
#     end
#     B = copy(get_deviations(integ)) # for use in buffer
#     if Ttr > 0
#         t0 = integ.t
#         while integ.t ≤ t0 + Ttr
#             step!(integ, Δt, true)
#             Q, R = _buffered_qr(B, get_deviations(integ))
#             set_deviations!(integ, Q)
#         end
#     end

#     k = size(get_deviations(integ))[2]
#     T = stateeltype(integ)
#     t0 = integ.t; t = zeros(T, N); t[1] = t0
#     λs = [zeros(T, k) for i in 1:N];
#     tr = [zeros(T, k) for i in 1:N];

#     for i in 2:N
#         step!(integ, Δt, true)
#         Q, R = _buffered_qr(B, get_deviations(integ))
#         for j in 1:k
#             @inbounds λs[i][j] = log(abs(R[j,j]))
#         end
#         t[i] = integ.t
#         tr[i] = integ.u
#         set_deviations!(integ, Q)
#         show_progress && ProgressMeter.update!(progress, i)
#     end
#     popfirst!(λs); popfirst!(t); popfirst!(tr)
#     return λs, tr, t
# end



"""
Finds laminar periods by identifying the phases in which the windowed standard deviation
of the f^3 map is small (er than a threshold). In the laminar period, f^3 is a fp,
constant, so std is very small. In chaos, it's nonconstant so higher. Works quite well,
but not perfectly, and needs some fine tuning of the parameters.
"""
function logistic_laminarperiods(tr; win_size=4, std_th)
    std_tr3 = moving_std(tr, 3*win_size, 3);
    laminar_period_bool = std_tr3 .< std_th #1s are laminar, 0s are nonlaminar
end
# # bool_fp_order = repetition_every_order([4, 1, 2, 3, 1,2,3, 1,2, 3, 3, 2, 1, 3, 2, 1], 3)
# a,b=length_samevalues_allowfluctuations(bool_fp_order, 0)
"""
Finds laminar periods by identifying the phases that repeat every order=3 iterations
"""
function logistic_laminarperiods(tr, order; kwargs...)
    bool_fp_order = repetition_every_order(tr, order; kwargs...)
    T, bool_laminarperiods = length_samevalues_allowfluctuations(bool_fp_order, order)
    return T[1], bool_laminarperiods
end

function logistic_estimate_laminarperiod_duration(tr; win_size=4, std_th, num_allowed_fluctuations=3)
    laminar_period_bool = logistic_laminarperiods(tr; win_size, std_th)
    T, corr_lpb = length_samevalues_allowfluctuations(laminar_period_bool, num_allowed_fluctuations)
    return T[1], corr_lpb ##laminar period only for 1
end

function logistic_chaoticperiods(tr, order; kwargs...)
    bool_fp_order = repetition_every_order(tr, order; kwargs...)
    T, bool_laminarperiods = length_samevalues_allowfluctuations(bool_fp_order, order)
    return T[0], bool_laminarperiods
end

# --------------------------------- PLOTTING --------------------------------- #
# using GLMakie
using DynamicalSystems
function plot_RM!(fig, t, tr, ϵ; tsmode="scatterline", recurrencetimes=false, logy=false, Δt=1.0)
    RM = RecurrenceMatrix(tr, ϵ)
    RM = skeletonize(RM)
    axs = []
    ax = Axis(fig[1:2,1], ylabel="j", xlabel="i"); push!(axs, ax)
    xs, ys = coordinates(RM)
    scatter!(ax, xs, ys; color = :black, markersize = 3)
    ax.limits = ((0, size(RM, 1)+1), (0, size(RM, 2)+1));
    ax.aspect = 1
    ax = Axis(fig[1:2,2], ylabel="x", xlabel="t"); push!(axs, ax)
    if tsmode == "scatterline" scatterlines!(ax, t, tr[:,1], color=:black)
    else lines!(ax, t, tr[:,1], color=:black) end
    if recurrencetimes
        times, count = histrecurrencetimes(RM; Δt)
        if logy==true count .+= 1; ax = Axis(fig[3,1:2], ylabel="Count", xlabel="recurrence times", yscale=log10)
        else ax = Axis(fig[3,1:2], ylabel="Count", xlabel="recurrence times") end
        if tsmode == "scatterline" scatterlines!(ax, times, count, color=:black)
        else lines!(ax, times, count, color=:black) end
        hlines!(ax, mean(count))
        push!(axs, ax)
    end
    return axs
end


 function animationdata(sol, Tplot, Δt, Δtanimation)
	speeds = norm.(timederivative(sol))
	tplot = Int64(Tplot/Δt);
	Δt_plot = Int64(Δtanimation/Δt)
    t = sol.t; tr = sol[:,:]'
	t_plot = t[1:Δt_plot:tplot];
	tr_plot = tr[1:Δt_plot:tplot, :];
	speeds_plot = log10.(1 ./ speeds[1:Δt_plot:tplot])
	frames = 2:length(t_plot)
    return t_plot, tr_plot, speeds_plot, frames
 end

 function animationdata_tr(tr, t, Tplot, Δt, Δtanimation)
	tplot = Int64(Tplot/Δt);
	Δt_plot = Int64(Δtanimation/Δt)
	t_plot = t[1:Δt_plot:tplot];
	tr_plot = tr[1:Δt_plot:tplot, :];
	frames = 2:length(t_plot)
    return t_plot, tr_plot, frames
 end


function pointspeed_as_colors(speeds_plot)
	#transform the vector with the info for the colors onto a Int vector going from 1 to 264; this is used to index the colormap (wihch has 264 colors); basically transforming it into an vector of indices
	v = (speeds_plot .- minimum(speeds_plot)) ./ (maximum(speeds_plot) .- minimum(speeds_plot)) .* 255 .+ 1
	v = round.(Int, v)
	colors = ColorSchemes.viridis[v]
	return colors
end


# ---------------------------------------------------------------------------- #
#                                Chaotic Saddle                                #
# ---------------------------------------------------------------------------- #
iswithinneighborhood(point, fp, threshold) = evaluate(Euclidean(), point, fp) ≤ threshold ? true : false
# iswithinneighborhood([3.0, 3.0], fp, 0.1)

using DynamicalSystemsBase:obtain_access, get_state
function trajectory_discrete(integ, t, u0 = nothing;
        Δt::Int = 1, Ttr = 0, a = nothing, diffeq = nothing, fp, threshold
    )
    !isnothing(u0) && reinit!(integ, u0)
    Δt = round(Int, Δt)
    T = eltype(get_state(integ))
    t0 = current_time(integ)
    tvec = (t0+Ttr):Δt:(t0+t+Ttr)
    L = length(tvec)
    T = eltype(get_state(integ))
    X = isnothing(a) ? dimension(integ) : length(a)
    data = Vector{SVector{X, T}}(undef, L)
    Ttr ≠ 0 && step!(integ, Ttr)
    data[1] = obtain_access(get_state(integ), a)
	t_convergence = 0
    for i in 2:L
        step!(integ, Δt)
		if iswithinneighborhood(integ.u, fp, threshold) for j=i:L data[j] = SVector{X, T}(fp) end; t_convergence=i; break end
        data[i] = SVector{X, T}(obtain_access(get_state(integ), a))
    end
    return Dataset(data), t_convergence
end

## just by trajectories entering neighborhood of FP; its NOT VERY ELEGANT! haha but a quick way to do it
#All right for getting dwell times, bad for coloring the saddle + fp in ikeda (transition is assigned to saddle!)
function time_to_converge(integ, t, u0 = nothing;
	Δt::Int = 1, Ttr = 0, a = nothing, diffeq = nothing, fp, threshold
)
	!isnothing(u0) && reinit!(integ, u0)
	Δt = round(Int, Δt)
	t0 = current_time(integ)
	tvec = (t0+Ttr):Δt:(t0+t+Ttr)
	L = length(tvec)
	Ttr ≠ 0 && step!(integ, Ttr)
	t_convergence = 0
	for i in 2:L
		step!(integ, Δt)
		if iswithinneighborhood(integ.u, fp, threshold)  return i; break end
	end
	return L
end

function distribution_times_chaotic_saddle(filename::String, numbins; mode=:pdf, T=1e6)
    τs = readdlm(filename)[:,1]
	filter!(x->x>=10, τs)
	filter!(x->x<T, τs)
	bins = collect(range(minimum(τs), maximum(τs), length=numbins+1))
	weights, bins = histogram(τs, bins; mode=mode);
	return weights, bins
end

function distribution_times_chaotic_saddle(numbins; mode=:pdf, T=1e6, Ttr=1e2, threshold=1, icsperdim=20)
	τs = dwell_times_chaotic_saddle(icsperdim; T, Ttr, threshold)
	filter!(x->x>=10, τs)
	filter!(x->x<T, τs)
	bins = collect(range(minimum(τs), maximum(τs), length=numbins+1))
	weights, bins = histogram(τs, bins; mode=mode);
	return weights, bins
end

function dwell_times_chaotic_saddle(icsperdim; T=1e6, Ttr=1e2, threshold=1)
	xs = range(-0.0, 1.0, length=icsperdim)
	ys = range(-2.0, 0,   length=icsperdim)
	u0s = [[x, y] for x ∈ xs for y ∈ ys]
	τs = zeros(Int, length(u0s))
	for (i, u0) ∈ enumerate(u0s)
		integ = integrator(ik, u0)
		t_conv = time_to_converge(integ, T, u0; Ttr, fp, threshold)
		τs[i] = t_conv
	end
	return τs
end


# ---------------------------------------------------------------------------- #
#                                Noisy bistable                                #
# ---------------------------------------------------------------------------- #
bistable_laminarperiods(v) = v .> 0 #1 is positive, 0 is negative