using DataStructures
using ColorSchemes
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
	bins = range(minimum(v), maximum(v), length=numbins)
		a = fit(Histogram, v, bins)
		normalize(a, mode=:pdf)
		return a.weights, bins
end


timederivative(sol) = [sol(t, Val{1}) for t ∈ sol.t]
norm(v) = sum(v.^2)

findbox(x, boxes) = searchsortedfirst.(Ref(boxes), x) #index of value inside (box[i], box[i+1]]

"""
count the number of adjacent (subsequent) permanences of each point in xs inside a box in boxes 
"""
function durationinboxes(xs, boxes)
    boxes_durations = [Float64[] for i=1:length(boxes)+1]
    #start counting once a new box is reached; so first need to find the first transition
    initial_box_id = findbox(xs[1], boxes)
    secondbox_time_id = 0; box_id = 1; box_dur = 0; prev_box_id = 1;
    for i=2:length(xs)-1
        box_id = findbox(xs[i], boxes)
        if box_id != initial_box_id 
            secondbox_time_id = i 
            prev_box_id = box_id
            box_dur = 1
            break 
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
    push!(boxes_extended, boxes_extended[end] + (boxes_extended[end]-boxes_extended[end-1]))
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




function peakfrequency(F, freqs)
    _, idxmax = findmax(abs.(F))
    return 1 / abs(freqs[idxmax])
end

function getspectrum(signal, ts, Ts)
    newsize = Int(nextprod((2,), length(signal))/2)
    signal = signal[1:newsize]
    ts = ts[1:newsize]
    signal .-= mean(signal)
    F = fft(signal) |> fftshift
    freqs = fftfreq(length(ts), 1.0/Ts) |> fftshift
    periods = 1 ./ freqs
    maxperiod = peakfrequency(F, freqs)
    return F, freqs, periods, maxperiod
end

# estimate_period(signal.-mean(signal), :mt)

# --------------------------------- PLOTTING --------------------------------- #
using GLMakie
using DynamicalSystems
function plot_RM!(fig, t, tr, ϵ; tsmode="scatterline")
    RM = RecurrenceMatrix(tr, ϵ)
    ax = Axis(fig[1,1], ylabel="j", xlabel="i")
    xs, ys = coordinates(RM)
    scatter!(ax, xs, ys; color = :black, markersize = 3)
    ax.limits = ((0, size(RM, 1)+1), (0, size(RM, 2)+1));
    ax.aspect = 1
    ax = Axis(fig[1,2], ylabel="x", xlabel="t")
    if tsmode == "scatterline" scatterlines!(ax, t, tr[:,1], color=:black)
    else lines!(ax, t, tr[:,1], color=:black) end
end 

