# ---------------------------------------------------------------------------- #
#                                Chaotic Saddle                                #
# ---------------------------------------------------------------------------- #
iswithinneighborhood(point, fp, threshold) = evaluate(Euclidean(), point, fp) ≤ threshold ? true : false
# iswithinneighborhood([3.0, 3.0], fp, 0.1)

using DynamicalSystemsBase:obtain_access, get_state
# function trajectory_discrete(integ, t, u0 = nothing;
#         Δt::Int = 1, Ttr = 0, a = nothing, diffeq = nothing, fp, threshold
#     )
#     !isnothing(u0) && reinit!(integ, u0)
#     Δt = round(Int, Δt)
#     T = eltype(get_state(integ))
#     t0 = current_time(integ)
#     tvec = (t0+Ttr):Δt:(t0+t+Ttr)
#     L = length(tvec)
#     T = eltype(get_state(integ))
#     X = isnothing(a) ? dimension(integ) : length(a)
#     data = Vector{SVector{X, T}}(undef, L)
#     Ttr ≠ 0 && step!(integ, Ttr)
#     data[1] = obtain_access(get_state(integ), a)
# 	t_convergence = 0
#     for i in 2:L
#         step!(integ, Δt)
# 		if iswithinneighborhood(integ.u, fp, threshold) for j=i:L data[j] = SVector{X, T}(fp) end; t_convergence=i; break end
#         data[i] = SVector{X, T}(obtain_access(get_state(integ), a))
#     end
#     return Dataset(data), t_convergence
# end

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

function distribution_times_chaotic_saddle(τs, numbins; mode=:pdf)
	# filter!(x->x>=10, τs)
	bins = collect(range(minimum(τs), maximum(τs), length=numbins+1))
	weights, bins = histogram(τs, bins; mode=mode);
	return weights, bins
end

# function distribution_times_chaotic_saddle(numbins; mode=:pdf, T=1e6, Ttr=1e2, threshold=1, icsperdim=20)
# 	τs = dwell_times_chaotic_saddle(icsperdim; T, Ttr, threshold)
# 	filter!(x->x>=10, τs)
# 	filter!(x->x<T, τs)
# 	bins = collect(range(minimum(τs), maximum(τs), length=numbins+1))
# 	weights, bins = histogram(τs, bins; mode=mode);
# 	return weights, bins
# end

function dwell_times_chaotic_saddle(ik, u0, fp; icsperdim=20, T=1e6, Ttr=1e2, threshold=1)
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

