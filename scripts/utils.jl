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
using StatsBase, LinearAlgebra
function histogram(v, numbins)
	bins = range(minimum(v), maximum(v), length=numbins)
		a = fit(Histogram, v, bins)
		normalize(a, mode=:pdf)
		return a.weights, bins
end

timederivative(sol) = [sol(t, Val{1}) for t ∈ sol.t]
norm(v) = sum(v.^2)