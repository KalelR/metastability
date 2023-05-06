function recurrencetimes(M; Δt=1)
    RTs = [typeof(Δt)[] for i=1:size(M,1)]
    for i=1:size(M,1)
        rts = diff(findall(x->x==true, M[i, :]))
        RTs[i] = rts .* Δt
    end
    return RTs
end

using StatsBase
function histrecurrencetimes(M; Δt=1.0)
    rts = recurrencetimes(M; Δt)
    rts = reduce(vcat, rts)
    bins = minimum(rts):Δt:maximum(rts)+Δt
    return bins[1:end-1], fit(Histogram, rts, bins, closed=:left).weights
end