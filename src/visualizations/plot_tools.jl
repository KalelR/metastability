
const columnsize = 3.54 #inches
inches_to_pts(l) = 72*l
const columnsize_pt = inches_to_pts(3.54) #inches

function plot_raster!(tr, t, th; ax=nothing, fig=nothing, mksize=5, alpha=1.0, plotbursts=true, plotspikes=false)
	spiketimes = spiketimes_sos(tr, t, th)
	bursttimes = [   spiketimes[i][findall(x->x≥burst_threshold, diff(spiketimes[i])) .+ 1] for i=1:length(spiketimes)  ]
	if isnothing(ax) fig = Figure();	ax = Axis(fig[1, 1], xgridvisible=false, ygridvisible=false)  end
	for i=1:N
		plotspikes && scatter!(ax, spiketimes[i], [i for ii=1:length(spiketimes[i])], color=:grey, markersize=mksize, alpha=1.0)
		plotbursts && scatter!(ax, bursttimes[i], [i for ii=1:length(bursttimes[i])], color=:black, markersize=mksize, alpha=1.0)
	end
	ax.ylabel = "unit idx"
	# ax.xlabel = "t"
	return fig, ax
end

function makiesave(filename, fig; kwargs...)
    mkpath(dirname(filename))
    save(filename, fig; kwargs...)
end


function plane!(ax, coords; alpha=1.0)
    ps = [Point3(coord) for coord in coords]
    mesh!(ax, [ps[1], ps[2], ps[3]], color=(:red, alpha), shading=false, transparency=true)
    mesh!(ax, [ps[4], ps[2], ps[3]], color=(:red, alpha), shading=false, transparency=true)
end

Text = Union{Symbol, <: AbstractString}

"""
    subplotgrid(m, n; sharex = false, sharey = false, kwargs...) -> fig, axs
Create a grid of `m` rows and `n` columns axes in a new figure and return the figure and the
matrix of axis. Optionally make every row share the y axis, and/or every column
to share the x axis. In this case, tick labels are hidden from the shared axes.
"""
function subplotgrid(m, n;
        sharex = false, sharey = false, coltitles = nothing, rowtitles = nothing,
        xlabels = nothing, ylabels = nothing, title = nothing, kwargs...
    )
    fig = Makie.Figure(;kwargs...)
    axs = Matrix{Axis}(undef, m, n)
    for i in 1:m
        for j in 1:n
            axs[i,j] = Axis(fig[i,j])
        end
    end
    if sharex
        for j in 1:n
            Makie.linkxaxes!(axs[:,j]...)
            for i in 1:m-1; Makie.hidexdecorations!(axs[i,j]; grid = false); end
        end
    end
    if sharey
        for i in 1:m # iterate through rows
            Makie.linkyaxes!(axs[i,:]...)
            for j in 2:n; Makie.hideydecorations!(axs[i,j]; grid = false); end
        end
    end
    if !isnothing(coltitles)
        for j in 1:n
            axs[1, j].title = coltitles[j]
        end
    end
    if !isnothing(rowtitles)
        for j in 1:m
            axs[j, 1].title = rowtitles[j]
        end
    end
    if !isnothing(xlabels)
        for j in 1:n
            axs[end, j].xlabel = xlabels isa Text ? xlabels : xlabels[j]
        end
    end
    if !isnothing(ylabels)
        for i in 1:m
            axs[i, 1].ylabel = ylabels isa Text ? ylabels : ylabels[i]
        end
    end
    if !isnothing(title)
        Label(fig[0, :], title, valign = :bottom,
            padding = (0, 0, 5, 0), tellheight = true, tellwidth = false,
            font = "TeX Gyre Heros Bold", # same font as Axis titles
        )
    end
    return fig, axs
end


function nearest_powers(v)
    vcat(floor.(Int64, log10.( max.(1e-8, minimum.(v[1:end-1])) )), ceil(Int64, log10(maximum(v))))
end

function superscriptnumber(i::Int)
    if i < 0
        c = [Char(0x207B)]
    else
        c = []
    end
    for d in reverse(digits(abs(i)))
        if d == 0 push!(c, Char(0x2070)) end
        if d == 1 push!(c, Char(0x00B9)) end
        if d == 2 push!(c, Char(0x00B2)) end
        if d == 3 push!(c, Char(0x00B3)) end
        if d > 3 push!(c, Char(0x2070+d)) end
    end
    return join(c)
end

function get_ticks_in_powers(yvalues)
    powers = nearest_powers(yvalues)
    yticks = sort(10.0 .^ powers); ytickstring = ["10$(superscriptnumber(power))" for power in powers]
    return (yticks, ytickstring)
end

function fix_y_ticks_limits!(ax, yvals)
    yticks = get_ticks_in_powers(yvals); 
    ax.yticks = yticks; 
    ylims!(ax, yticks[1][1], yticks[1][end])
end

function fix_x_ticks_limits!(ax, vals; powers=true)
    xticks = powers ? get_ticks_in_powers(vals) : (vals, string.(vals))
    ax.xticks = xticks; 
    xlims!(ax, xticks[1][1], xticks[1][end])
end

function plot_arrow_following_data!(fig, ax, pos_rough, data)
    mindist, idxmin = findmin(mapslices(x->evaluate(Euclidean(), pos_rough, x), data, dims=2))
    idx_closest_point = idxmin[1]
    p1 = data[idx_closest_point, :] #data point closest to my choic e
    p2 = data[idx_closest_point+1, :] #next data point

    ps = [Point3f(p1)] #origin of arrow 
    ns = [Point3f(p2 .- p1)] #component of arrow (diff)
    arrows!(ax, ps, ns, lengthscale=0.3, arrowsize=Vec3f(6e-3, 6e-3, 1e-2), align = :center, fxaa=true) 
    return fig, ax
end