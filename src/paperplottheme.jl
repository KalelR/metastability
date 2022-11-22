
inches_to_pts(l) = 72*l
const columnsize_pt = inches_to_pts(3.54) #inches
width_pt = 2*columnsize_pt;
height_pt = 1.2*width_pt;

_fontsize = 16
_labelsize = 20
_spinewidth = 1.8
papertheme = Theme(
    fontsize = _fontsize,
    Axis=(
        spinewidth=_spinewidth, rightspinevisible=false, topspinevisible=false, xtickalign=1, ytickalign=1, ygridvisible=false, xgridvisible=false,
        xlabelsize = _labelsize, ylabelsize = _labelsize
        ),
    Axis3=(
        xspinewidth=_spinewidth, yspinewidth=_spinewidth, zspinewidth=_spinewidth,
        xtickalign=1, ytickalign=1,
        xlabelsize = _labelsize, ylabelsize = _labelsize, zlabelsize = _labelsize
        ),
    )
set_theme!(papertheme)

#Thjanks to George Datseries for the code below!
if isdefined(Main, :DrWatson)
    # Extension of DrWatson's save functionality for default CairoMakie saving
    function DrWatson._wsave(filename, fig::Makie.Figure, args...; kwargs...)
        if filename[end-3] != '.'; filename *= ".png"; end
        if isdefined(Main, :CairoMakie)
            CairoMakie.activate!()
            CairoMakie.save(filename, fig, args...; px_per_unit = 4, kwargs...)
        elseif isdefined(Main, :GLMakie)
            GLMakie.activate!()
            GLMakie.save(filename, fig, args...; px_per_unit = 4, kwargs...)
        else
            Makie.save(filename, fig, args...; kwargs...)
        end
    end
end