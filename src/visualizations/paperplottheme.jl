
inches_to_pts(l) = 72*l
const columnsize_pt = inches_to_pts(3.54) #inches
width_pt = 2*columnsize_pt;
height_pt = 1.2*width_pt;

_FONTSIZE = 20
_LABELSIZE = 24
_SPINEWIDTH = 1.8
papertheme = Theme(
    fontsize = _FONTSIZE,
    Axis=(
        spinewidth=_SPINEWIDTH, rightspinevisible=false, topspinevisible=false, xtickalign=1, ytickalign=1, ygridvisible=false, xgridvisible=false,
        xlabelsize = _LABELSIZE, ylabelsize = _LABELSIZE
        ),
    Axis3=(
        xspinewidth=_SPINEWIDTH, yspinewidth=_SPINEWIDTH, zspinewidth=_SPINEWIDTH,
        xtickalign=1, ytickalign=1,
        xlabelsize = _LABELSIZE, ylabelsize = _LABELSIZE, zlabelsize = _LABELSIZE
        ),
    )
# set_theme!(papertheme)

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

_FONTSIZE = 9
_LABELSIZE = 9
_SPINEWIDTH = 1.0
figure3 = Theme(
    fontsize = _FONTSIZE,
    Axis=(
        spinewidth=_SPINEWIDTH, rightspinevisible=false, topspinevisible=false, xtickalign=1, ytickalign=1, ygridvisible=false, xgridvisible=false,
        xlabelsize = _LABELSIZE, ylabelsize = _LABELSIZE, 
        xlabelpadding=2, ylabelpadding=2,
        ),
    Axis3=(
        xspinewidth=_SPINEWIDTH, yspinewidth=_SPINEWIDTH, zspinewidth=_SPINEWIDTH,
        xtickalign=1, ytickalign=1,
        xlabelsize = _LABELSIZE, ylabelsize = _LABELSIZE, zlabelsize = _LABELSIZE,
        xlabeloffset=2, ylabeloffset=2, zlabeloffset=1, protrusions=0, viewmode=:stretch,
        ),
)
# set_theme!(papertheme)
