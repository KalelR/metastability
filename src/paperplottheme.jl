
inches_to_pts(l) = 72*l
const columnsize_pt = inches_to_pts(3.54) #inches
width_pt = 2*columnsize_pt;
height_pt = 1.2*width_pt;

papertheme = Theme(
    fontsize = 12,
    Axis=(spinewidth=1.4, rightspinevisible=false, topspinevisible=false, xtickalign=1, ytickalign=1, ygridvisible=false, xgridvisible=false),
    )
set_theme!(papertheme)