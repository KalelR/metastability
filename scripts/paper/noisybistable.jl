using DrWatson
@quickactivate "metastability"
using GLMakie
using DifferentialEquations, CurveFit, DelimitedFiles
include("$(srcdir())/paperplottheme.jl")
include("$(scriptsdir())/utils.jl")
include("$(srcdir())/systems/duffing.jl")




# fig, axs = integrate_and_plot_noisybistable(T=1e8; obtaindwelltime="generate"); #for the paper, generating the dwell times
# fig, axs = integrate_and_plot_noisybistable(T=1e4; obtaindwelltime="readfromfile"); #quick version: integrates for a short time and reads dwell times from long integration made previously
# fig
# safesave("$(plotsdir())/mechanisms/paper-version/noisybistable.png", fig)

