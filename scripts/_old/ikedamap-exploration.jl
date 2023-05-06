using DrWatson
@quickactivate "metastability"
using GLMakie, DynamicalSystems

# ----------------------------------------- Logistic map has interior crisis at the end of periodic windows ---------------------------------------- #
dirname = "interiorcrisis"
include("$(scriptsdir())/utils.jl")

a = 0.84
b = 0.9 
c = 0.4 
# d = -7.24
d = 7.1
dvalues = 3:0.01:8
# dvalues = [7.1, 7.3]

T = 1e6
Ttr = Int(1e3)
n = 100


ds = Systems.ikedamap(;a, b, c, d)
output = orbitdiagram(ds, 1, 4, dvalues; Ttr, n, u0=rand(2))

# x = dvalues
y = [Float64[] for i=1:length(output)]
x = [Float64[] for i=1:length(output)]
for i=1:length(output)
    y[i] = unique(output[i])
    x[i] = [dvalues[i] for j=1:length(y[i])]
end
x = reduce(vcat, x)
y = reduce(vcat, y)

fig = Figure()
ax = Axis(fig[1,1])
scatter!(ax, x, y, markersize=4)