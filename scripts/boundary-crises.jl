using GLMakie, DynamicalSystems

# ----------------------------------------- Logistic map has boundary crisis at r = 4 ---------------------------------------- #

# Cobweb interactive
# the second range is a convenience for intermittency example of logistic
using InteractiveDynamics 
n = 2000
Ttr = 2000
rrange = 3.4:0.0005:4.0; #interesting part of the diagram 
rrange = 3.82:0.00001:3.88; #period 3
rrange = 3.855:0.000001:3.859; #to find the boundary crisis point; critical value around 3.8568; can pick before 3.8567 and after at 3.8569
lo = Systems.logistic(0.4; r = rrange[1]);
# interactive_cobweb(lo, rrange, 5)
output = orbitdiagram(lo, 1, 1, rrange; n, Ttr);

L = length(rrange);
x = Vector{Float64}(undef, n*L);
y = copy(x);
for j in 1:L
    x[(1 + (j-1)*n):j*n] .= rrange[j]
    y[(1 + (j-1)*n):j*n] .= output[j]
end

fig, ax = scatter(x, y; axis = (xlabel = L"r", ylabel = L"x"),
    markersize = 0.8, color = ("black", 0.05),
)

T = 2000
fig = Figure(resolution=(1920, 1080)) 
axs = []
for (i, r) ∈ enumerate([4.0, 4.00001])
	lo = Systems.logistic(0.4; r);
	tr = trajectory(lo, T; Ttr)
	x = tr[:,1];  t = Ttr:Ttr+T
	ax = Axis(fig[1,i], ylabel="x_n", xlabel="n", title=["Before interior crisis (r = $(r))", "After interior crisis (r = $(r))"][i]); push!(axs, ax)
	# scatterlines!(ax, t, x, color=(:black, 1.0), markersize=2)
	scatter!(ax, t, x, color=(:black, 1.0), markersize=2)
	ylims!(-6, 1)
end
save("boundarycrisis-logisticmap.png", fig, px_per_unit=4)


# ---------------------------------------------------- Ikeda map also famous for boundary crisis --------------------------------------------------- #
# translation from ds to my own notation: t = χ, c=C, d=K (= a in Alligood), a=A(=r in Alligood), b=B 
b = 0.9 
c = 0.4 
d = 6.0
a = 1.0 # a=A=r varying
u0 = [0.875, -1.1]
ik = Systems.ikedamap(u0;a, b, c, d)

T = 1e6
Ttr = 0
# tr = trajectory(ik, T; Ttr)
# x = tr[:,1]; y = tr[:,2]

# fig = Figure() 
# ax = Axis(fig[1,1])
# scatter!(ax, x, y, markersize=2, color=:black)


using DataStructures
function histmeasure(tr, numdigits)
	v = [[tr[i,1], tr[i,2]] for i=1:length(tr)]
	vround = [round.(el, digits=numdigits) for el in v]
	c = counter(vround)
	measure = [c[elround] for elround in vround ]
	measure ./ length(v)
end


fig = Figure(resolution=(1920, 1080)) 
axs = []
for (i, a) ∈ enumerate([0.997, 1.003])
	ik = Systems.ikedamap(;a, b, c, d)
	tr = trajectory(ik, T, u0; Ttr)
	x = tr[:,1]; y = tr[:,2]; t = Ttr:Ttr+T
	measure = histmeasure(tr, 2)
	ax = Axis(fig[1:3,i], ylabel="y", xlabel="x", title=["Before boundary crisis (a = $(a))", "After boundary crisis (a = $(a))"][i]); push!(axs, ax)
	scatter!(ax, x, y, markersize=3, color=measure)
	xlims!(-0.6, 6.2)
	ylims!(-2.8, 6)
	ax = Axis(fig[4, i], ylabel="x")
	lines!(t, x, color=:black)
	xlims!(0, 20000)
	ylims!(-0.6, 6.2)
	ax = Axis(fig[5, i], ylabel="y", xlabel="t")
	lines!(t, y, color=:black)
	xlims!(0, 20000)
	ylims!(-0.6, 6.2)
end
save("boundarycrisis-ikedamap.png", fig, px_per_unit=4)
