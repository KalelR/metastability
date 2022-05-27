using InteractiveDynamics, GLMakie, DynamicalSystems

# the second range is a convenience for intermittency example of logistic
rrange = 1:0.001:4.0
# rrange = (rc = 1 + sqrt(8); [rc, rc - 1e-5, rc - 1e-3])
# rrange = 3.8248:0.00001:3.8249
# rrange = 3.7:0.0001:3.845
lo = Systems.logistic(0.4; r = rrange[1])
interactive_cobweb(lo, rrange, 5)