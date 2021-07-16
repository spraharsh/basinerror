"""
Quench that solves the differential equation
\frac{dx}{dt} = - \grad{V(x)} to find the attractor corresponding to the point
"""



include("one_iteration.jl")



# function one_iteration!(integrator)
#     step!(integrator)
#     (integrator.t, integrator.u, get_du(integrator))
# end

# """
# Initializes the ODE problem for our situation
# """
# function initialize_ode(get_negative_grad, x_0, tmax)
#     tspan = (0.0, 1)
#     prob = ODEProblem(get_negative_grad, x_0, tspan)
#     init(prob, QNDF())
# end
