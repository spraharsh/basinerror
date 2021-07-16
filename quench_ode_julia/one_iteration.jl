using DifferentialEquations
using PyCall
using LSODA
using Sundials
pot =pyimport("pele.potentials")

"""
This function ideally defines an interface to the step

the goal of the function is to take a single step with the ode solver and useful information
that we can use to define a new step
Parameters:
-------------

integrator: integrator to be modified with a new step being taken
primary point of modification

# grad_func: callable
# grad_func(x) = - \\grad{E(x)} where E is the energy. and x is the position

# position: array
# flat array of positions

# optional: Jacobian
# needed for implicit solvers?

dt: defines the step
Returns
-------------
tnew: parametrized time
xnew: array
New positions after taking a single ode step
gradnew: new gradient

we're not using the exclamation point because we want to import this in python
"""


function one_iteration!(integrator)
    step!(integrator)
    (integrator.t, integrator.u, get_du(integrator))
end




"""
Initializes the ODE problem for our situation
"""
function initialize_ode(get_negative_grad, x_0, tmax)
    tspan = (0.0, 1)
    prob = ODEProblem(get_negative_grad, x_0, tspan)
    init(prob, QNDF())
end


global f_eval = 0


function f(u, p, t)
    global f_eval += 1
    - (2*u)*sech(t)^2
end



th = pot.InversePower(2.5, 2, [1.0, 1.0], ndim=2)

energy = th.getEnergy([1.0, 1.0, 2.0, 1.0])


println(energy)


if abspath(PROGRAM_FILE) == @__FILE__
    u0 = [1, 1]
    integrator = initialize_ode(f, u0,10000)
    println(one_iteration!(integrator))
    println(integrator.destats)
end

# tspan = (0, 10)
# x0 = [1, 1]
# f(u, p, t) = -2*u
# prob = ODEProblem(f, x0, tspan)
# init(prob, CVODE_BDF())
# step!()
# sol = solve(prob, lsoda())






