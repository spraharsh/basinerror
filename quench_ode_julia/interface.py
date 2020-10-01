#!/home/praharsh/anaconda3/envs/changebranch3 python
"""
This file manages quenches from the DifferentialEquations package in julia.
"""
from pele.optimize import Result
import numpy as np
import julia
from diffeqpy import de
from scipy.integrate import LSODA
from pele.potentials import Harmonic
j = julia.Julia()


# def get_negative_grad(x, p, t):
#     return - 2*x

class potclass:
    def __init__(self):
        self.f_eval = 0
    def get_negative_grad(self, x, p, t):
        self.f_eval += 1
        return -2*x
def get_negative_grad(x, p, t):
    return -2*x


        



def ode_julia_naive(coords, pot, tol=1e-4, nsteps=20000,
                    convergence_check=None, solver_type=de.AutoTsit5(de.Rosenbrock23(autodiff=False)), **kwargs):
    class feval_pot:
        """ wrapper class that interfaces base potential to functions for ode solver
        """
        def __init__(self):
            self.nfev = 0
            self.nhev = 0
        def get_negative_grad(self, x, p, t):
            """
            negative grad is f(u, p, t)
            """
            self.nfev +=1
            return -pot.getEnergyGradient(x.copy())[1]
        def get_energy_gradient(self, x):
            self.nfev +=1
            return pot.getEnergyGradient(x.copy())
        # def get_jacobian(jac, x, p, t):
        #     self.nhev += 1
        #     return -pot.getEnergyGradientHessian(x.copy())[2]
            
            
    
    function_evaluate_pot = feval_pot()
    converged = False
    n = 0
    if convergence_check == None:
        convergence_check = lambda g: np.linalg.norm(g)<tol
    # odefunc = de.ODEFunction(function_evaluate_pot.get_negative_grad)
    # initialize ode problem
    tspan = (0, float('inf'))
    prob = de.ODEProblem(function_evaluate_pot.get_negative_grad, coords, tspan)
    integrator = de.init(prob, solver_type, rtol = 1e-2, atol = 1e-2)
    x_ = np.full(len(coords), np.nan)
    while not converged and n<nsteps:
        xold = x_
        de.step_b(integrator)
        x_ = integrator.u
        n+=1
        converged = convergence_check(de.get_du(integrator))
    res = Result()
    res.coords = x_
    res.energy = pot.getEnergy(x_)
    res.rms = 0
    res.grad = 0
    res.nfev = function_evaluate_pot.nfev
    res.nsteps = n
    res.success = converged
    return res
# u0 = 0.5
# pot = potclass()
# tspan = (0, 1.0)
# print(pot.get_negative_grad(1, 1, 1))
# prob = de.ODEProblem(pot.get_negative_grad,u0, tspan)
# de.init(prob, de.Tsit5())


if __name__=="__main__":
    pot = Harmonic(np.zeros(3), 2)
    res = ode_julia_naive([1.0, 1.0, 1.0], pot)
    print(res.nfev)
    print(res.coords)
    print(res.nfev)





# def initialize_ode(get_negative_grad, x0):
#     tspan = (0, float('inf'))
#     prob = de.ODEProblem(get_negative_grad,x0, tspan)
#     return de.init(prob, de.Tsit5())


# def one_iteration(integrator):
#     de.step_b(integrator)
#     return (integrator.t, integrator.u, de.get_du(integrator))


# u0 = np.array([1, 1])
# pot = potclass()
# integrator = initialize_ode(pot.get_negative_grad, u0)
# ans = one_iteration(integrator)
# print(ans)
# print(pot.f_eval)
# de.step_b(integrator)            # error9












