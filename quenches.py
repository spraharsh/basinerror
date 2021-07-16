#!/usr/bin/env python3
"""
Simulates a 2/3d system of particles in a periodic chooses a single
 particles $N=0$ and checks the projection of the basin of attraction,
 in that particles coordinates with the steepest descent data. This is
 done for the $LJ$ cut potential.
"""
from numpy.core.arrayprint import ComplexFloatingFormat
from pele.potentials import HS_WCA
# , CVODEBDFOptimizer
from pele.optimize import MixedOptimizer, GradientDescent_CPP, LBFGS_CPP
from pele.utils.cell_scale import get_ncellsx_scale
from pele.utils.cell_scale import get_box_length
from pele.distance import Distance
from pele.optimize import CVODEBDFOptimizer
from PyCG_DESCENT import CGDescent
import numpy as np


np.random.seed(0)  # check how to use seed sequences later on
QUENCH_FOLDER_NAME = 'lbfgs_m1_final'

def quench_mixed_optimizer(x0, pot, nsteps=2000, **kwargs):
    """ 
        \"Subroutine\" for quenching mixed optimizer, add subtract yada yada
        to control how information gets returned, basically simply passing
        pot, x0 with these default parameters should give identical results
        between different pieces.
    """
    # really need to make sure that this is correctly chosen
    mx_opt = MixedOptimizer(pot, x0, **kwargs)
    mx_opt.run(nsteps)
    res = mx_opt.get_result()
    return res


def quench_LBFGS(x0, pot, steps=2000, **kwargs):
    """ \"Subroutine\" for quenching LBFGS, add subtract yada yada
        to control how information gets returned, basically simply passing
        pot, x0 with these default parameters should give identical results
        between different pieces.
    """
    lbfgs = LBFGS_CPP(x0, pot, steps=2000, **kwargs)
    lbfgs.run(nsteps)
    res = LBFGS_CPP.get_result()
    return res


def quench_steepest(x0, pot, tol=1e-4, nsteps=10000, stepsize=1e-4):
    """ "Subroutine" for quenching steepest descent, add subtract yada yada
        to control how information gets returned, basically simply passing
        pot, x0 with these default parameters should give identical results
        between different pieces.
    """
    steepest = GradientDescent_CPP(x0, pot ,tol=tol, stepsize=stepsize)
    steepest.run(nsteps)
    res = steepest.get_result()
    return res


def quench_cvode_opt(x0, pot, nsteps=1000000, **kwargs):
    """ "Subroutine" for quenching steepest descent, add subtract yada yada
        to control how information gets returned, basically simply passing
        pot, x0 with these default parameters should give identical results
        between different pieces.
    """
    print(kwargs)
    cvode = CVODEBDFOptimizer(pot, x0, **kwargs)
    cvode.run(nsteps)
    res = cvode.get_result()
    return res


def quench_pycg_descent(x0, pot, nsteps=2000, **kwargs):
    """ "subroutine" for quenching using cvode.

    Parameters
    ----------
    pot: BasePotential
        potential to quench
    x0: ndarray[ndim*nparticles]
        initial coordinates
    nsteps: int
        maximum number of steps
    """
    cgd = CGDescent(x0, pot, **kwargs)
    res= cgd.run(nsteps)
    return res

if __name__ == "__main__":
    # general parameters
    ndim = 2
    phi = 0.7
    # particle parameters
    nparticles = 8
    radius_mean = 1.0
    radius_std = 0.05

    # potential parameters
    use_cell_lists = False
    pot_sca = 0.1
    radius_sca = 0.9
    eps = 1.0
    rcut = 2.5

    # optimizer parameters (Mixed Optimizer)
    tol = 1e-4
    T = 1
    step = 1
    conv_tol = 1e-2
    conv_factor = 2
    nsteps = 2000
    H0 = 1e-10
    

    # Set up the coordinates and potential

    # choose a particular particle location we want to vary (particle 0,
    # without loss of generality)
