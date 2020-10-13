#!/usr/bin/env python3
"""
Simulates a 2/3d system of particles in a periodic chooses a single
 particles $N=0$ and checks the projection of the basin of attraction,
 in that particles coordinates with the steepest descent data. This is
 done for the $LJ$ cut potential.
"""
from pele.potentials import HS_WCA
from pele.optimize import MixedOptimizer, GradientDescent_CPP, LBFGS_CPP # , CVODEBDFOptimizer
from pele.utils.cell_scale import get_ncellsx_scale
from pele.utils.cell_scale import get_box_length
from pele.distance import Distance
from pele.optimize import CVODEBDFOptimizer
import numpy as np
np.random.seed(0)  # check how to use seed sequences later on

QUENCH_FOLDER_NAME = 'cvdoeopt'


def quench_mixed_optimizer(pot,
                           x0,
                           tol=1e-4,
                           T=1,
                           step=1,
                           conv_tol=1e-2,
                           conv_factor=2,
                           nsteps=2000,
                           H0=1e-10,
                           rtol=1e-4,
                           atol=1e-4):
    """ "Subroutine" for quenching mixed optimizer, add subtract yada yada
        to control how information gets returned, basically simply passing
        pot, x0 with these default parameters should give identical results
        between different pieces.
    """
    mx_opt = MixedOptimizer(pot, x0, tol, T, step, conv_tol=conv_tol, conv_factor=conv_factor,
                            nsteps=nsteps, rtol =rtol, atol=atol)
    mx_opt.run(nsteps)
    res = mx_opt.get_result()
    return res


def quench_LBFGS(pot, x0, tol=1e-4, M=4, nsteps=2000):
    """ "Subroutine" for quenching LBFGS, add subtract yada yada
        to control how information gets returned, basically simply passing
        pot, x0 with these default parameters should give identical results
        between different pieces.
    """
    lbfgs = LBFGS_CPP(pot, x0, tol=tol, M=M)
    lbfgs.run(nsteps)
    res = LBFGS_CPP.get_result()
    return res


def quench_steepest(pot, x0, tol=1e-4, nsteps=10000, stepsize=1e-4):
    """ "Subroutine" for quenching steepest descent, add subtract yada yada
        to control how information gets returned, basically simply passing
        pot, x0 with these default parameters should give identical results
        between different pieces.
    """
    steepest = GradientDescent_CPP(pot, x0, tol=tol, stepsize=stepsize)
    steepest.run(nsteps)
    res = steepest.get_result()
    return res




def quench_cvode_opt(pot, x0, tol=1e-4, nsteps=10000, atol=1e-7, rtol=1e-7):
    """ "Subroutine" for quenching steepest descent, add subtract yada yada
        to control how information gets returned, basically simply passing
        pot, x0 with these default parameters should give identical results
        between different pieces.
    """
    cvode = CVODEBDFOptimizer(pot, x0, tol=tol, atol=atol, rtol=rtol)
    cvode.run(nsteps)
    res = cvode.get_result()
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
