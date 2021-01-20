import numpy as np
from quenches import quench_steepest
from quenches import quench_mixed_optimizer
from quenches import quench_cvode_opt
from params import load_params, load_secondary_params
from params import BASE_DIRECTORY
from pele.potentials import HS_WCA, InversePower
from pele.distance import Distance
from pele.optimize._quench import fire, steepest_descent, modifiedfire_cpp
from scipy.integrate import ode


def FindMinimumHSWCA(foldname):
    """ Finds the true minimum corresponding:
        Note Gradient Descent step should be adequately small for this to work
    """
    foldpath = BASE_DIRECTORY + '/' + foldname
    sysparams = load_params(foldpath)
    (hs_radii, initial_coords, box_length) = load_secondary_params(foldpath)
    box_length = float(box_length)
    boxv = [box_length] * sysparams.ndim.value
    potential = HS_WCA(use_cell_lists=False,
                       eps=sysparams.eps.value,
                       sca=sysparams.pot_sca.value,
                       radii=hs_radii * sysparams.radius_sca.value,
                       boxvec=boxv,
                       ndim=sysparams.ndim.value,
                       distance_method=Distance.PERIODIC)
    # ret = steepest_descent(initial_coords, potential)
    # ret = fire(initial_coords, potential, iprint=1)
    # E, V = print(potential.getEnergyGradient(initial_coords))
    # ret = quench_mixed_optimizer(potential, initial_coords, conv_tol=1e-100)
    ret = quench_mixed_optimizer(potential,
                                 initial_coords,
                                 conv_tol=0,
                                 nsteps=1000)
    # ret = quench_steepest(potential, initial_coords, stepsize=0.05, nsteps=2000)
    print(ret.nsteps)
    print(ret.nfev)
    np.savetxt(foldpath + '/trueminimum.txt', ret.coords, delimiter=',')


def FindMinimumInversePower(foldname):
    """ Finds the true minimum corresponding:
        Note Gradient Descent step should be adequately small for this to work
    """
    foldpath = BASE_DIRECTORY + '/' + foldname
    sysparams = load_params(foldpath)
    (hs_radii, initial_coords, box_length) = load_secondary_params(foldpath)
    box_length = float(box_length)
    boxv = [box_length] * sysparams.ndim.value
    potential = InversePower(sysparams.power.value,
                             sysparams.eps.value,
                             use_cell_lists=False,
                             ndim=sysparams.ndim.value,
                             radii=hs_radii * 1.0,
                             boxvec=boxv)
    print(box_length, "box_length")
    print(len(initial_coords))
    # initial_coords = np.loadtxt(foldpath + '/coords_of_minimum.txt',
    #                             delimiter=',')
    print(initial_coords, 'initial_coords')
    # ret = steepest_descent(initial_coords, potential, dx=1e-4)
    # ret = modifiedfire_cpp(initial_coords, potential, iprint=1, tol=1e-4)
    # E, V = print(potential.getEnergyGradient(initial_coords))
    # ret = quench_mixed_optimizer(potential, initial_coords, conv_tol=1e-100)
    # ret = quench_mixed_optimizer(potential, initial_coords, conv_tol=1e-100, nsteps=3000)
    # ret = quench_steepest(
    #     potential,
    #     initial_coords,  # make sure right coords are being passed
    #     stepsize=0.0001,
    #     nsteps=1000,
    #     tol=1e-7)
    # ret = quench_cvode_opt(initial_coords, potential, )
    ret = quench_cvode_opt(initial_coords, potential, tol=1e-9, rtol=1e-10, atol=1e-10)

    finalcoords = ret.coords
    print(ret.coords)
    # np.savetxt(foldpath + '/coords_of_minimum.txt', finalcoords, delimiter=',')
    print(ret)
    # print(ret.coords)
    # print(ret2.coords)
    # print(ret2.coords-ret.coords)
    print(finalcoords)
    E, V, H = potential.getEnergyGradientHessian(finalcoords)
    print(np.linalg.eigvals(H), 'all eigenvalues')
    # print(potential.getEnergyGradient(initial_coords))


if __name__ == "__main__":
    # foldnameHSWCA = ("ndim=2phi=0.7seed=0n_part=8radius_mean=1.0"
    #             + "radius_std=0.05use_cell_lists=0pot_sca=0.1radius_sca=0.9"
    #             + "eps=1.0rcut=2.5")
    # FindMinimumHSWCA(foldnameHSWCA)
    foldnameInversePower = "ndim=2phi=0.9seed=0n_part=32r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0"
    FindMinimumInversePower(foldnameInversePower)
