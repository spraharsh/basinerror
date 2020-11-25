from typing import List
import numpy as np
from numpy.lib.arraysetops import unique
from basinerror import quench_steepest
from basinerror import quench_mixed_optimizer
from basinerror import quench_cvode_opt
from params import load_params, load_secondary_params
from params import BASE_DIRECTORY
from pele.potentials import HS_WCA, InversePower
from pele.distance import Distance
from pele.optimize._quench import fire, steepest_descent, modifiedfire_cpp
from scipy.integrate import ode
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee
from hilbertcurve.hilbertcurve import HilbertCurve



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
                             use_cell_lists=True,
                             ndim=sysparams.ndim.value,
                             radii=hs_radii * 1.0,
                             boxvec=boxv)
    print(box_length, "box_length")
    # initial_coords = np.loadtxt(foldpath + '/coords_of_minimum.txt',
    #                             delimiter=',')
    print(initial_coords, 'initial_coords')
    hcurvepoints = np.floor(initial_coords
                            / box_length
                            * 4*8).reshape(int(len(initial_coords)/sysparams.ndim.value),
                                            int(sysparams.ndim.value)).astype(int).tolist()
    # initiate hilbert curve filling
    print(len(hcurvepoints), "hcurvepoints")

    p = 6
    N = sysparams.ndim.value
    print(hcurvepoints)
    hilbert_curve = HilbertCurve(p, N)
    hilb_dist_arr = [hilbert_curve.distance_from_coordinates(coords) for coords in hcurvepoints]
    print(hilb_dist_arr, 'hilbdistarr')
    unique_dist = np.sort(np.unique(hilb_dist_arr))
    print(unique_dist, 'unique hilb dist')
    masks = hilb_dist_arr == unique_dist[3]
    print(np.array(hcurvepoints)[masks], 'starting point')
    print(np.count_nonzero(hilb_dist_arr==unique_dist[0]), 'this')
    local_particle_args = np.argsort(hilb_dist_arr)
    print(sysparams.ndim.value)
    print(np.array(hcurvepoints)[local_particle_args], 'hcurvepoints[localparticleargs]')
    if sysparams.ndim.value == 2:
        hilbertorder = np.array([local_particle_args*2,
                                 local_particle_args*2 +1]).T.flatten()
    else:
        hilbertorder = np.array([local_particle_args*3,
                                 local_particle_args*3 + 1,
                                 local_particle_args*3 + 2]).T.flatten()
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
    s, x, hess = potential.getEnergyGradientHessian(initial_coords)
    import matplotlib.pyplot as plt
    # np.set_printoptions(threshold=4000)
    order = (potential.getAtomOrder(initial_coords))
    order3d = np.array([order*3, order*3 +1, order*3+2]).T.flatten()
    order2d = np.array([order*2, order*2 +1]).T.flatten()
    # hess = hess[hilbertorder][:, hilbertorder]
    # print(order)
    # print(np.max(order))
    # print(order)
    # print(initial_coords)
    bl = (initial_coords[hilbertorder]/
          box_length).reshape(int(len(initial_coords)/
                                  sysparams.ndim.value),
                              int(sysparams.ndim.value))
    # print(bl, 'hilbertorder')
    consecutive_dist_arr = bl[1:]-bl[:-1]
    th =np.array([np.linalg.norm(x) for x in consecutive_dist_arr])
    maxdist = np.max(th)
    # print(th)
    # print(maxdist)
    # print(np.array(hilb_dist_arr)[local_particle_args])
    # hcurvepointsordered = np.floor(initial_coords[hilbertorder]
    #                                / box_length
    #                                * 4).reshape(int(len(initial_coords)/sysparams.ndim.value),
    #                                             int(sysparams.ndim.value)).astype(int).tolist()
    # print(hcurvepointsordered)
    # hess =hess[order3d][:,order3d]
    # plt.imshow(hess[order][:,order], cmap='hot')
    # plt.spy(hess[order][:,order])
    sparse_hess = csr_matrix(hess)
    cmorder = reverse_cuthill_mckee(sparse_hess)
    orderedhess = hess[cmorder][:, cmorder]
    # plt.spy(hess[cmorder][:, cmorder])
    # plt.spy(hess[hilbertorder][:, hilbertorder])
    bandlist = []
    for diagarg, hessrow in enumerate(orderedhess):
        nonzeroarg = np.nonzero(hessrow)
        if len(nonzeroarg[0]) > 0:
            maxband = np.max(nonzeroarg) - diagarg
            bandlist.append(maxband)
    print(np.max(bandlist), 'bandwidth')
    plt.spy(orderedhess)
    print(np.count_nonzero(hess),'count')
    # plt.spy(hess[order3d][:, order3d])
    plt.show()
    # np.savetxt('/home/praharsh/h1.txt', hess, delimiter=',')
    quit()
    ret = quench_cvode_opt(potential, initial_coords, tol=1e-9, rtol=1e-10, atol=1e-10)


    finalcoords = ret.coords
    # print(ret.coords)
    np.savetxt(foldpath + '/coords_of_minimum.txt', finalcoords, delimiter=',')
    # print(ret)
    # print(ret.coords)
    # print(ret2.coords)
    # print(ret2.coords-ret.coords)
    # print(finalcoords)

    E, V, H = potential.getEnergyGradientHessian(finalcoords)
    # print(potential.getEnergyGradient(initial_coords))


if __name__ == "__main__":
    # foldnameHSWCA = ("ndim=2phi=0.7seed=0n_part=8radius_mean=1.0"
    #             + "radius_std=0.05use_cell_lists=0pot_sca=0.1radius_sca=0.9"
    #             + "eps=1.0rcut=2.5")
    # FindMinimumHSWCA(foldnameHSWCA)
    foldnameInversePower ='ndim=2phi=0.9seed=0n_part=256r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=1.0power=2.5eps=1.0'
    FindMinimumInversePower(foldnameInversePower)
