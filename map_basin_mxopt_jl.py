import julia

from julia import Main
Main.include("quench_ode_julia/mxopt_int.jl")
import numpy as np


from pele.potentials import InversePower
from utils.cell_scale import get_box_length, get_ncellsx_scale
from params import BASE_DIRECTORY, load_params, load_secondary_params
from random_vectors import VEC_8_0, VEC_8_1, VEC_8_2,  VEC_16_0, VEC_16_1, VEC_16_2, VEC_32_0, VEC_32_1, VEC_32_2
from utils.cell_scale import get_box_length, get_ncellsx_scale


def map_binary_inversepower_mxopt_jl(foldername,
                                     particle_coords,
                                     optimizer, opt_param_dict, random_coord_0=0, random_coord_1=-1, z=0):
    """
    Finds whether a point defined by particle_coord
    on the meshgrid correspond to a minimum or not for a 2d
    case.
    """
    foldpath = BASE_DIRECTORY + '/' + foldername
    # import params
    sysparams = load_params(foldpath)
    (hs_radii, initial_coords, box_length) = load_secondary_params(foldpath)
    assert (sysparams.ndim.value == 2)
    minimum_coords = np.loadtxt(foldpath + '/coords_of_minimum.txt',
                                delimiter=',')
    minimum_coords = np.loadtxt(foldpath + '/coords_of_minimum.txt',
                                delimiter=',')
    quench_coords = initial_coords.copy()

    quench_coords = quench_coords + \
    particle_coords[0]*VEC_16_0 + particle_coords[1]*VEC_16_1 + z*VEC_16_2

    # quench_coords = quench_coords + \
    #     particle_coords[0]*VEC_8_0 + particle_coords[1]*VEC_8_1 + z*VEC_8_2
    # print(quench_coords, 'quench coords')
    # box length
    box_length = float(box_length)
    boxv = [box_length] * sysparams.ndim.value
    ncellx_scale = get_ncellsx_scale(hs_radii, boxv)
    print(hs_radii, 'hs_radii')
    print(quench_coords, 'quench_coords')
    print(boxv)
    potential = Main.pot.InversePower(sysparams.power.value,
                             sysparams.eps.value,
                             use_cell_lists=False,
                             ndim=sysparams.ndim.value,
                             radii=hs_radii * 1.0,
                             boxvec=boxv)

    ppot = Main.PythonPotential(potential)
    nls = Main.NewtonLinesearch(ppot, quench_coords, opt_param_dict['tol'])
    mxd = Main.Mixed_Descent(ppot, Main.CVODE_BDF(), nls, quench_coords, opt_param_dict['T'], opt_param_dict['rtol'], opt_param_dict['conv_tol'], opt_param_dict['tol'])
    try:
        Main.run_b(mxd, 10000)
    except:
        print(quench_coords, 'failed here')
        print(initial_coords, 'coords')
        print(len(quench_coords))
    # ret = lbfgs_cpp(quench_coords, potential, tol=1e-8, M=1)

    # This exists because some runs don't have hessian evaluations
    coordarg = 0
    print(mxd.converged, 'converged')
    results = (mxd.optimizer.x0, mxd.converged, 0, mxd.n_g_evals, mxd.iter_number, mxd.n_h_evals)
    # the reason the potential is being passed is because quench coords needs the potential to figure out what to do
    return results