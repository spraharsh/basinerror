import julia

from julia import Main
Main.include("../../basins.jl/src/minimumassign/minimumassign.jl")
import numpy as np


from utils.cell_scale import get_box_length, get_ncellsx_scale
from generate_points_random import SUB_FOLD_NAME, get_hs_radii
from params import BASE_DIRECTORY, load_params


def quench_single_CVODE_inverse_power_julia(coord_file_name, foldpath, sub_fold_name, optimizer, opt_param_dict):
    """
    quenches a single system through mxopt
    Parameters
    ----------
    coord_file_name: string
        name of the path to the coordinates
    foldername: str
        folder definining the run
    sub_fold_name:
        name of subfolder where the run data is stored
    optimizer: optimizer
        quench
    opt_param_dict: dict
        dictionary of parameters for the optimizer
    """
    
    sysparams = load_params(foldpath)

    # path to quench coords
    quench_coords_path = foldpath + '/' + sub_fold_name + '/' + 'ensemble/' + coord_file_name
    quench_coords = np.loadtxt(quench_coords_path)
    radii = get_hs_radii(foldpath, sub_fold_name)

    box_length = get_box_length(radii, sysparams.ndim.value, sysparams.phi.value)
    box_length
    
    boxv = np.array([box_length]*sysparams.ndim.value)
    
    ncellx_scale = get_ncellsx_scale(radii, boxv)

    # potential = InversePower(sysparams.power.value,
    #                          sysparams.eps.value,
    #                          use_cell_lists=False,
    #                          ndim=sysparams.ndim.value,
    #                          radii=radii * 1.0,
    #                          boxvec=boxv)
    
    print(radii)
    pot = Main.pot.InversePower(
        sysparams.power.value,
        sysparams.eps.value,
        radii,
        ndim = sysparams.ndim.value,
        boxvec = boxv,
        use_cell_lists = True,
        ncellx_scale = ncellx_scale,
    )
    
    ppot = Main.PythonPotential(pot)

    func = Main.gradient_problem_function_pele_b(ppot)
    ba = Main.BasinAssigner(opt_param_dict['rtol'], opt_param_dict['atol'], opt_param_dict['tol'])
    res = Main.find_corresponding_minimum(ba, func, quench_coords, int(10**8), ppot)

    try:
        None
    except:
        print("exception occured")
        # if exception occurs, treat as failure. this is for rattlers
        # not enough failures occur that it would make a difference to not just assume this never happens
        # but we shoudl switch this out
        # but jic
        return (quench_coords, False, 0, 0, 0, 0, 0)

    results = (res.coords, res.success, res.ngeval, res.nsolve, res.nheval, 0, 0)
    
    print(quench_coords_path)
    return results