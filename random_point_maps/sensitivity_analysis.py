from utils.cell_scale import get_box_length, get_ncellsx_scale
from random_point_maps.generate_points_random import get_hs_radii
import numpy as np



def quench_single_inverse_power(particle_coords, foldpath, sub_fold_name
                                optimizer, opt_param_dict):
    """
    figures out the minimum correspoding to a set of particle coords
    Parameters
    ----------
    particle_coords: coords
        initial coords of the particle to run
    foldername: str
        folder definining the run
    sub_fold_name:
        name of subfolder where the run data is stored
    optimizer: optimizer
        quench
    opt_param_dict: dict
        dictionary of parameters for the optimizer
    """
    
    foldpath = BASE_DIRECTORY + '/' + foldername
    sysparams = load_params(foldpath)
    
    radii = get_hs_radii(foldpath, sub_fold_name)
    box_length = get_box_length(radii, sysparams.ndim.value, sysparams.phi.value)

    boxv = [box_length]*sysparams.ndim.value
    ncellx_scale = get_ncellsx_scale(radii, boxv)

    potential = InversePower(sysparams.power.value,
                             sysparams.eps.value,
                             use_cell_lists=True,
                             ndim=sysparams.ndim.value,
                             radii=hs_radii * 1.0,
                             boxvec=boxv)
    
    try:
        ret = optimizer(quench_coords, potential, **parameter_dict)
    except:
        print(quench_coords, 'failed here')
        print(initial_coords, 'coords')
        print(len(quench_coords))
        raise Exception('failure')

    # This exists because some runs don't have hessian evaluations
    try:
        ret['nhev']
    except:
        ret['nhev']=0

    results = (ret.coords, ret.success, ret.nfev, ret.nsteps, ret.nhev)
    return results



def quench_multiple(particle_coords, foldpath, sub_fold_name
                    optimizer, opt_param_dict):
    """
    does multiple quenches and returns results for each quench

    Returns:
        

    """
    


