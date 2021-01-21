import numpy as np






def quench_single_inverse_power(foldpath, sub_fold_name
                                particle_coords, optimizer, opt_param_dict):
    """
    Parameters
    ----------
    foldername: str
        folder definining the run
    sub_fold_name:
        name of subfolder where the run data is stored
    particle_coords: coords
        initial coords of the particle to run
    optimizer: optimizer
        quench
    opt_param_dict: dict
        dictionary of parameters for the optimizer
    """
    foldpath = BASE_DIRECTORY + '/' + foldername
    
    sysparams = load_params(foldpath)
    


