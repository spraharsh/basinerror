#!/usr/bin/env python3
from pickle import TRUE
import numpy as np
import os
from pele.utils.cell_scale import get_box_length, get_ncellsx_scale
from enum import Enum
import json
from global_folder_vars import BASE_DIRECTORY
import re





class SystemParamHSWCA(Enum):
    """ Parameters that uniquely define a system along with
        a start coordinate set for minimization. As done for HS_WCA
    """
    # general parameters
    ndim = 3
    phi = 0.7
    # The seed should ideally determine every randomly generated
    # extra parameters
    seed = 0
    # particle parameters
    n_part = 8
    radius_mean = 1.0
    radius_std = 0.05

    # potential parameters
    use_cell_lists = False
    pot_sca = 0.1
    radius_sca = 0.9
    eps = 1.0
    rcut = 2.5


class SystemParamInversePower(Enum):
    """ Parameters that uniquely define a system along a start coordinate for minimization.
        As done for inverse powe
    """
    # general parameters
    ndim = 3
    phi = 0.9
    # The seed should ideally determine every randomly generated
    # extra parameters
    seed = 0

    # particle parameters
    n_part = 8
    radius_mean = 1.0
    radius_std = 0.05

    # potential parameters
    use_cell_lists = False
    power = 2.5  # hertz power
    eps = 1


class SystemParamInversePowerBinary(Enum):
    """ Parameters that uniquely define a system along a start coordinate for minimization.
        As done for inverse powe
    """
    # general parameters
    ndim = 2
    phi = 0.9
    # The seed should ideally determine every randomly generated
    # extra parameters
    seed = 0

    # particle parameters
    n_part = 1024
    r1 = 1.0
    r2 = 1.4
    rstd1 = 0.05
    rstd2 = 0.05 * 1.4

    # potential parameters
    use_cell_lists = False
    power = 2.5  # hertz power
    eps = 1


class SystemParamBLJ(Enum):
    """
    Parameters that uniquely define a system along a start coordinate for minimization.
    As done for Lennard jones binary
    """
    # general parameters
    ndim = 2
    phi = 0.9
    # The seed should ideally determine every randomly generated
    # extra parameters
    seed = 0


def save_params(params_to_save, folder):
    """Saves the parameters
    """
    filename = folder + '/systemdata.json'
    paramdict = {}
    for name, member in params_to_save.__members__.items():
        paramdict[name] = member.value

    with open(filename, 'w') as f:
        json.dump(paramdict, f)


def load_params(folder):
    """Loads the parameters as an enum object
    """
    filename = folder + '/systemdata.json'
    with open(filename, 'r') as f:
        syspar = json.load(f)
    return Enum('SystemParam', syspar)


def generate_save_run_params_ip_binary(par, folder, seed=0):
    """ Generates and saves the run parameters.
        Creates a new folder called sec_params if it does not exist and a subfolder defined by the seed

    Parameters
    ----------
    par : SystemParamInversePowerBinary
        System defining the parameters
    folder : str
        Folder defining our system
    seed : int, optional
        seed for generating random parameters
    """
    np.random.seed(seed)

    # for binary mixture
    npart_by_2 = par.n_part.value // 2
    hs_radii = np.array(
        list(par.r1.value + par.rstd1.value * np.random.randn(npart_by_2)) +
        list(par.r2.value +
             par.rstd2.value * np.random.randn(par.n_part.value - npart_by_2)))
    print(hs_radii, 'hs_radii')
    print(len(hs_radii), 'nparts')
    box_length = get_box_length(hs_radii, par.ndim.value, par.phi.value)
    print(box_length, 'box_length')
    # initial coords are present, to keep continuity with previous code. but they are not necessary
    initial_coords = np.random.rand(
        par.n_part.value * par.ndim.value) * box_length
    print(initial_coords, 'initial coords')
    print(get_ncellsx_scale(hs_radii, [box_length, box_length, box_length]), 'cell scale')
    print(len(initial_coords), "initial_coords length")
    os.makedirs(folder + '/' + 'sec_params' + '/' + str(seed), exist_ok=True)
    np.savetxt(folder + '/' + 'sec_params' + '/' + str(seed) + '/hs_radii.txt', hs_radii, delimiter=',')
    np.savetxt(folder + '/' + 'sec_params' + '/' + str(seed) + '/initial_coords.txt', initial_coords, delimiter=',')
    np.savetxt(folder + '/' + 'sec_params' + '/' + str(seed) + '/box_length.txt', [box_length], delimiter=',')

def load_secondary_params_ip_binary(folder, seed=0):
    """Loads secondary parameters, should write this into a class
    """


    hs_radii = np.loadtxt(folder + '/' + 'sec_params' + '/' + str(seed)  + '/hs_radii.txt', delimiter=',')
    initial_coords = np.loadtxt(folder+ '/' + 'sec_params' + '/' + str(seed)  + '/initial_coords.txt', delimiter=',')
    box_length = np.loadtxt(folder + '/' + 'sec_params' + '/' + str(seed) + '/box_length.txt', delimiter=',')
    return hs_radii, initial_coords, float(box_length)


def make_system_folder(param, datadir):
    """ Makes a folder for the system

    Parameters
    ----------
    param : Enum
        SystemParam defining the system
    datadir : str
        Directory where the system folder is to be saved

    Returns
    -------
    str
        Folder name
    """
    namestr = ''

    for name, member in param.__members__.items():
        namestr += name
        namestr += "="
        namestr += str(member.value)
    foldername = datadir + '/' + namestr
    os.makedirs(foldername, exist_ok=True)
    return foldername

def get_system_parameters_from_folder(folder):
    """ Loads system parameters from the folder string

    Parameters
    ----------
    folder : str
        Folder name
    Returns
    -------
    SystemParam
        SystemParam defining the system
    """
    foldname = "ndim=2phi=0.9seed=0n_part=16r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0"
    # split by floatin point number
    split_foldname = re.split(r"=([-+]?\d*\.\d+|\d+)", foldname)
    split_foldname = split_foldname[:-1]
    param_dict = {}
    for i in range(len(split_foldname)):
         if i % 2 == 0:
             param_dict[split_foldname[i]] = float(split_foldname[i + 1])
    print(param_dict)
    return Enum('SystemParam', param_dict)





if __name__=="__main__":
    # save system parameters in a unique directory
    datadir = BASE_DIRECTORY
    os.makedirs(datadir, exist_ok=True)
    param = SystemParamInversePowerBinary
    namestr = ''
    paramdict = {}
    for name, member in param.__members__.items():
        namestr += name
        namestr += "="
        namestr += str(member.value)
        paramdict[name] = member.value
    print("unique directory name " + namestr)
    # make the directory
    foldername = datadir + '/' + namestr
    os.makedirs(foldername, exist_ok=True)
    
    # Inverse power saving
    save_params(param, foldername)
    generate_save_run_params_ip_binary(param, foldername)
    # load the parameters
    th = get_system_parameters_from_folder(foldername)
    print(th)
    print(th.ndim.value, th.phi.value, th.seed.value, th.n_part.value, th.r1.value, th.r2.value, th.rstd1.value, th.rstd2.value, th.use_cell_lists.value, th.power.value, th.eps.value)
