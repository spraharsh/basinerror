#!/usr/bin/env python3
import numpy as np
import os
from pele.utils.cell_scale import get_box_length, get_ncellsx_scale
from enum import Enum
import json

BASE_DIRECTORY = '/home/praharsh/Dropbox/research/bv-libraries/basinerror/datainv'


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
    n_part = 16
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


def generate_save_secondary_params(par, folder):
    """ 
    Generates a list of secondary parameters hs_radii,
    boxlength and scale and initial coords and saves them to the folder
    """
    np.random.seed(0)
    # hs_radii = (par.radius_std.value *
    #             np.random.randn(par.n_part.value)
    #             + par.radius_mean.value)

    # for binary mixture
    #
    npart_by_2 = par.n_part.value // 2
    hs_radii = np.array(
        list(par.r1.value + par.rstd1.value * np.random.randn(npart_by_2)) +
        list(par.r2.value +
             par.rstd2.value * np.random.randn(par.n_part.value - npart_by_2)))
    print(hs_radii)
    print(len(hs_radii))
    box_length = get_box_length(hs_radii, par.ndim.value, par.phi.value)
    initial_coords = np.random.rand(
        par.n_part.value * par.ndim.value) * box_length
    print(len(initial_coords), "initial_coords length")
    np.savetxt(folder + '/hs_radii.txt', hs_radii, delimiter=',')
    np.savetxt(folder + '/initial_coords.txt', initial_coords, delimiter=',')
    np.savetxt(folder + '/box_length.txt', [box_length], delimiter=',')


def load_secondary_params(folder):
    """Loads secondary parameters, should write this into a class
    """
    hs_radii = np.loadtxt(folder + '/hs_radii.txt', delimiter=',')
    initial_coords = np.loadtxt(folder + '/initial_coords.txt', delimiter=',')
    box_length = np.loadtxt(folder + '/box_length.txt', delimiter=',')
    return hs_radii, initial_coords, float(box_length)


def generate_params_from_foldername(foldername):
    """ Automatically generates necessary params from the foldername
    """


if True:
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
    # fname = foldername + '/systemdata.json'
    # Inverse power saving
    save_params(param, foldername)
    generate_save_secondary_params(param, foldername)
