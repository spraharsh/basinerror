""" 
    Generates and saves a completely random set of configurations
    """

import numpy as np

from numpy.core.fromnumeric import size
from utils.cell_scale import get_box_length
from params import BASE_DIRECTORY, load_params
import numpy as np
from numpy.random import Generator, PCG64, SeedSequence
import os


RANDOM_POINTS_FOLDER_NAME = 'random_points'

# SUB_FOLD_NAME = 'random_points_run'


def generate_random_radii_inverse_power(params, seed_radii, foldpath, subfoldname=SUB_FOLD_NAME):
    """ Generates radii for inverse potential and saves accordingly
    Args:
       foldname folder name
       seed seed
    Returns:
        Points
    """
    n_part_by_2 = params.n_part.value // 2

    # default generator passed with the seed
    rng = Generator(PCG64(seed_radii))
    # first generate N samples with mean 0 and std 1
    samples = rng.standard_normal(params.n_part.value)

    samples_1 = samples[:n_part_by_2]
    samples_2 = samples[n_part_by_2:]

    print(samples_2)
    radii_1 = list(params.r1.value + params.rstd1.value*samples_1)
    radii_2 = list(params.r2.value + params.rstd2.value*samples_2)
    radii = np.array(radii_1+radii_2)
    os.makedirs(foldpath + '/' + subfoldname, exist_ok=True)
    np.savetxt(foldpath + '/' + subfoldname + '/radii.txt', radii)
    return radii


def get_hs_radii(foldpath, subfoldname, delimiter=','):
    return np.loadtxt(foldpath + '/' + subfoldname + '/radii.txt')


def generate_random_configuration_single(box_length, n_part, n_dim, seed_coords):
    """ Given a box length, radii and seeds, generates a set of coords

    Args:
       radii radii of particles in the box
       box_length box length
       seed_coords seed for random number generator

    Returns:
        generates radii
    """

    rng = Generator(PCG64(seed_coords))
    initial_coords = rng.uniform(low=0, high=box_length, size=n_part*n_dim)
    return initial_coords


def generate_random_configurations(box_length, n_part, n_dim, n_ensemble, seed_entropy):
    """ Generates a set of configurations in an ensemble

    Parameters
    ----------
    box_length : float
        length of the box
    n_part : int
        number of particles
    n_dim : int
        dimension of the system
    seed_coords : int
        seed for random number generator
    n_ensemble : int
        number of configurations to generate

    Returns
    -------
    list of configurations : list of numpy arrays
        list of configurations in the ensemble
    seeds : list of ints
        list of seeds for random number generator
        for each configuration
    """
    sq = SeedSequence(seed_entropy)
    seeds = sq.spawn(n_ensemble)
    coords_list = []
    for seed in seeds:
        print(seed)
        initial_coords = generate_random_configuration_single(
            box_length, n_part, n_dim, seed)
        coords_list.append(initial_coords)
    return (coords_list, seeds)




def save_random_configurations(coords_list, seeds, folder):
    """
    Saves random configurations labelled by their spawn key
    """
    assert len(coords_list) == len(seeds)
    config_base_path = folder + '/' + RANDOM_POINTS_FOLDER_NAME
    for i, coords in enumerate(coords_list):
        config_filename = config_base_path + '/' + \
            str(seeds[i].spawn_key[0]) + '.txt'
        np.savetxt(config_filename, coords, delimiter=',')


def construct_save_random_configurations(base_folder, n_ensemble, n_part, n_dim, box_length, seed_entropy=None):
    """
    Generates a set of configurations in an ensemble and saves them

    Parameters
    ----------
    base_folder : str
        folder where the ensemble is stored
    n_ensemble : int
        number of configurations to generate
    n_part : int
        number of particles
    n_dim : int
        dimension of the system
    box_length : float
        length of the box
    seed_entropy : int, optional
        entropy of the seed, by default None

    Returns
    -------
    coords_list : list of numpy arrays
        list of coordinates in the ensemble
    seeds : list of ints
        list of seeds for random number generator
    """
    base_folder = str(base_folder)
    os.makedirs(base_folder, exist_ok=True)
    entropy_file = base_folder + '/' + 'entropy.txt'
    if seed_entropy is None:
        # check if entropy file exists
        if os.path.isfile(entropy_file):
            seed_entropy = np.loadtxt(entropy_file)
            seed_entropy = int(seed_entropy)
        else:
            sq = SeedSequence()
            seed_entropy = sq.entropy
            np.savetxt(entropy_file, [seed_entropy],
                       delimiter=',', fmt="%1.1i")
    else:
        # save user defined entropy anyway
        entropy_file_user_defined = base_folder + '/' + 'entropy_user_defined.txt'
        np.savetxt(entropy_file_user_defined, [
                   seed_entropy], delimiter=',', fmt="%1.1i")
    # generate random configurations
    coords_list, seeds = generate_random_configurations(
        box_length, n_part, n_dim, n_ensemble, seed_entropy)
    save_random_configurations(coords_list, seeds, base_folder)
    return coords_list, seeds


def load_configuration(foldpath, subfoldname, spawn_key):
    """
    loads configuration given the ensemble directory and spawn key 
    """
    random_config_filename = foldpath + '/' + \
        subfoldname + '/' + 'ensemble/' + str(spawn_key)
    return np.loadtxt(random_config_filename, delimiter=',')