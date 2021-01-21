"""
Generates a set of starting random configurations
and stores them in a folder with unique name after the parameters.

The current version generates a set of random configurations in 2d and saves them to file
"""



from numpy.core.fromnumeric import size
from utils.cell_scale import get_box_length
from params import BASE_DIRECTORY, load_params
import numpy as np
from numpy.random import Generator, PCG64, SeedSequence
import os

SUB_FOLD_NAME= 'random_configuration_run'

def generate_random_radii_inverse_power(params, seed_radii, foldpath, subfoldname=SUB_FOLD_NAME):
    """ Generates radii for inverse potential and saves accordingly
    Args:
       foldname folder name
       seed seed
    Returns:
        Points
    """
    n_part_by_2 = params.n_part.value //2

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
    np.savetxt(foldpath + '/' + subfoldname + '/radii.txt', radii)
    return radii

def get_hs_radii(foldpath, subfoldname, delimiter=','):
    return np.loadtxt(foldpath + '/' + subfoldname + '/radii.txt')



def generate_random_configuration_single(box_length,n_part, n_dim, seed_coords):
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



def generate_random_configurations(foldpath, ensemble_size, seed_entropy=None, seed_radii=0, subfoldname=SUB_FOLD_NAME):
    """ generates and saves random configurations 
        for particles given a folder with particular configurations
    
    Args:
       foldname ${1:arg1}
       seed_radii: seed for radii
       seed_coords_list: seeds list for random configuration if none generates from seedsequence
       subfoldname ${2:arg2} (default 'random_configurations_run')
    """

    # print(foldpath)
    # print(str(foldpath)+'/systemdata.json')
    foldpath=str(foldpath)
    params = load_params(foldpath)
    radii = generate_random_radii_inverse_power(params, seed_radii, foldpath, subfoldname=SUB_FOLD_NAME)
    box_length = get_box_length(radii, params.ndim.value, params.phi.value)

    run_foldpath= foldpath + '/' + subfoldname + '/' 
    os.makedirs(run_foldpath, exist_ok=True)
    entropy_file = run_foldpath+ '/' + 'seed_entropy.txt'
    
    if seed_entropy==None:
        # check whether entropy already exists in the folder
        previous_entropy_exists = os.path.isfile(entropy_file)
        
        if previous_entropy_exists:
            # hacky way to deal with large integer issues with seeds
            seed_entropy = np.loadtxt(entropy_file, delimiter=',',dtype=str)
            seed_entropy = int(seed_entropy)
        else:
            # generate entropy and save to file
            sq = SeedSequence()
            seed_entropy = sq.entropy
            np.savetxt(run_foldpath+ '/' + 'seed_entropy.txt', [seed_entropy], delimiter=',', fmt="%1.1i")
    else:
        # save user defined entropy/seed to file any ways
        np.savetxt(run_foldpath+ '/' + 'seed_entropy_user_defined.txt', [seed_entropy], delimiter=',', fmt="%1.64i")

    print(seed_entropy)
    sq1 = SeedSequence(seed_entropy)
    
    seeds = sq1.spawn(ensemble_size)
    assert(len(seeds) == ensemble_size)
    print(seed_entropy)
    
    config_list = []
    for seed_coords in seeds:
        config_list.append(generate_random_configuration_single(box_length, params.n_part.value, params.ndim.value, seed_coords))
    return np.array(config_list)


def save_random_configurations(random_configs, foldpath, subfoldname=SUB_FOLD_NAME):
    """ saves random configurations to a subfolder as a numpy array
    Args:
       random_configs
    """
    ensemblesize = len(random_configs)
    np.savetxt(foldpath+'/' + subfoldname +'/'
               + random_config_filename(ensemble_size), random_configs, delimiter=',')


def random_config_filename(ensemble_size):
    return  'random_configs' + str(ensemble_size) + '.txt'


def load_configurations(foldpath, ensemble_size, subfoldname=SUB_FOLD_NAME):
    return np.loadtxt(foldpath+'/' + subfoldname +'/'
                      + random_config_filename(ensemble_size), delimiter=',')







if __name__=="__main__":
    foldname = "ndim=2phi=0.9seed=0n_part=8r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0"
    foldpath = str(BASE_DIRECTORY+'/' + foldname)
    ensemble_size = int(5e4)
    configs = generate_random_configurations(foldpath, ensemble_size)
    save_random_configurations(configs, foldpath)
    configs2 =load_configurations(foldpath, ensemble_size)












