"""
Script to plot minima configurations from file
"""

# from map_basin_steepest import QUENCH_FOLDER_NAME, MINIMA_DATABASE_NAME



from ssl import ALERT_DESCRIPTION_BAD_CERTIFICATE
from map_basin_steepest import MINIMA_DATABASE_NAME
import numpy as np
import matplotlib.pyplot as plt
from params import load_params, load_secondary_params


import colorcet as cc
# for overlap intersections
import shapely.geometry as sg
import descartes

glasbey = cc.glasbey





#                            global variables: rewrite

# Database path for radii



BASE_DIRECTORY =  "/home/praharsh/Dropbox/research/bv-libraries/basinerror/datainv/ndim=2phi=0.9seed=0n_part=32r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0/"
MINIMA_DATABASE_PATH = BASE_DIRECTORY + 'minima_database.npy'




# reshape for cycling by particle




def get_configuration(folder_path, minima_database_name, minima_number):
    """ Gets a minimum configuration i.e the coordinates and the radii

    Args:
       folder_path path
       minima_database_name folder name
       minima_number index of minima to be plotted

    Returns:
        minimum configuration, radii
    """
    mimima_database_path = folder_path + minima_database_name
    minima = np.load(mimima_database_path)
    (radii, initial_coords, box_length) = load_secondary_params(folder_path)
    nparticles = len(radii)
    dimensions = 2
    configuration = np.reshape(minima[minima_number], (nparticles, dimensions))
    return radii, configuration, box_length


    


def plot_configuration_2d(radii, configuration, box_length, figname):
    """ Plots a given minima configuration in 2d as defined by the minimum number and saves it as pdf
    Args:
     radii
     configuration
     box length
     figname
    """
    dimensions=2
    nparticles= len(radii)
    assert(nparticles == (len(configuration)))
    
    # take columnwise mean. gives mean displacement
    # puts in box
    configuration = configuration % box_length
    # center of mass assuming the masses are the same
    # calculated by taking the mean along the column axis
    center_of_mass = np.mean(configuration, axis=0)
    
    fig, ax = plt.subplots(figsize=(5,5))
    for particle in range(nparticles):
        # non dimensonalized w.r.t box length
        radius = radii[particle]/box_length
        coordinates = (configuration[particle] -center_of_mass)/box_length + 0.5
        
        # discs
        disc = sg.Point(coordinates[0],coordinates[1]).buffer(radius)

        # add images for periodic boundary conditions
        disc_up = sg.Point(coordinates[0],coordinates[1]+1.).buffer(radius)
        disc_down = sg.Point(coordinates[0],coordinates[1]-1.).buffer(radius)
        disc_right = sg.Point(coordinates[0]+1,coordinates[1]).buffer(radius)
        disc_left = sg.Point(coordinates[0]-1,coordinates[1]).buffer(radius)
        disc_left_up = sg.Point(coordinates[0]-1.,coordinates[1]+1.).buffer(radius)
        disc_left_down = sg.Point(coordinates[0]-1.,coordinates[1]-1.).buffer(radius)
        disc_right_up = sg.Point(coordinates[0]+1,coordinates[1]+1.).buffer(radius)
        disc_right_down = sg.Point(coordinates[0]+1,coordinates[1]-1.).buffer(radius)
        
        alpha = 0.3
        disc_color = glasbey[particle]
        border_color = 'k'
        ax.add_patch(descartes.PolygonPatch(disc, fc =disc_color, ec= border_color, alpha=alpha))
        ax.add_patch(descartes.PolygonPatch(disc_up, fc =disc_color, ec= border_color, alpha=alpha))
        ax.add_patch(descartes.PolygonPatch(disc_down, fc =disc_color, ec= border_color, alpha=alpha))
        ax.add_patch(descartes.PolygonPatch(disc_right, fc =disc_color, ec= border_color, alpha=alpha))
        ax.add_patch(descartes.PolygonPatch(disc_left, fc =disc_color, ec= border_color, alpha=alpha))
        ax.add_patch(descartes.PolygonPatch(disc_left_up, fc =disc_color, ec= border_color, alpha=alpha))
        ax.add_patch(descartes.PolygonPatch(disc_left_down, fc =disc_color, ec= border_color, alpha=alpha))
        ax.add_patch(descartes.PolygonPatch(disc_right_up, fc =disc_color, ec= border_color, alpha=alpha))
        ax.add_patch(descartes.PolygonPatch(disc_right_down, fc =disc_color, ec= border_color, alpha=alpha))

    plt.axis('off')
    
    plt.savefig(figname + '.svg', dpi=300, bbox_inches='tight',pad_inches = 0)
    plt.show()
    
    


if __name__=="__main__":
    mdn = 'minima_database.npy'
    minima_number = 0
    minima_8 = [14, 42, 7]
    minima_16 = [106]
    minima_32 = [30, 47]

    BASE_DIRECTORY_8 = "/home/praharsh/Dropbox/research/bv-libraries/basinerror/datainv/ndim=2phi=0.9seed=0n_part=8r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0/"
    configuration_list = []
    for minima_number in minima_8:
        radii, configuration, box_length = get_configuration(BASE_DIRECTORY_8, mdn, minima_number)
        configuration_list.append(configuration)
        print(radii)
        plot_configuration_2d(radii, configuration, box_length, 'configuration' + str(minima_number) + 'n8fixed')

    BASE_DIRECTORY_16 = "/home/praharsh/Dropbox/research/bv-libraries/basinerror/datainv/ndim=2phi=0.9seed=0n_part=16r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0/"
    configuration_list = []
    for minima_number in minima_16:
        radii, configuration, box_length = get_configuration(BASE_DIRECTORY_16, mdn, minima_number)
        configuration_list.append(configuration)
        print(radii)
        plot_configuration_2d(radii, configuration, box_length, 'configuration' + str(minima_number) + 'n16fixed')


    BASE_DIRECTORY_32 = "/home/praharsh/Dropbox/research/bv-libraries/basinerror/datainv/ndim=2phi=0.9seed=0n_part=32r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0/"
    configuration_list = []
    for minima_number in minima_32:
        radii, configuration, box_length = get_configuration(BASE_DIRECTORY_32, mdn, minima_number)
        configuration_list.append(configuration)
        print(radii)
        plot_configuration_2d(radii, configuration, box_length, 'configuration' + str(minima_number) + 'n32fixed')
    
    



