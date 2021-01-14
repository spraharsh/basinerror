"""
Script to plot minima configurations from file
"""

# from map_basin_steepest import QUENCH_FOLDER_NAME, MINIMA_DATABASE_NAME



from ssl import ALERT_DESCRIPTION_BAD_CERTIFICATE
from map_basin_steepest import MINIMA_DATABASE_NAME
import numpy as np
import matplotlib.pyplot as plt
from params import load_params, load_secondary_params

# for overlap intersections
import shapely.geometry as sg
import descartes





#                            global variables: rewrite

# Database path for radii



BASE_DIRECTORY =  "/home/praharsh/Dropbox/research/bv-libraries/basinerror/datainv/ndim=2phi=0.9seed=0n_part=16r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0/"
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
    print(radii)
    print(configuration)
    print(len(radii))
    print(len(configuration))
    assert(nparticles == (len(configuration)))

    fig, ax = plt.subplots(figsize=(5,5))
    for particle in range(nparticles):
        # non dimensonalized w.r.t box length
        radius = radii[particle]/box_length
        coordinates = configuration[particle]/box_length
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
        disc_color = 'b'
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

    plt.xlabel(r'x($L$)')
    plt.ylabel(r'y($L$)')
    plt.savefig(figname + '.pdf')
    plt.show()
    
    


if __name__=="__main__":
    mdn = 'minima_database.npy'
    minima_number = 0
    radii, configuration, box_length = get_configuration(BASE_DIRECTORY, mdn, minima_number)
    plot_configuration_2d(radii, configuration, box_length, 'configuration' + str(minima_number))




