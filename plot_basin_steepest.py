from tokenize import Double
from types import SimpleNamespace
from matplotlib.pyplot import show
from basinerror import quench_steepest
import numpy as np
import matplotlib.pyplot as plt
from params import BASE_DIRECTORY
from params import load_secondary_params
from checksameminimum import CheckSameMinimum
from types import SimpleNamespace
from matplotlib import colors
import colorsys
import random
import colorcet as cc
from map_basin_steepest import QUENCH_FOLDER_NAME, MINIMA_DATABASE_NAME
from accuracy_check import correct_minima_check
np.random.seed(0)


def extract_min_max_spacing(coordslist):
    """ Extracts spacing data from a 2d coordinate list
    Parameters
    ----------
    coordslist: list(points)
        List of points in a grid
    Returns
    ----------
    Namespace with ymin ymax, xmin, xmax and range
    """
    # extract and format : should be slow, choose better in
    # case it's too slow
    coords_x, coords_y = coordslist.T
    ucoordsx, ucoordsy = np.unique(coords_x), np.unique(coords_y)
    print(ucoordsy)
    print(ucoordsx)
    print(ucoordsx.shape)
    print(ucoordsy.shape)
    res = SimpleNamespace()
    res.xmin = np.min(ucoordsx)
    res.xmax = np.max(ucoordsx)
    res.ymin = np.min(ucoordsy)
    res.ymax = np.max(ucoordsy)
    res.xspace = ucoordsx[1] - ucoordsx[0]
    res.yspace = ucoordsy[1] - ucoordsy[0]
    res.xlen = len(ucoordsx)
    res.ylen = len(ucoordsy)
    return res


def _get_colors(num_colors):
    """ generate a unique colormap for use with matplotlib
    """
    colors_list = []
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i / 360.
        lightness = (50 + np.random.rand() * 10) / 100.
        saturation = (90 + np.random.rand() * 10) / 100.
        colors_list.append(colorsys.hls_to_rgb(hue, lightness, saturation))
        random.shuffle(colors_list)
    return colors.LinearSegmentedColormap(colors_list)


if __name__ == "__main__":

    # load data
    foldnameInversePower = "ndim=2phi=0.9seed=0n_part=32r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0"
    minima_database_path = BASE_DIRECTORY + '/' + foldnameInversePower + '/' + MINIMA_DATABASE_NAME
    # quench_type = "cvodeopt"
    # quench_type = "Fire"
    quench_type = QUENCH_FOLDER_NAME # if you want to plot the last one you get from map_basin_steepest
    # quench_type = 'correct_minima'
    data_fold_path = BASE_DIRECTORY + '/' + foldnameInversePower + '/' + quench_type
    print('loading data from')
    print(data_fold_path)

    data = CheckSameMinimum.load_map(data_fold_path,
                                     max_minima_l=2000,
                                     minima_database_path=minima_database_path)
    initial_coords = data.initial_coords
    order_params = data.order_params
    minimalist = data.minimalist
    

    # set min max values
    vmax = 255
    vmin = 0
    res = extract_min_max_spacing(initial_coords)
    
    xlen = int(np.sqrt(len(order_params)))
    (hs_radii, initial_coords,
     box_length) = load_secondary_params(BASE_DIRECTORY + '/' +
                                         foldnameInversePower)
    op_2d = np.reshape(order_params, (xlen, xlen))
    print(np.max(order_params), 'order parameters')
    # np.set_printoptions(threshold=np.inf)
    d = lambda x: x / box_length
    print(op_2d)
    cmap = colors.ListedColormap(cc.glasbey)
    print(len(cc.glasbey_bw), 'length of glasbey')
    cmap2 = 'tab20c'
    # cmapsmall =
    plt.imshow(op_2d,
               extent=[d(res.ymin),
                       d(res.ymax),
                       d(res.xmin),
                       d(res.xmax)],
               cmap=cmap,
               vmin=vmin,
               vmax=vmax)
    print(np.max(order_params))
    plt.xlabel('x (L)')
    plt.ylabel('y (L)')
    plt.title(QUENCH_FOLDER_NAME + ' (L = length of box)')
    plt.savefig(BASE_DIRECTORY + '/' + foldnameInversePower + '/' +
                quench_type + '.pdf')
    plt.show()
    print(correct_minima_check(quench_type, foldnameInversePower), 'correct_minima_heavy')
# if __name__ == "__main__":
#     ncolors = 30
#     gradient = np.linspace(0, 1, ncolors)
#     cmap = _get_colors(ncolors)
#     gradient = np.vstack((gradient, gradient))
#     plt.imshow(gradient, cmap=cmap)
#     plt.show()

# deprecated: old save data
# if __name__ == "__main__":
#     foldnameInversePower = "ndim=2phi=0.9seed=0n_part=8r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0"
#     # boollist = np.loadtxt(BASE_DIRECTORY + '/' + foldnameInversePower + 'quench_results.txt')
#     boollist = np.loadtxt(BASE_DIRECTORY + '/'
#                           + foldnameInversePower
#                           + '/' + 'quench_results_mxopt.txt')
#     yrange = np.loadtxt(BASE_DIRECTORY + '/' + foldnameInversePower + 'yrange.txt', delimiter = ',')
#     (hs_radii, initial_coords, box_length) = load_secondary_params(BASE_DIRECTORY + '/' + foldnameInversePower)
#     box_length = float(box_length)
#     boollistreshaped = np.reshape(boollist, (21, 21))
#     print(boollistreshaped)
#     plt.imshow(boollistreshaped, extent= (yrange[0]/box_length, yrange[-1]/box_length, yrange[0]/box_length, yrange[-1]/box_length))
#     plt.show()
