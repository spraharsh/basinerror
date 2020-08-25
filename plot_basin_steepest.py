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
from map_basin_steepest import QUENCH_FOLDER_NAME
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
    res.xspace = ucoordsx[1] -ucoordsx[0]
    res.yspace = ucoordsy[1] - ucoordsy[0]
    res.xlen = len(ucoordsx)
    res.ylen = len(ucoordsy)
    return res

def _get_colors(num_colors):
    """ generate a unique colormap for use with matplotlib
    """
    colors_list=[]
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i/360.
        lightness = (50 + np.random.rand() * 10)/100.
        saturation = (90 + np.random.rand() * 10)/100.
        colors_list.append(colorsys.hls_to_rgb(hue, lightness, saturation))
        random.shuffle(colors_list)
    return colors.ListedColormap(colors_list)




if __name__ == "__main__":
    foldnameInversePower = "ndim=2phi=0.9seed=0n_part=8r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0"
    quench_type = QUENCH_FOLDER_NAME
    data_fold_path = BASE_DIRECTORY + '/' + foldnameInversePower + '/' + QUENCH_FOLDER_NAME
    data = CheckSameMinimum.load_map(data_fold_path, max_minima_l=2000)
    initial_coords = data.initial_coords
    order_params = data.order_params
    minimalist = data.minimalist
    res = extract_min_max_spacing(initial_coords)
    print(res.xspace, "xspace")
    print(res.yspace, "yspace")
    (hs_radii, initial_coords, box_length) = load_secondary_params(BASE_DIRECTORY
                                                                   + '/' +
                                                                   foldnameInversePower)
    op_2d = np.reshape(order_params, (res.xlen, res.ylen))
    print(op_2d)
    d = lambda x: x/box_length
    
    cmap = colors.ListedColormap(cc.glasbey_bw)
    # cmapsmall = 
    plt.imshow(op_2d,
               extent=[d(res.ymin), d(res.ymax),
                       d(res.xmin), d(res.xmax)],
               cmap=cmap)
    print(np.max(order_params))
    plt.xlabel('x (L)')
    plt.ylabel('y (L)')
    plt.title(quench_type + ' (L = length of box)')
    plt.savefig(BASE_DIRECTORY + '/' + foldnameInversePower + '/' + quench_type + '.pdf')
    plt.show()
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