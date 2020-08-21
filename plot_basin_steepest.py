import numpy as np
import matplotlib.pyplot as plt
from params import BASE_DIRECTORY
from params import load_secondary_params












if __name__ == "__main__":
    foldnameInversePower = "ndim=2phi=0.9seed=0n_part=8r1=1.0r2=1.4rstd1=0.05rstd2=0.06999999999999999use_cell_lists=0power=2.5eps=1.0"
    # boollist = np.loadtxt(BASE_DIRECTORY + '/' + foldnameInversePower + 'quench_results.txt')
    boollist = np.loadtxt(BASE_DIRECTORY + '/'
                          + foldnameInversePower
                          + '/' + 'quench_results_mxopt.txt')
    yrange = np.loadtxt(BASE_DIRECTORY + '/' + foldnameInversePower + 'yrange.txt', delimiter = ',')
    (hs_radii, initial_coords, box_length) = load_secondary_params(BASE_DIRECTORY + '/' + foldnameInversePower)
    box_length = float(box_length)
    boollistreshaped = np.reshape(boollist, (21, 21))
    print(boollistreshaped)
    plt.imshow(boollistreshaped, extent= (yrange[0]/box_length, yrange[-1]/box_length, yrange[0]/box_length, yrange[-1]/box_length))
    plt.xlabel('x (L)')
    plt.ylabel('y (L)')
    plt.title('Mixed optimizer convtol= 10000 (L = length of box)')
    plt.savefig(BASE_DIRECTORY + '/' + foldnameInversePower + '/' + 'mixed_optimizerconvtol10000.pdf')
    plt.show()