"""
Generates an ensemble of initial points
for which we want to find the corresponding basins of attraction.
This file contains functuions for generating planes
"""

import numpy as np


def construct_point_set_nd(center_coords, nmesh_list, vec_list, n_scale_list, box_length):
    """
    Constructs a set of points in n dimensions around a central point.

    Parameters
    ----------
    nmesh : Int
        number of points in the meshed grid in each dimension
    vec_list : List of np.ndarray
        list of the vectors defining the mesh in each dimension
    n_scale_list : List of float
        list of the scale of the mesh in each dimension
    box_length : float
        length of the box. Necessary to ensure points stay within the box

    Returns
    -------
        set of points in n dimensions around a center_coords
    """

    
    
    # minimum_coords = np.loadtxt(foldpath + '/coords_of_minimum.txt',
    #                             delimiter=',')

    center = center_coords
    center = [-box_length/2, -box_length/2]
    vec_arr = np.array(vec_list)
    n_scale_arr = np.array(n_scale_list)

    scaled_vecs = n_scale_arr * vec_arr

    # check for orthogonality among vec_arr
    orthogonal_vecs = np.zeros(vec_arr.shape)
    for i in range(vec_arr.shape[0]):
        for j in range(vec_arr.shape[0]):
            if i != j:
                if np.dot(vec_arr[i], vec_arr[j]) != 0:
                    print(vec_arr[i], vec_arr[j])
                    raise ValueError("Vectors are not orthogonal")




    if (nmesh != 1):
        # construct a unit mesh for each dimension
        unit_mesh = construct_plain_centered_mesh(nmesh_list)

        # calculate displacement vectors for each mesh point
        displacement_vecs = np.einsum('ij,jk->ik', unit_mesh, scaled_vecs)

        # add displacement vectors to the center
        displacement_vecs = displacement_vecs + center
        return displacement_vecs
    else:
        return np.array([center_coords])



def construct_plain_centered_mesh(nmesh_list):
    """ Constructs a unit mesh of nmesh points in each dimension
        that is centered around the origin.
    Parameters
    ----------
    nmesh_list : list of int
        list of the number of points in the mesh in each dimension
    """
    mesh_points = []
    for i in range(len(nmesh_list)):
        mesh_points.append(np.linspace(0, 1, nmesh_list[i]))
    return np.array(np.meshgrid(*(np.array(mesh_points)- 0.5))).T.reshape(-1, len(nmesh_list))



if __name__ == '__main__':
    pass
