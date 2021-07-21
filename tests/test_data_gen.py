from basinerror import generate_plane


def test_generate_plane_mesh_3d():
    Vec_x = [1, 0, 0]
    Vec_y = [0, 1, 0]
    Vec_z = [0, 0, 1]
    vec_list = [Vec_x, Vec_y, Vec_z]
    nmesh_list = [3, 3, 1]
    point = [0.5, 0.5, 0.5]
    n_scale_list = [0.1, 0.1, 0.1]
    points = generate_plane.construct_point_set_nd(
        point, nmesh_list, vec_list, n_scale_list, 1)
    # just check the shape of the output
    assert points.shape == (9, 3)
    assert points[0, 0] == 0.45
    assert points[8, 2] == 0.5

def test_generate_plane_mesh_3d_2_vecs():
    Vec_x = [1, 0, 0]
    Vec_y = [0, 1, 0]
    Vec_z = [0, 0, 1]
    vec_list = [Vec_x, Vec_y]
    nmesh_list = [3, 3]
    point = [0.5, 0.5, 0.5]
    n_scale_list = [0.1, 0.1, 0.1]
    points = generate_plane.construct_point_set_nd(
        point, nmesh_list, vec_list, n_scale_list, 1)
    assert points.shape == (9, 3)
    assert points[0, 0] == 0.45
    assert points[8, 2] == 0.5


if __name__ == '__main__':
    test_generate_plane_mesh_3d()
    test_generate_plane_mesh_3d_2_vecs()
