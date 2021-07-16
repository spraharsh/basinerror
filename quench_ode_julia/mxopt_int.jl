include("../../basins.jl/src/optimizer/newton.jl")
include("../../basins.jl/src/minimumassign/mxopt.jl")


# natoms = 64
# radii_arr = generate_radii(0, natoms, 1.0, 1.4, 0.05, 0.05 * 1.4)
# dim = 2
# phi = 0.9
# power = 2.5
# eps = 1


# length_arr = get_box_length(radii_arr, phi, dim)

# coords = generate_random_coordinates(length_arr, natoms, dim)

# using PyCall
# boxvec = [length_arr, length_arr]

# utils = pyimport("pele.utils.cell_scale")
# cell_scale = utils.get_ncellsx_scale(radii_arr, boxvec)
# println(cell_scale)

# pele_wrapped_pot_2 = pot.InversePower(
#     2.5,
#     1.0,
#     radii_arr,
#     ndim = 2,
#     boxvec = boxvec,
#     use_cell_lists = false,
#     ncellx_scale = cell_scale,
# )


# pele_wrapped_python_pot_2 = PythonPotential(pele_wrapped_pot_2)

# using Sundials
# solver = CVODE_BDF()

# println(pele_wrapped_pot)
# println(coords)
# println(boxvec)
# println(radii_arr)
# mxd = Mixed_Descent(pele_wrapped_python_pot_2, solver, nls, coords,  30, 10^-5, 0.0 , 10^-3)
# println(coords)


# run!(mxd, 2000)

# println(coords)
# println(mxd.integrator.u .- coords)
# println("hellow world")
# println(p_energy(coords))



# function func2(x)
#    return 2*x
# end




