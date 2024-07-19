


using MID
using MIDParallel

using Plots; plotlyjs()


#first we define the problem and write to file.
#this is identical to MID.
#see what kind of mem is required for this!
#think the matrix is v sparse now, so hopefully memory is not as big a problemo.
Nr=50;
Nθ=10
#rgrid = collect(LinRange(0, 1, N));
geo = GeoParamsT(R0=10.0);
prob = init_problem(q=Axel_q, geo=geo, dens=axel_dens, δ=-0e-7)#, met=diagonal_toroidal_metric!); 
rgrid = init_fem_grid(N=Nr)
θgrid = init_fem_grid(N=Nθ, pf=2)
ζgrid = init_sm_grid(start=-2, count = 1)
grids = init_grids(rgrid, θgrid, ζgrid);


#looks like full path is needed... a bit annoying tbh.
#dir_base = "/home/149/mt3516/island_damping/MIDParallel/data/example/"
dir_base = "/scratch/y08/mt3516/ffs_verification/axel/"

mkpath(dir_base)

dir = dir_base * "50x10/"

mkdir(dir)

inputs_to_file(prob=prob, grids=grids, dir=dir);



ω = par_vals_from_file(dir*"vals.dat", geo.R0);

#with diagonal met 800x50
#0.38761929918593957 - 0.0013047632651995039im 
#with normal met 800x50
#0.39439328047979666 - 0.001752252913284776im

display(ω[1])

display(imag(ω[1]^2))
(0.38761929918593957 - 0.0013047632651995039im)^2

#then the efuncs
#ϕ = par_func_from_file(dir*"funcs.dat", 1, grids);
ϕ = par_funcs_from_file(dir*"funcs.dat", length(ω), grids);

ϕms = mode_structure(ϕ, grids);

reconstruct_continuum(ω, ϕms, grids)

#maybe we want an imag flag for this? i.e to match fig 5 of Bowden.
plot_potential(ϕms, grids, 1, 1)



plot_phi_surface(ϕ, grids, 1)