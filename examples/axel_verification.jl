


using MID
using MIDParallel

using Plots; plotlyjs()


#first we define the problem and write to file.
#this is identical to MID.
#see what kind of mem is required for this!
#think the matrix is v sparse now, so hopefully memory is not as big a problemo.
Nr=500;
Nθ=20
Nζ=3
#rgrid = collect(LinRange(0, 1, N));
geo = GeoParamsT(R0=10.0);
flr = FLRT(δ = -4.0e-7)
prob = init_problem(q=Axel_q, geo=geo, dens=axel_dens, flr=flr)
rgrid = rfem_grid(N=Nr, sep1=0.8, sep2=0.95, frac=0.4)
θgrid = afem_grid(N=Nθ, pf=2)
#ζgrid = init_sm_grid(start=-2, count = 1)
ζgrid = afem_grid(N=Nζ, pf=-2)
grids = init_grids(rgrid, θgrid, ζgrid);


#looks like full path is needed... a bit annoying tbh.
#dir_base = "/home/149/mt3516/island_damping/MIDParallel/data/example/"
#dir_base = "/scratch/y08/mt3516/ffs_verification/axel/"

dir = "/scratch/y08/mt3516/fff_axel_damping/"

#mkpath(dir_base)

#dir = dir_base * "50x10/"

mkpath(dir)

inputs_to_file(prob=prob, grids=grids, dir=dir);



prob, grids = inputs_from_file(dir=dir);

process_hdf5(dir)

#now we can read the data in. first the eigenvalues,
evals = evals_from_file(dir=dir);

continuum_plot(evals)

tae_ind = find_ind(evals, 0.394)

display(evals.ω[tae_ind])