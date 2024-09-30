
#real juicy one!


using MID 
using MIDParallel

using Plots; plotlyjs()


#first we define the problem and write to file.
#this is identical to MID.
Nr=30;
Nθ=6;
Nζ=1;
#rgrid = collect(LinRange(0, 1, N));
geo = GeoParamsT(R0=10.0);
prob = init_problem(q=Axel_q, geo=geo); 
rgrid = rfem_grid(N=Nr)
θgrid = afem_grid(N=Nθ, pf=2);
ζgrid = afem_grid(N=Nζ, pf=-2);
grids = init_grids(rgrid, θgrid, ζgrid);
#tae_freq = 0.396 #/ geo.R0)^2; #previously found tae_freq.

#looks like full path is needed... a bit annoying tbh.
#dir_base = "/home/149/mt3516/island_damping/MIDParallel/data/example/"
dir_base = "data/example/"

#dir_base = "/scratch/y08/mt3516/fff/fu_dam/300x20x8/"

mkpath(dir_base)


inputs_to_file(prob=prob, grids=grids, dir=dir_base);

"""
Now execute the command in parallel from terminal/bash script.
eg run from MIDParallel/
>>mpiexecjl -n 2 julia -e 'using MIDParallel; using MID; par_spectrum_from_file(dir="data/example/", target_freq=0.31)'
>>mpiexec -n 2 julia -e 'using MIDParallel; using MID; par_spectrum_from_file(dir="/home/149/mt3516/island_damping/MIDParallel/data/example/", target_freq=0.301)'

See convergence.sh for other examples of running in parallel.

"""

process_hdf5(dir_base)

#now we can read the data in. first the eigenvalues,
evals = evals_from_file(dir=dir_base);

#wot is this, this is cooked af
plot_continuum(evals);

#found the tae, freq is significantly different for this low res example.
tae_ind = find_ind(evals, 0.3835)

display(length(evals.ω))

tae_ind = find_ind(evals, 0.764298)


#the tae for this small res is ind=14.
ϕft = efunc_from_file(dir = dir_base, ind=tae_ind);
#potential plotting is fked af.
plot_potential(ϕft, grids);

ϕ = efunc_from_file(dir = dir_base, ind=tae_ind, ft=false);
#potential plotting is fked af.
plot_potential(ϕ, grids);

#looks good, other than mode labels.
contour_plot(ϕ, grids)

surface_plot(ϕ, grids)


prob, grids = inputs_from_file(dir=dir_base);

display(grids)