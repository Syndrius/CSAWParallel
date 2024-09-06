


using MID 
using MIDParallel
using Plots; plotlyjs()


#first we define the problem and write to file.
#this is identical to MID.
Nr=40;
Nθ=6;
#rgrid = collect(LinRange(0, 1, N));
geo = GeoParamsT(R0=10.0);
prob = init_problem(q=Axel_q, geo=geo); 
rgrid = init_fem_grid(N=Nr);
θgrid = init_fem_grid(N=Nθ, pf=2);
ζgrid = init_sm_grid(start=-2, count = 1);
grids = init_grids(rgrid, θgrid, ζgrid);
#tae_freq = 0.396 #/ geo.R0)^2; #previously found tae_freq.

#looks like full path is needed... a bit annoying tbh.
dir_base = "/home/149/mt3516/island_damping/MIDParallel/data/example/"
#dir_base = "data/example/"

inputs_to_file(prob=prob, grids=grids, dir=dir_base);

"""
Now execute the command in parallel from terminal/bash script.
eg run from MIDParallel/
>>mpiexecjl -n 2 julia -e 'using MIDParallel; using MID; par_spectrum_from_file(dir="data/example/", target_freq=0.31)'

See convergence.sh for other examples of running in parallel.

"""


process_hdf5(dir_base) #unfort have we have to do this!
#now we can read the data in. first the eigenvalues,
evals = evals_from_file(dir=dir_base);


plot_continuum(evals);

tae_ind = find_ind(evals, 0.7039555)


#mode labels are cooked af again... why. continuum seems to be fine.
#not plot_potential though..
#perhaps it is both??? V unsure though.
ϕft = efunc_from_file(dir=dir_base, ind=tae_ind);
plot_potential(ϕft, grids);
ϕ = efunc_from_file(dir=dir_base, ind=tae_ind, ft=false);

#so this seems to have a heap of θ lines of z=0... v worrying...
contour_plot(ϕ, grids, ind=1)
surface_plot(ϕ, grids, ind=1)