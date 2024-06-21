"""

Basic use case of MIDParallel. This file is designed for running in the repl, not in parallel.
Gives example of seting up the problem, then problme should be solved outside the repl in parallel, then this file can read the outputs.

Probably want to change this to a different q-profile eventually.

"""

using MID 
using MIDParallel


#first we define the problem and write to file.
#this is identical to MID.
N=100;
#rgrid = collect(LinRange(0, 1, N));
geo = GeoParamsT(R0=10.0);
prob = init_problem(q=Axel_q, geo=geo); 
grids = init_grids(N=N, mstart=2, mcount=2, nstart=-2, ncount=1);
tae_freq = (0.381 / geo.R0)^2; #previously found tae_freq.

#looks like full path is needed... a bit annoying tbh.
#dir_base = "/home/149/mt3516/island_damping/MIDParallel/data/example/"
dir_base = "data/example/"

inputs_to_file(prob=prob, grids=grids, dir=dir_base);

"""
Now execute the command in parallel from terminal/bash script.
eg run from MIDParallel/
>>mpiexecjl -n 2 julia -e 'using MIDParallel; using MID; par_spectrum_from_file(dir="data/example/", freq=0.00145)'

See convergence.sh for other examples of running in parallel.

"""

#now we can read the data in. first the eigenvalues,
ω = par_vals_from_file(dir_base*"vals.dat", geo.R0);

#then the efuncs
ϕ = par_funcs_from_file(dir_base*"funcs.dat", length(ω), N, grids.pmd.count, grids.tmd.count);

#or we can read a single eigenfunction, which is useful for larger datasets, here we read the third eigenvalue.
ϕ3 = par_func_from_file(dir_base*"funcs.dat", 3, N, grids.pmd.count, grids.tmd.count);

#now we can plot the TAE, which will be the first efuncs as we have specified the desired frequency to solve for.
display(ω[1]); #should find a tae with normalised frequency 0.3812
plot_potential(grids=grids, ϕ=ϕ, ind=1, n=1);

#the individual eigenfunction only contains a single index so we pass in ind=1
plot_potential(grids=grids, ϕ=ϕ3, ind=1, n=1);

