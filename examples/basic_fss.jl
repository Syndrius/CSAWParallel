"""

Basic use case of MIDParallel. This file is designed for running in the repl, not in parallel.
Gives example of seting up the problem, then problme should be solved outside the repl in parallel, then this file can read the outputs.

Probably want to change this to a different q-profile eventually.

"""

using MID 
using MIDParallel

using Plots; plotlyjs()


#first we define the problem and write to file.
#this is identical to MID.
Nr=100;
#rgrid = collect(LinRange(0, 1, N));
geo = GeoParamsT(R0=4.0);
prob = init_problem(q=fu_dam_q, geo=geo); 
rgrid = init_fem_grid(N=Nr)
θgrid = init_sm_grid(start=1, count = 2)
ζgrid = init_sm_grid(start=-1, count = 1)
grids = init_grids(rgrid, θgrid, ζgrid);


#looks like full path is needed... a bit annoying tbh.
#dir_base = "/home/149/mt3516/island_damping/MIDParallel/data/example/"
dir_base = "data/example/"

inputs_to_file(prob=prob, grids=grids, dir=dir_base);

"""
Now execute the command in parallel from terminal/bash script.
eg run from MIDParallel/
>>mpiexecjl -n 2 julia -e 'using MIDParallel; using MID; par_spectrum_from_file(dir="data/example/", target_freq=0.29)'

See convergence.sh for other examples of running in parallel.

"""

#now we can read the data in. first the eigenvalues,
evals = evals_from_file(dir=dir_base);

plot_continuum(evals);

display(evals.ω[1]);

ϕ1 = efunc_from_file(dir=dir_base, ind=1);
plot_potential(ϕ1, grids);



