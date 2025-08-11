"""

Basic use case of MIDParallel. This file is designed for running in the repl, not in parallel.
Gives example of seting up the problem, then problme should be solved outside the repl in parallel, then this file can read the outputs.

Probably want to change this to a different q-profile eventually. 

"""

using MID 
using MIDParallel
using MIDViz
using Plots; plotlyjs()
#%%

#first we define the problem and write to file.
#this is identical to MID.
#rgrid = collect(LinRange(0, 1, N));
geo = init_geo(R0=4.0);
prob = init_problem(q=fu_dam_q, geo=geo); 
#%%

Nr=20;
rgrid = init_grid(type=:rf, N=Nr)
θgrid = init_grid(type=:as, start=1, N = 2)
ζgrid = init_grid(type=:as, start=-1, N = 1)
grids = init_grids(rgrid, θgrid, ζgrid);
#%%

#solver = init_solver(nev=100, targets=[0.0, 0.33, 0.8], prob=prob)
#solver = init_solver(nev=3, target=0.33, prob=prob)
solver = init_solver(nev=10, targets=[0.2, 0.3], prob=prob)


#%%

#looks like full path is needed... a bit annoying tbh.
dir_base = "/scratch/y08/mt3516/test/"
dir_base = "/Users/matt/phd/MIDParallel/data/example/"
#%%

inputs_to_file(prob=prob, grids=grids, solver=solver, dir=dir_base);

"""
Now execute the command in parallel from terminal/bash script.
eg run from MIDParallel/
>>mpiexec -n 2 julia -e 'using MIDParallel; par_spectrum_from_file(dir="/Users/matt/phd/MIDParallel/data/example/")'
>>mpiexec -n 2 julia -e 'using MIDParallel; par_spectrum_from_file(dir="/scratch/y08/mt3516/test/")'

See convergence.sh for other examples of running in parallel.

"""



par_post_process(dir_base) #unfort have we have to do this!
#now we can read the data in. first the eigenvalues,
evals = evals_from_file(dir=dir_base);


continuum_plot(evals);


ind = find_ind(evals, 0.289)


ϕft = efunc_from_file(dir=dir_base, ind=ind);
potential_plot(ϕft, grids);
ϕ = efunc_from_file(dir=dir_base, ind=1, ft=false);
#plot_potential(ϕ, grids);

#inconsistent with phase of two modes, here we have the og case for some reason...
contour_plot(ϕ, grids, ind=1)
surface_plot(ϕ, grids, ind=1)
