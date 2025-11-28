"""

Basic use case of MIDParallel. This file is designed for running in the repl, not in parallel.
Gives example of seting up the problem, then problem should be solved outside the repl in parallel, then this file can read the outputs.

"""

using MID 
using MIDParallel
#using MIDViz
using Plots
#%%

#first we define the problem and write to file.
#this is identical to MID.
#rgrid = collect(LinRange(0, 1, N));
geo = init_geometry()#, met)
fields = init_fields()

prob = init_problem(geometry=geo, fields=fields)
#%%
x1grid = init_grid(:ψ, 30)
x2grid = init_grid(:sm, 2, start=1)
x3grid = init_grid(:spectral, 1, start=-1) #perhaps remove start from kwargs?

x1grid = init_grid(:ψ, 40)
x2grid = init_grid(:θ, 10, pf=1)
x3grid = init_grid(:spectral, 1, start=-1) #perhaps remove start from kwargs?

x1grid = init_grid(:ψ, 30)
x2grid = init_grid(:θ, 6, pf=1)
x3grid = init_grid(:φ, 1, pf=-1) #perhaps remove start from kwargs?

grids = init_grids(x1grid, x2grid, x3grid)

#%%

#solver = init_solver(nev=100, targets=[0.0, 0.33, 0.8], prob=prob)
solver = init_solver(nev=100, target=0.3, prob=prob)


#%%

#looks like full path is needed... a bit annoying tbh.
#dir_base = "/scratch/y08/mt3516/test/"
dir_base = "/Users/matt/phd/MIDParallel/data/example/"
par_compute_spectrum(dir_base)
#%%

inputs_to_file(prob=prob, grids=grids, solver=solver, dir=dir_base);

"""
#need to change this, instead we will probably just use the run script.
Now execute the command in parallel from terminal/bash script.
eg run from MIDParallel/
>>mpiexec -n 2 julia -e 'using MIDParallel; par_spectrum_from_file(dir="/Users/matt/phd/MIDParallel/data/example/")'
>>mpiexec -n 2 julia -e 'using MIDParallel; par_spectrum_from_file(dir="/scratch/y08/mt3516/test/")'

See convergence.sh for other examples of running in parallel.

"""



par_post_process(dir_base) #unfort have we have to do this!
#now we can read the data in. first the eigenvalues,
evals = evals_from_file(dir_base);


continuum_plot(evals);

scatter(evals.x1, real.(evals.ω))

ind = find_ind(evals, 0.3)


ϕft = efunc_from_file(dir_base, ind);

pgrid = MID.inst_grid(x1grid)
#pretty fkn annoying tbh!
plot(pgrid, real.(ϕft[:, 1, 1]))
plot!(pgrid, real.(ϕft[:, 2, 1]))
plot!(pgrid, real.(ϕft[:, 3, 1]))


harmonic_plot(ϕft, grids);
ϕ = efunc_from_file(dir=dir_base, ind=1, ft=false);

contour_plot(ϕ, grids, ind=1)
surface_plot(ϕ, grids, ind=1)
