#finding additional information about the harmonic structure of the computed spectrum
using MID 
using MIDParallel
using MIDViz
#using Plots; plotlyjs()
#%%

#first we define the problem and write to file.
#this is identical to MID.
#rgrid = collect(LinRange(0, 1, N));
geo = init_geo(R0=4.0);
prob = init_problem(q=fu_dam_q, geo=geo); 
#%%

Nr=100;
rgrid = init_grid(type=:rf, N=Nr)
θgrid = init_grid(type=:as, start=1, N = 3)
ζgrid = init_grid(type=:as, start=-2, N = 2)
grids = init_grids(rgrid, θgrid, ζgrid);
#%%

#solver = init_solver(nev=100, targets=[0.0, 0.33, 0.8], prob=prob)
#solver = init_solver(nev=3, target=0.33, prob=prob)
solver = init_solver(nev=80, targets=[0.2, 0.3, 0.4], prob=prob)


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
#this assumes post_processing has been done for evals.jld2
#this is pretty slow, perhaps expected.
MID.Mapping.harmonic_info(dir_base)
#%%
using JLD2
using Plots
ht = load_object(joinpath(dir_base, "ht.jld2"))

scatter(ht.x1, ht.fwhm)
scatter(ht.x1, ht.harm2)
length(ht.x1)


