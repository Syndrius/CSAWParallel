
using MID 
using MIDParallel
using MIDViz
#%%

Nr=10;
Nθ=5;
Nζ=3;
#rgrid = collect(LinRange(0, 1, N));
geo = init_geo(R0=10.0);
prob = init_problem(q=fu_dam_q, geo=geo); 

isl = init_island(m0=2, n0=-1, w=0.03, r0=0.5, qp=2.0)
#prob = MID.Structures.init_isl_problem(geo=geo, isl=isl)
prob =  init_problem(geo=geo, q = cantori_q)

rgrid = init_grid(type=:rf, N=Nr, stop=0.999)
θgrid = init_grid(type=:af, N=Nθ, pf=0)
ζgrid = init_grid(type=:af, N=Nζ, pf=0)

grids = init_grids(rgrid, θgrid, ζgrid);
#%%
#tae_freq = 0.396 #/ geo.R0)^2; #previously found tae_freq.

solver = init_solver(nev=50, target=0.3, prob=prob)
#%%
#looks like full path is needed... a bit annoying tbh.
#dir_base = "/home/149/mt3516/island_damping/MIDParallel/data/example/"
dir_base = "/scratch/y08/mt3516/test/"
#dir_base = "/Users/matt/phd/MIDParallel/data/example/"
#%%
#dir_base = "/scratch/y08/mt3516/fff/fu_dam/300x20x8/"



inputs_to_file(prob=prob, grids=grids, solver=solver, dir=dir_base);
#%%
"""
Now execute the command in parallel from terminal/bash script.
eg run from MIDParallel/
>>mpiexec -n 2 julia -e 'using MIDParallel; par_spectrum_from_file(dir="/Users/matt/phd/MIDParallel/data/example/")'
>>mpiexec -n 2 julia -e 'using MIDParallel; par_spectrum_from_file(dir="/home/149/mt3516/island_damping/MIDParallel/data/example/")'

See convergence.sh for other examples of running in parallel.

"""


par_post_process(dir_base) #unfort have we have to do this!
#now we can read the data in. first the eigenvalues,
evals = evals_from_file(dir=dir_base);


continuum_plot(evals, n=-1);


ind = find_ind(evals, 0.30)


ϕft = efunc_from_file(dir=dir_base, ind=ind);
potential_plot(ϕft, grids);
