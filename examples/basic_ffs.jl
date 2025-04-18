
using MID 
using MIDParallel
using MIDViz
#%%

#first we define the problem and write to file.
#this is identical to MID.
Nr=40;
Nθ=6;
#rgrid = collect(LinRange(0, 1, N));
geo = init_geo(R0=4.0);
prob = init_problem(q=fu_dam_q, geo=geo); 
rgrid = MID.Structures.rfem_grid(N=Nr);
θgrid = MID.Structures.afem_grid(N=Nθ, pf=1);
ζgrid = MID.Structures.asm_grid(start=-1, N = 1);
grids = init_grids(rgrid, θgrid, ζgrid);
#%%
#tae_freq = 0.396 #/ geo.R0)^2; #previously found tae_freq.


solver = init_solver(nev=100, targets=[0.0, 0.33, 0.8], prob=prob)
#%%
#looks like full path is needed... a bit annoying tbh.
#dir_base = "/home/149/mt3516/island_damping/MIDParallel/data/example/"
dir_base = "/Users/matt/phd/MIDParallel/data/example/"

#%%
inputs_to_file(prob=prob, grids=grids, solver=solver, dir=dir_base);
#%%
"""
Now execute the command in parallel from terminal/bash script.
eg run from MIDParallel/
>>mpiexec -n 2 julia -e 'using MIDParallel; par_spectrum_from_file(dir="/Users/matt/phd/MIDParallel/data/example/")'

"""


par_post_process(dir_base) #unfort have we have to do this!
#now we can read the data in. first the eigenvalues,
evals = evals_from_file(dir=dir_base);


continuum_plot(evals);


ind = find_ind(evals, 0.3364)


ϕft = efunc_from_file(dir=dir_base, ind=ind);
potential_plot(ϕft, grids);
