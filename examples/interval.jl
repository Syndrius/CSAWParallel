
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
solver = init_solver(prob=prob, left=0.15, right=0.35)

0.2^2 / geo.R0^2
0.3^2 / geo.R0^2
0.35^2 / geo.R0^2
#%%

#looks like full path is needed... a bit annoying tbh.
#dir_base = "/scratch/y08/mt3516/test/"
dir_base = "/Users/matt/phd/MIDParallel/data/example/"
#%%

inputs_to_file(prob=prob, grids=grids, solver=solver, dir=dir_base);

#%%
par_post_process(dir_base)

evals = evals_from_file(dir=dir_base)

continuum_plot(evals)
