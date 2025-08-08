
#testing the eigenvalue tolerances
#ideally we want to remove overlapping eigenvalues, and be more confident in the ones we get.
using MIDParallel
using MID
using Plots
using MIDViz
using Plots; plotlyjs()
#%%

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
solver = init_solver(nev=10, targets=[0.23, 0.29, 0.3, 0.31, 0.5], prob=prob)
#solver = init_solver(nev=13, target=0.31, prob=prob)



#%%

#looks like full path is needed... a bit annoying tbh.
#dir_base = "/scratch/y08/mt3516/test/"
dir_base = "/Users/matt/phd/MIDParallel/data/example/"

inputs_to_file(prob=prob, grids=grids, solver=solver, dir=dir_base)


#%%

par_post_process(dir_base)

evals = evals_from_file(dir=dir_base)
continuum_plot(evals)
length(evals.ω)


display(evals.ω[1:10])
display(evals.ω[11:20])

evals.ω[2]
evals.ω[12]

ϕft1 = efunc_from_file(dir=dir_base, ind=8)
potential_plot(ϕft1, grids)
ϕft2 = efunc_from_file(dir=dir_base, ind=18)
potential_plot(ϕft2, grids)
#%%
rg = MID.inst_grid(rgrid)
plot(rg, abs.(ϕft1[:, 2, 1]))
plot!(rg, abs.(ϕft2[:, 2, 1]))
