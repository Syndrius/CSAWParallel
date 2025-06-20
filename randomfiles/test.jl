

using MID
using MIDParallel
using MIDViz
using JLD2
using Plots; plotlyjs()
#%%
rgrid = init_grid(type=:rf, N = 20, start=0.05, stop=0.95, sep1=0.5, sep2=0.66, frac=0.5)
θgrid = init_grid(type=:af, N = 5, pf=3) 
ζgrid = init_grid(type=:af, N = 2, pf=-2)
grids = init_grids(rgrid, θgrid, ζgrid)
#%%

k1 = 0.00055
#k2 = 0.0003
#k3 = 0.0000
geo = init_geo(R0=4.0)
isl1 = init_island(m0=3, n0=-2, A=k1, flux=true)
isl2 = init_island(m0=5, n0=-3, A=k1, flux=true)
#isl3 = init_island(m0=8, n0=-5, A=k3/8) #unsure if we will want this one as well

isls = MID.Geometry.IslandT[isl1, isl2]

prob = init_problem(geo=geo, q=MID.Equilibrium.cyl_qfm_q, isls=isls, met=:cylinder, type=:flux)
#%%
solver = init_solver(prob=prob, target=0.0, nev=100)
#%%
dir_base = "/Users/matt/phd/MIDParallel/data/qfm/"

inputs_to_file(prob=prob, grids=grids, solver=solver, dir=dir_base)

#%%
par_post_process(dir_base)

evals = evals_from_file(dir=dir_base);

continuum_plot(evals)
