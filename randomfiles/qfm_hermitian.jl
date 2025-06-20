
#looks like the problem is that our matrices are not exactly positive definite, which seems to become a srs problemo for qfm
#see if we can replicate it in a smaller sample
using MID
using MIDParallel

#%%
rgrid = init_grid(type=:rf, N = 3, gp=4)
θgrid = init_grid(type=:af, N = 3, gp=4, pf=0)
ζgrid = init_grid(type=:af, N = 3, gp=4, pf=0)
#θgrid = init_grid(type=:as, N = 2, start=1)
#ζgrid = init_grid(type=:as, N = 1, start=-1)

grids = init_grids(rgrid, θgrid, ζgrid)
#%%
#surfs = load_object("/Users/matt/phd/MID/data/surfaces/island_21a.jld2");
#%%
geo = init_geo(R0=4.0)

isl21a = init_island(m0=2, n0=-1, w=0.05, r0=0.5, qp=2.0)
#isl21a = init_island(m0=2, n0=-1, A=0.001)
#start with no islands
prob = init_problem(geo=geo, q=MID.Equilibrium.island_equiv_q, met=:cylinder, isl=isl21a)

#%%
solver = init_solver(nev=10, targets=[0.2, 0.3, 0.4], prob=prob)
#%%
dir_base = "/Users/matt/phd/MIDParallel/data/qfm/"
#%%
inputs_to_file(prob=prob, grids=grids, solver=solver, dir=dir_base)
