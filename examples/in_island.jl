#TODO
#probably requires that MID is fixed first
using MID
using MIDParallel
#%%

geo = init_geo(R0=10.0)
isl = init_island(m0=2, n0=-1, w=0.03, qp=2.0, r0=0.5)
prob = init_problem(q=MID.Equilibrium.island_coords_q, met=:island, geo=geo, isl=isl)
#%%

Nκ = 20
Nᾱ = 6
Nτ = 2

κgrid = init_grid(type=:rf, N=Nκ, stop=0.999, left_bc=false)
ᾱgrid = init_grid(type=:af, N=Nᾱ, pf=0)
τgrid = init_grid(type=:af, N=Nτ, pf=0)

grids = init_grids(κgrid, ᾱgrid, τgrid)
#%%
solver = init_solver(nev=100, targets=[0.0, 0.33, 0.8], prob=prob)
#%%
dir_base = "/scratch/y08/mt3516/test/"
#%%
inputs_to_file(prob=prob, grids=grids, solver=solver, dir=dir_base);
#%%
"""
Now execute the command in parallel from terminal/bash script.
eg run from MIDParallel/
>>mpiexec -n 2 julia -e 'using MIDParallel; par_spectrum_from_file(dir="/Users/matt/phd/MIDParallel/data/example/")'
>>mpiexec -n 2 julia -e 'using MIDParallel; par_spectrum_from_file(dir="/scratch/y08/mt3516/test/")'

See convergence.sh for other examples of running in parallel.

"""
#%%
