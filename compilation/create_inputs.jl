#creates a small test set for the warmup case to run.
using MID
using JLD2 #not ideal, save surfs should be in MID>

dir_base = ARGS[1]

Nr=4;
Nθ=2;
Nζ=2;
geo = init_geo(R0=4.0);
rgrid = init_grid(type=:rf, N=Nr, start=0.4, stop=0.6)
θgrid = init_grid(type=:af, N=Nθ, pf=1)
ζgrid = init_grid(type=:af, N=Nζ, pf=-1)

grids = init_grids(rgrid, θgrid, ζgrid);
isl = init_island(m0=3, n0=-2, A=0.01/3, flux=true)

prob = init_problem(q=cantori_q, geo=geo, type=:flux, isls=MID.IslandT[isl]); 

solver = init_solver(nev=2, targets=[0.2, 0.5], prob=prob)
dir = joinpath(dir_base, "data/")

inputs_to_file(prob=prob, grids=grids, solver=solver, dir=dir);

#now testing a qfm case!
dir = joinpath(dir_base, "qfm_data/")
surf_dir = joinpath(dir_base, "qfm_data/warmup_surfs.jld2")

surfs = construct_surfaces([(4, 3), (3, 2), (7, 5), (2, 1), (5, 4)], [0.5, 0.6, 0.58, 0.5, 0.5], prob, M=4, N=2)

save_object(surf_dir, surfs)
inputs_to_file(prob=prob, grids=grids, solver=solver, dir=dir);

