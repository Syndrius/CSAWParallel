

geo = init_geometry(:tor, R0=4.0)

isl = init_island(m0=3, n0=1, A=0.1)

fields = init_fields(:ψ, q=cantori_q, isl=isl)

prob = init_problem(geometry=geo, fields=fields)

sgrid = init_grid(:s, 15, start=0.15, stop=0.9)
ϑgrid = init_grid(:ϑ, 4, pf=1)
ζgrid = init_grid(:ζ, 1, pf=-1)

grids = init_grids(sgrid, ϑgrid, ζgrid)


solver = init_solver(prob=prob, targets=[0.25, 0.4], nev=5)


dir = abspath(joinpath(pathof(CSAWParallel), "../../test/data/"))
inputs_to_file(dir=dir, prob=prob, grids=grids, solver=solver)

surf_dir = abspath(joinpath(pathof(CSAWParallel), "../../test/data/benchmark_surfaces.jld2"))
par_compute_spectrum(dir, surf_dir)

par_post_process(dir)

evals = evals_from_file(dir)

tae_ind = find_ind(evals, 0.25)
tae_freq = evals.ω[tae_ind]

@test tae_ind == 1
@test tae_freq ≈ 0.24188 atol=0.001

rm(joinpath(dir, "evals.jld2"))
rm(joinpath(dir, "grids.jld2"))
rm(joinpath(dir, "solver.jld2"))
rm(joinpath(dir, "prob.jld2"))
rm(joinpath(dir, "errs.jld2"))
rm(joinpath(dir, "unique_inds.jld2"))
rm(joinpath(dir, "vals_raw.jld2"))
rm(joinpath(dir, "efuncs_raw/"), recursive=true)
rm(joinpath(dir, "efuncs/"), recursive=true)
rm(joinpath(dir, "efuncs_ft/"), recursive=true)

