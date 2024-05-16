
using MID
using MIDParallel

#this is now working ish. we will have one more crack at getting petsc and shit to work.
#solving from file to avoid many issues.

N = 100
rgrid = collect(LinRange(0, 1, N));

#this is fkn annoying
geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo); 
grids = init_grids(rgrid=rgrid, mstart=2, mcount=2, nstart=-2, ncount=1);
dir = "data/serial/"
tae_freq = (0.381 / geo.R0)^2

par_solve_from_file(dir=dir)

#=

#solve_from_file(dir=dir, grids=grids, R0=prob.geo.R0, σ=tae_freq)

dir = "data/par/"
tae_freq = (0.381 / geo.R0)^2
solve_from_file(dir=dir, grids=grids, R0=prob.geo.R0, σ=tae_freq)



reconstruct_continuum(ω = ω, ϕ = ϕ, grids = grids)

tae_ind = find_ind(ω, 0.381)
plot_potential(r=rgrid, ϕ=ϕ, ind=tae_ind, pmd=grids.pmd, n=1)

=#
