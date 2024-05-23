

#slepc only works with 1 core atm
#test if it is faster than arpack.
#this is significantly faster than slepc with only 1 core.

using MID

#this is now working ish. we will have one more crack at getting petsc and shit to work.
#solving from file to avoid many issues.

N = 1000;
#the collect is a bit annoying, but ok because we will typically use a clustered grid.
#rgrid = collect(LinRange(0, 1, N));
rgrid = clustered_grid(N, 0.4, 0.6, 0.5)

geo = GeoParamsT(R0=10.0)
isl = IslandT(A=4e-5, m0=5, n0=4)
prob = init_problem(q=Axel_q, geo=geo, isl=isl, δ=-4e-9); 
grids = init_grids(rgrid=rgrid, mstart=-3, mcount=12, nstart=-6, ncount=3, nincr=4);

dir = "data/"
tae_freq = (0.381 / geo.R0)^2



solve_from_file(dir=dir, grids=grids, R0=prob.geo.R0, σ=tae_freq)


