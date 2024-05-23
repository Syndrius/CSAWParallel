
#julia file to be called from terminal, creates the different inputs based on varying amplitude
using MID;


qprof = island_damping_q;
R0 = 10.0;
isl = IslandT(A=parse(Float64, ARGS[1]), m0=5, n0=4)
geo = GeoParamsT(R0=R0);
prob = init_problem(q=qprof, geo=geo, isl=isl, δ=-4e-9);
N = 1000;
rgrid = clustered_grid(N, 0.4, 0.6, 0.5)
grids = init_grids(rgrid=rgrid, mstart=-3, mcount=12, nstart=-6, ncount=3, nincr=4);
inputs_to_file(prob=prob, grids=grids, dir=ARGS[2]);