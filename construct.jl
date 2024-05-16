

#here we consider the case where we just construct straight to file.
using MID #not sure if this is needed but seems un-ideal.
using MIDParallel



N = 100;
#the collect is a bit annoying, but ok because we will typically use a clustered grid.
rgrid = collect(LinRange(0, 1, N));

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo); 
grids = init_grids(rgrid=rgrid, mstart=2, mcount=2, nstart=-2, ncount=1);


par_construct_to_file(prob=prob, grids=grids, dir="data/par/")