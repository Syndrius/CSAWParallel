

#using MPI
using MID
using MIDParallel

#this is now working, it is more fical than we would like, but it does work
#Easiest way to ensure this works (on mac) is to install MPI@0.19.2, which is the latest version compatible with slepcwrap.

#seems like this version will be very difficult to get working because Petsc seems to have a different default of MPI than julia's MPI. Proper installation may get this to work who knows.
#still not sure which is actually using MPICH, either way it is heking annoying.
#I think we can be reasnably confident that julia is using OpenMPI now, so not sure why Petsc cant work with MPICH yet seems to be running it?
#so now this seems to be working, but Slepc is hadning out garbage or segfaulting. so eps is clearly cooked.
#works perfectly with a single proc. Something cooked is going on.
#hopefully segfaults are just due to matrices being cooked for more procs, as evals seem to be getting closer to zero.

N = 100; 
#the collect is a bit annoying, but ok because we will typically use a clustered grid.
rgrid = collect(LinRange(0, 1, N));

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=Axel_q, geo=geo); #probbaly should use geo if it is part of prob,
#even if it is not really used.
grids = init_grids(rgrid=rgrid, mstart=2, mcount=2, nstart=-2, ncount=1);
tae_freq = (0.381 / geo.R0)^2

par_construct_and_solve(prob=prob, grids=grids, σ=tae_freq, dir="data/")

#test(0.1)