"""

Slepc can be tricky to get working. This file, combined with solve.jl, show how we can still make use of MPI to construct the matrices but then return to Arpack to solve the matrices without MPI.

This file should be called with 
>>mpiexecjl -n (num_procs) julia examples/NoSlepc/construct.jl

Then the file solve.jl should be called in serial with
>>julia examples/NoSlepc/solve.jl

Solutions will be written to data/example/

"""


using MID 
using MIDParallel

#Define the problem and grid, alternativly can be read from file.
N = 100;
rgrid = collect(LinRange(0, 1, N));
geo = GeoParamsT(R0=10.0)
prob = init_problem(q=Axel_q, geo=geo); 
grids = init_grids(rgrid=rgrid, mstart=2, mcount=2, nstart=-2, ncount=1);


#construct the matrices in parallel and write them to file.
par_construct_to_file(prob=prob, grids=grids, dir="data/example/")