"""

Basic use case of MIDParallel. This file is designed for running in the repl, not in parallel.
Gives example of seting up the problem, then problme should be solved outside the repl in parallel, then this file can read the outputs.

Probably want to change this to a different q-profile eventually.

"""

using MID 
using MIDParallel


#first we define the problem and write to file.
#this is identical to MID.
N=100;
rgrid = collect(LinRange(0, 1, N));
geo = GeoParamsT(R0=10.0);
prob = init_problem(q=Axel_q, geo=geo); 
grids = init_grids(rgrid=rgrid, mstart=2, mcount=2, nstart=-2, ncount=1);
tae_freq = (0.381 / geo.R0)^2; #previously found tae_freq.

inputs_to_file(prob=prob, grids=grids, dir="data/example/");

"""
Now execute the command in parallel from terminal/bash script.
eg run from MIDParallel/
>>mpiexecjl -n (num_procs) julia -e 'using MIDParallel; Using MID; par_spectrum_from_file(dir="data/example/", freq=0.00145)'

See convergence.sh for other examples of running in parallel.

"""

#now we can read the data in. first the eigenvalues,
ω = par_vals_from_file("data/example/vals.dat", geo.R0);

#then the efuncs
ϕ = par_funcs_from_file("data/example/funcs.dat", length(ω), N, grids.pmd.count, grids.tmd.count);

#now we can plot the TAE, which will be the first efuncs as we have specified the desired frequency to solve for.
display(ω[1]); #should find a tae with normalised frequency 0.3812
plot_potential(r = rgrid, ϕ=ϕ, ind=1, pmd = grids.pmd, n=1);