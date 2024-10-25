"""

Basic use case of MIDParallel. This file is designed for running in the repl, not in parallel.
Gives example of seting up the problem, then problme should be solved outside the repl in parallel, then this file can read the outputs.

Probably want to change this to a different q-profile eventually. 

"""

using MID 
using MIDParallel
using Plots; plotlyjs()


#first we define the problem and write to file.
#this is identical to MID.
Nr=100;
#rgrid = collect(LinRange(0, 1, N));
geo = GeoParamsT(R0=10.0);
prob = init_problem(q=Axel_q, geo=geo); 
rgrid = rfem_grid(N=Nr)
θgrid = asm_grid(start=2, N = 2)
ζgrid = asm_grid(start=-2, N = 1)
grids = init_grids(rgrid, θgrid, ζgrid);


#looks like full path is needed... a bit annoying tbh.
dir_base = "/home/149/mt3516/island_damping/MIDParallel/data/example/"
#dir_base = "data/example/"

inputs_to_file(prob=prob, grids=grids, dir=dir_base);

"""
Now execute the command in parallel from terminal/bash script.
eg run from MIDParallel/
>>mpiexecjl -n 2 julia -e 'using MIDParallel; using MID; par_spectrum_from_file(dir="data/example/", target_freq=0.29)'
>>mpiexec -n 2 julia -e 'using MIDParallel; using MID; par_spectrum_from_file(dir="/home/149/mt3516/island_damping/MIDParallel/data/example/", target_freq=0.38)'

See convergence.sh for other examples of running in parallel.

"""

#first 10 before being normed.
"""
vpr = 0.0014167655630368723 - 5.816043703656595e-20im
vpr = 0.0013798042087976256 - 5.303827685943013e-21im
vpr = 0.001376898964057661 - 1.3085633161023262e-19im
vpr = 0.001373290765449895 - 5.785241153220936e-20im
vpr = 0.0013669913435559468 + 1.3537179963861004e-20im
vpr = 0.0013599287375589363 - 6.325193351626356e-18im
vpr = 0.001352397218309042 + 4.7560503917302217e-17im
vpr = 0.001339380207669367 + 6.173880765132081e-16im
vpr = 0.0013344694061573198 + 1.1221046793288665e-16im
vpr = 0.0013156435132545142 + 1.0892707095796296e-16im
vpr = 0.0013124374116757099 + 1.62357826663483e-15im
"""


process_hdf5(dir_base) #unfort have we have to do this!
#now we can read the data in. first the eigenvalues,
evals = evals_from_file(dir=dir_base);


continuum_plot(evals);

println(evals.ω[1:20])

ind = find_ind(evals, 0.3764)


ϕft = efunc_from_file(dir=dir_base, ind=ind);
potential_plot(ϕft, grids);
ϕ = efunc_from_file(dir=dir_base, ind=1, ft=false);
#plot_potential(ϕ, grids);

#inconsistent with phase of two modes, here we have the og case for some reason...
contour_plot(ϕ, grids, ind=1)
surface_plot(ϕ, grids, ind=1)