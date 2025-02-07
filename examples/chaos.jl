
#testing if the chaos example will work in parallel. 
#may be a good comparison with fortran. We will also want to test different boundary conditions.
using MID
using MIDParallel

Nr = 100;

geo = GeoParamsT(R0=10.0)


isl = IslandT(m0=2.0, n0=-1.0, A=0.0e-5)
isl2 = IslandT(m0=3.0, n0=-2.0, A=0.0e-5)
#rgrid = collect(LinRange(0, 1, N));
prob = init_problem(q=chaos_q, geo=geo, isl=isl, isl2=isl2, met=cylindrical_metric!); 

prob = init_problem(q=Axel_q, geo=geo, met=cylindrical_metric!); 
#rgrid = rfem_grid(N=Nr)
#θgrid = asm_grid(start=2, N = 2)
#ζgrid = asm_grid(start=-2, N = 2)
#grids = init_grids(rgrid, θgrid, ζgrid);

Nr = 30
Nθ = 6
Nζ = 3
rgrid = rfem_grid(N=Nr, gp=4, sep1=0.4, sep2=0.6, frac=0.0);
θgrid = afem_grid(N=Nθ, pf=2, gp=4);
ζgrid = afem_grid(N=Nζ, pf=-2, gp=4);

#looks like full path is needed... a bit annoying tbh.
dir_base = "/home/149/mt3516/island_damping/MIDParallel/data/example/"
#dir_base = "data/example/"

inputs_to_file(prob=prob, grids=grids, dir=dir_base);


"""
Now execute the command in parallel from terminal/bash script.
eg run from MIDParallel/
>>mpiexecjl -n 2 julia -e 'using MIDParallel; using MID; par_spectrum_from_file(dir="data/example/", target_freq=0.29)'
>>mpiexec -n 2 julia -e 'using MIDParallel; using MID; par_spectrum_from_file(dir="/home/149/mt3516/island_damping/MIDParallel/data/example/", target_freq=0.30)'

See convergence.sh for other examples of running in parallel.

"""


process_hdf5(dir_base) #unfort have we have to do this!
#now we can read the data in. first the eigenvalues,
evals = evals_from_file(dir=dir_base);


continuum_plot(evals);

tae_ind = find_ind(evals, 0.268)


#the tae for this small res is ind=14.
ϕft = efunc_from_file(dir = dir_base, ind=tae_ind);
#potential plotting is fked af.
potential_plot(ϕft, grids);