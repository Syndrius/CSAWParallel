
#real juicy one!


using MID 
using MIDParallel


#first we define the problem and write to file.
#this is identical to MID.
Nr=20;
Nθ=5;
Nζ=2;
#rgrid = collect(LinRange(0, 1, N));
geo = GeoParamsT(R0=10.0);
prob = init_problem(q=Axel_q, geo=geo); 
rgrid = init_fem_grid(N=Nr)
θgrid = init_fem_grid(N=Nθ, pf=2)
ζgrid = init_fem_grid(N=Nζ, pf=-2)
grids = init_grids(rgrid, θgrid, ζgrid);
#tae_freq = 0.396 #/ geo.R0)^2; #previously found tae_freq.

#looks like full path is needed... a bit annoying tbh.
dir_base = "/home/149/mt3516/island_damping/MIDParallel/data/example/"
#dir_base = "data/example/"

inputs_to_file(prob=prob, grids=grids, dir=dir_base);

"""
Now execute the command in parallel from terminal/bash script.
eg run from MIDParallel/
>>mpiexecjl -n 2 julia -e 'using MIDParallel; using MID; par_spectrum_from_file(dir="data/example/", freq=0.396)'
>>mpiexec -n 2 julia -e 'using MIDParallel; using MID; par_spectrum_from_file(dir="/home/149/mt3516/island_damping/MIDParallel/data/example/", freq=0.396)'

See convergence.sh for other examples of running in parallel.

"""

#now we can read the data in. first the eigenvalues,
ω = par_vals_from_file(dir_base*"vals.dat", geo.R0);

#then the efuncs
ϕ = par_funcs_from_file(dir_base*"funcs.dat", length(ω), grids);

#or we can read a single eigenfunction, which is useful for larger datasets, here we read the third eigenvalue.
ϕ3 = par_func_from_file(dir_base*"funcs.dat", 7, grids);

#now we can plot the TAE, which will be the first efuncs as we have specified the desired frequency to solve for.
display(ω[1]); #should find something that is kinda like a tae with ω=0.396
#plot_potential(ϕ, grids, 1, 1);


ϕms = mode_structure(ϕ, grids);
#why does this look way better than the serial case??? Different solver??
#pretty wild difference though!
#wot is going on here??? -> eval is same to ~10 places. How can efunc be so different?
#we are using the same function to read them in??
plot_potential(ϕms, grids, 1, 1);


#cool, this still work for a single func.
ϕ3ms = mode_structure(ϕ3, grids);
#the individual eigenfunction only contains a single index so we pass in ind=1
plot_potential(ϕ3ms, grids, 1, 2);


haskey(ENV, "SLEPC_DIR")

print(ENV)