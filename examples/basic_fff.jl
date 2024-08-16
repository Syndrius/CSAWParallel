
#real juicy one!


using MID 
using MIDParallel

using Plots; plotlyjs()


#first we define the problem and write to file.
#this is identical to MID.
Nr=50;
Nθ=5;
Nζ=2;
#rgrid = collect(LinRange(0, 1, N));
geo = GeoParamsT(R0=4.0);
prob = init_problem(q=fu_dam_q, geo=geo); 
rgrid = init_fem_grid(N=Nr)
θgrid = init_fem_grid(N=Nθ, pf=1);
ζgrid = init_fem_grid(N=Nζ, pf=-1);
grids = init_grids(rgrid, θgrid, ζgrid);
#tae_freq = 0.396 #/ geo.R0)^2; #previously found tae_freq.

#looks like full path is needed... a bit annoying tbh.
dir_base = "/home/149/mt3516/island_damping/MIDParallel/data/example/"
#dir_base = "data/example/"

inputs_to_file(prob=prob, grids=grids, dir=dir_base);

"""
Now execute the command in parallel from terminal/bash script.
eg run from MIDParallel/
>>mpiexecjl -n 2 julia -e 'using MIDParallel; using MID; par_spectrum_from_file(dir="data/example/", target_freq=0.31)'
>>mpiexec -n 2 julia -e 'using MIDParallel; using MID; par_spectrum_from_file(dir="/home/149/mt3516/island_damping/MIDParallel/data/example/", target_freq=0.396)'

See convergence.sh for other examples of running in parallel.

"""

#now we can read the data in. first the eigenvalues,
evals = evals_from_file(dir=dir_base);

plot_continuum(evals);

#found the tae, freq is significantly different for this low res example.
tae_ind = find_ind(evals, 0.3267)

display(evals.ω[tae_ind])



ϕft = efunc_from_file(dir = dir_base, ind=tae_ind);
#potential plotting is fked af.
plot_potential(ϕft, grids, 1);

#or choosing a specfic n to focus on.
plot_potential(ϕft, grids, 1, n=-1)

