

using MID
using MIDParallel
using Plots; plotlyjs()


Nr=500;
#rgrid = collect(LinRange(0, 1, N));
geo = GeoParamsT(R0=10.0);
prob = init_problem(q=island_damping_q, geo=geo); 
rgrid = init_fem_grid(N=Nr)
θgrid = init_sm_grid(start=2, count=2)
ζgrid = init_sm_grid(start=-2, count = 1)
grids = init_grids(rgrid, θgrid, ζgrid);
#tae_freq = 0.396 #/ geo.R0)^2; #previously found tae_freq.

#looks like full path is needed... a bit annoying tbh.
#dir_base = "/home/149/mt3516/island_damping/MIDParallel/data/example/"
#dir_base = "/scratch/y08/mt3516/fss/no_island/500x2x1/"
#mkpath(dir_base)

dir_base = "data/example/"

inputs_to_file(prob=prob, grids=grids, dir=dir_base);

evals = evals_from_file(dir=dir_base);

#wot is this, this is cooked af
plot_continuum(evals);


ϕft = efunc_from_file(dir = dir_base, ind=1);
#potential plotting is fked af.
plot_potential(ϕft, grids, 1);
