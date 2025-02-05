
#real juicy one!


using MID 
using MIDParallel

using Plots; plotlyjs()

#first we define the problem and write to file.
#this is identical to MID.
Nr=5;
Nθ=2;
Nζ=1;
#rgrid = collect(LinRange(0, 1, N));
geo = GeoParamsT(R0=10.0);
prob = init_problem(q=Axel_q, geo=geo); 

rgrid = rfem_grid(N=Nr, gp=5)
θgrid = afem_grid(N=Nθ, pf=0, gp=5);
ζgrid = afem_grid(N=Nζ, pf=-0, gp=5);
ae49d9024e63c680edb3c6ec3
grids = init_grids(rgrid, θgrid, ζgrid);
#tae_freq = 0.396 #/ geo.R0)^2; #previously found tae_freq.

#looks like full path is needed... a bit annoying tbh.
#dir_base = "/home/149/mt3516/island_damping/MIDParallel/data/example/"
dir_base = "data/example/"

#dir_base = "/scratch/y08/mt3516/fff/fu_dam/300x20x8/"

f = zeros(2, 3, 4)

display(size(f)...)

mkpath(dir_base)


inputs_to_file(prob=prob, grids=grids, dir=dir_base);

"""
Now execute the command in parallel from terminal/bash script.
eg run from MIDParallel/
>>mpiexecjl -n 2 julia -e 'using MIDParallel; using MID; par_spectrum_from_file(dir="data/example/", target_freq=0.31)'
>>mpiexec -n 2 julia -e 'using MIDParallel; using MID; par_spectrum_from_file(dir="/home/149/mt3516/island_damping/MIDParallel/data/example/", target_freq=0.301)'

See convergence.sh for other examples of running in parallel.

"""

process_hdf5(dir_base)

#now we can read the data in. first the eigenvalues,
evals = evals_from_file(dir=dir_base);

#wot is this, this is cooked af
continuum_plot(evals, ymax=1);

display(evals.ω[1:5])

#found the tae, freq is significantly different for this low res example.
tae_ind = find_ind(evals, 0.3835)

display(length(evals.ω))

tae_ind = find_ind(evals, 0.384)


#the tae for this small res is ind=14.
ϕft = efunc_from_file(dir = dir_base, ind=3);
#potential plotting is fked af.
potential_plot(ϕft, grids);

ϕ = efunc_from_file(dir = dir_base, ind=tae_ind, ft=false);
#potential plotting is fked af.
potential_plot(ϕ, grids);

#looks good, other than mode labels.
contour_plot(ϕ, grids)

surface_plot(ϕ, grids)


prob, grids = inputs_from_file(dir=dir_base);

display(grids)


grid_points = Tuple{Int, Int, Int}[]

push!(grid_points, (1, 2, 3))

display(size(grid_points)[1])

rgrid, θgrid, ζgrid = inst_grids(grids)

#adds back the periodicity.
z = zeros(ComplexF64, grids.r.N, grids.ζ.N+1)

z[:, 1:end-1] = ϕ[ :, 2, :]
z[:, end] = ϕ[ :, 2, 1]

ζgrid = range(0, 2π, grids.ζ.N+1)
p = contourf(ζgrid, rgrid, real.(z), levels=100, color=:turbo)

display(p)
