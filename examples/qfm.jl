"""
Example using Quadratic Flux Minimising (QFM) coordinates.
"""

using ChaoticShearAlfvenWaves
using CSAWParallel
using CSAWViz

# QFM surfaces are generated in the companion package CSAWCantori.
# We have some precomputed used for testing.
# In parallel we read from file
surf_dir = abspath(joinpath(pathof(CSAWParallel), "../../test/data/benchmark_surfaces.jld2"))

# This case uses a very large non-resonant perturbation
isl = init_island(m0=3, n0=1, A=0.1)
geo = init_geometry()
fields = init_fields(:ψ, q=cantori_q, isl=isl)
prob = init_problem(geometry=geo, fields=fields)

# We use finite elements for mapping below.
# The unrealistic perturbation causes issues at the boundaries so we shrink the size.
sgrid = init_grid(:s, 30, start=0.15, stop=0.9)
ϑgrid = init_grid(:ϑ, 5, pf=1)
ζgrid = init_grid(:ζ, 1, pf=-1)

grids = init_grids(sgrid, ϑgrid, ζgrid)

#in parallel with larger matrix sizes, we cannot solve the full spectrum.
#we can use shift and invert to target specific eigenvalues
#or if we want a larger part of the spectrum, we can slice with multiple shift and inverts
#this initiates the solver to target 0.2, 0.3 and 0.4, with 50 eigenvalues each
solver = init_solver(prob=prob, targets=[0.2, 0.3, 0.4], nev=50)

#this writes the input into the testing directory
dir = abspath(joinpath(pathof(CSAWParallel), "../../test/data/"))
inputs_to_file(dir=dir, prob=prob, grids=grids, solver=solver)

#compute the spectrum
#also passing in the location of the surfs
#see examples/run.sh to run this in parallel.
par_compute_spectrum(dir, surf_dir)

#the current implementation of PetscWrap.jl does not convert Petsc arrays to julia arrays.
#this means we have to post process the solutions in serial after solving
par_post_process(dir)

#load the resuls
evals = evals_from_file(dir)

#despite the perturbation the continuum appears normal
continuum_plot(evals)

#we look at specific solutions
#a TAE and a typical continuum solution
tae_ind = find_ind(evals, 0.25) 
cont_ind = find_ind(evals, 0.20)

#load the specific solutions
ϕft_tae = efunc_from_file(dir, tae_ind);
ϕft_cont = efunc_from_file(dir, cont_ind);
#and plot
harmonic_plot(ϕft_tae, grids) 
harmonic_plot(ϕft_cont, grids)

"""
We can also Map the solutions back to toroidal coordinates.
With CSAWParallel output this is done straight from file.
"""

# Note that this mapping can be slow for large grids.
# The radial grid is reduced again to prevent issues with interpolation at the edges.
ψgrid = init_grid(:ψ, 80, start=0.25, stop=0.8)
θgrid = init_grid(:θ, 30)
φgrid = init_grid(:φ, 10)
tor_grids = init_grids(ψgrid, θgrid, φgrid)

# Map the solutions
qfm_spectrum_to_tor(dir, tor_grids, surf_dir);

#results are written to a subfolder
map_dir = joinpath(dir, "tor_map/")

#we can load the mapped solutions as normal
ϕ_cont_tor = efunc_from_file(map_dir, cont_ind, ft=false);

#and plot
contour_plot(ϕ_cont_tor, tor_grids)
