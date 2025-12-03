"""
Example showing basic usage.
Usage is almost identical to CSAW, however, inputs must be read in from file.
"""

using ChaoticShearAlfvenWaves
using CSAWParallel
using CSAWViz

#inputs are generated with CSAW
#with all of the previous options 
geo = init_geometry()
fields = init_fields()
prob = init_problem(geometry=geo, fields=fields)
ψgrid = init_grid(:ψ, 100)
θgrid = init_grid(:sm, 2, start=1)
φgrid = init_grid(:sm, 1, start=-1) 
grids = init_grids(ψgrid, θgrid, φgrid)

#in parallel with larger matrix sizes, we cannot solve the full spectrum.
#we can use shift and invert to target specific eigenvalues
#or if we want a larger part of the spectrum, we can slice with multiple shift and inverts
#this initiates the solver to target 0.2, 0.3 and 0.4, with 50 eigenvalues each
solver = init_solver(prob=prob, targets=[0.2, 0.3, 0.4], nev=50)


#this writes the input into the testing directory
dir = abspath(joinpath(pathof(CSAWParallel), "../../test/data/"))
inputs_to_file(dir=dir, prob=prob, grids=grids, solver=solver)

#we can solve directly in the REPL
#to actually solve in parallel see examples/run.sh
par_compute_spectrum(dir)

#the current implementation of PetscWrap.jl does not convert Petsc arrays to julia arrays.
#this means we have to post process the solutions in serial after solving
par_post_process(dir)

#we can then read the solutions from file
evals = evals_from_file(dir)

#and visualise in the same way
continuum_plot(evals)

#looking at specific TAE solution
ind = find_ind(evals, 0.3)

#read the specific solution
ϕft = efunc_from_file(dir, ind);

#and plot
harmonic_plot(ϕft, grids)
