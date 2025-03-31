#this example with the chaos system seems to be working ok.
#actual example of qfm case.
#note that this assumes the surfaces have already been generated, example can be found in qfm_surfaces.jl


using MID
using MIDParallel
using MIDViz
using Plots; plotlyjs()
#%%
#generating qfm surfaces in parallel requires they are combined later.
#this assumes qfm_surfaces.jl has been run with 2 procs.
gather_surfs("/Users/matt/phd/MIDParallel/data/qfm/", 2)
#%%
Nr = 30


R0=10.0

#amp needs further thought!
#define the non-resonant island
k = 0.00022
isl = init_island(m0=5, n0=-2, A=k/5)
isl2 = init_island(m0=7, n0=-3, A=k/7)

geo = init_geo(R0=R0)

#to solve non-Hermitian
#flr = MID.Structures.FLRT(δ = 1e-18)
prob = init_problem(q=qfm_benchmark_q, met=:cylinder, geo=geo, isl=isl, isl2=isl2)#, flr=flr)

sgrid = init_grid(type=:rf, N=Nr, start=0.4, stop=0.7)
#ϑgrid = init_grid(type=:as, N = 4, start = 2)
#φgrid = init_grid(type=:as, N = 3, start = -2)
ϑgrid = init_grid(type=:af, N = 6, pf=2)
φgrid = init_grid(type=:af, N = 2, pf=-2)

grids = init_grids(sgrid, ϑgrid, φgrid)

dir_base = "/Users/matt/phd/MIDParallel/data/qfm/"
#dir_base = "/home/149/mt3516/island_damping/MIDParallel/data/qfm/"
#%%
inputs_to_file(prob=prob, grids=grids, dir=dir_base)

#note the use of surfaces here is currently not very clear!!
"""
Now we would use the terminal to actually run the process
This kind of already assumes the surfaces have been created
>>mpiexecjl -n 2 julia -e 'using MIDParallel; using MID; qfm_spectrum_from_file(dir="/Users/matt/phd/MIDParallel/data/qfm/", qfm_surfs="/Users/matt/phd/MIDParallel/data/qfm/qfm_benchmark_surfaces.jld2", target_freq=0.3, nev=200)'
>>mpiexec -n 2 julia -e 'using MIDParallel; using MID; qfm_spectrum_from_file(dir="/home/149/mt3516/island_damping/MIDParallel/data/qfm/", qfm_surfs="/home/149/mt3516/island_damping/MIDParallel/data/qfm/qfm_benchmark_surfaces.jld2", target_freq=0.3, nev=200)'
"""


#now we can actually have a look at this and see what is going on.
#first process the output
#%%
process_hdf5(dir_base)
#%%

evals = evals_from_file(dir=dir_base);

continuum_plot(evals)#, ymax=40.5)

length(evals.ω)
