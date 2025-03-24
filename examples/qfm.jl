
#actual example of qfm case.


using MID
using MIDParallel
#%%
Nr = 100

geo = init_geo(R0=4.0)

isl = init_island(m0=3, n0=2, A=0.005)

prob = init_problem(q=qfm_benchmark_q, geo=geo, isl=isl)

sgrid = init_grid(type=:rf, N=Nr, start=0.4, stop=0.7)
ϑgrid = init_grid(type=:as, N = 4, start = 2)
φgrid = init_grid(type=:as, N = 3, start = -2)

grids = init_grids(sgrid, ϑgrid, φgrid)

dir_base = "/Users/matt/phd/MIDParallel/data/qfm/"
#%%
inputs_to_file(prob=prob, grids=grids, dir=dir_base)

"""
Now we would use the terminal to actually run the process
This kind of already assumes the surfaces have been created
>>mpiexecjl -n 2 julia -e 'using MIDParallel; using MID; qfm_spectrum_from_file(dir="/Users/matt/phd/MIDParallel/data/qfm/", qfm_surfs="/Users/matt/phd/MIDParallel/data/qfm/qfm_benchmark_surfaces.jld2", target_freq=0.3, nev=200)'
"""

