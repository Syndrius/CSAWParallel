#this example with the chaos system seems to be working ok.
#actual example of qfm case.
#note that this assumes the surfaces have already been generated, example can be found in qfm_surfaces.jl
#TODO
#example of this is worthwhile to see how the surfaces are passed around.


using MID
using MIDParallel
using MIDViz
using JLD2
using Plots; plotlyjs()
#%%
#generating qfm surfaces in parallel requires they are combined later.
#this assumes qfm_surfaces.jl has been run with 2 procs.
#gather_surfs("/Users/matt/phd/MIDParallel/data/qfm/", 2)
#gather_surfs("/home/149/mt3516/island_damping/MIDParallel/data/qfm/", 2)
#%%
R0=4.0

k = 0.0005
isl = init_island(m0=3, n0=-2, A=k/3)

geo = init_geo(R0=R0)

prob = init_problem(q=cantori_q, geo=geo, isl=isl,type=:flux)

#qlist, plist = farey_tree(3, 2, 1, 3, 1)
rationals = lowest_rationals(6, prob.q(0.0)[1], prob.q(1.0)[1])

guess_list = surface_guess(rationals, prob.q)
#%%
surfs = construct_surfaces(rationals, guess_list, prob);

plot_surfs(surfs);
#%%
#%%
save_object("/Users/matt/phd/MIDParallel/data/qfm/surfaces.jld2", surfs)
#%%
#surfs = load_object("/Users/matt/phd/MIDParallel/data/qfm/surfaces.jld2");
#surfs = load_object("/Users/matt/phd/MID/data/surfaces/total_bench_surfs.jld2");
#%%
Nr = 30
sgrid = init_grid(type=:rf, N=Nr, start=0.4, stop=0.9)
ϑgrid = init_grid(type=:as, N = 2, start = 1)
φgrid = init_grid(type=:as, N = 1, start = -1)
grids = init_grids(sgrid, ϑgrid, φgrid)
#%%
Nr = 20
sgrid = init_grid(type=:rf, N=Nr, start=0.4, stop=0.9)
ϑgrid = init_grid(type=:af, N = 4, pf=1)
φgrid = init_grid(type=:as, N = 1, start = -1)
grids = init_grids(sgrid, ϑgrid, φgrid)
#%%
Nr = 20
sgrid = init_grid(type=:rf, N=Nr, start=0.4, stop=0.9)
ϑgrid = init_grid(type=:af, N = 3, pf=1)
φgrid = init_grid(type=:af, N = 1, pf=-1)

grids = init_grids(sgrid, ϑgrid, φgrid)
#%%

solver = init_solver(nev=100, targets=[0.2, 0.3, 0.4], prob=prob)
#%%
dir_base = "/Users/matt/phd/MIDParallel/data/qfm/"
#dir_base = "/home/149/mt3516/island_damping/MIDParallel/data/qfm/"
dir_base = "/scratch/y08/mt3516/test/"
#%%
inputs_to_file(prob=prob, grids=grids, solver=solver, dir=dir_base)

#note the use of surfaces here is currently not very clear!!
"""
Now we would use the terminal to actually run the process
This kind of already assumes the surfaces have been created
>>mpiexec -n 2 julia -e 'using MIDParallel; qfm_spectrum_from_file(dir="/Users/matt/phd/MIDParallel/data/qfm/", qfm_surfs="/Users/matt/phd/MIDParallel/data/qfm/surfaces.jld2")'
>>mpiexec -n 2 julia -e 'using MIDParallel; using MID; qfm_spectrum_from_file(dir="/home/149/mt3516/island_damping/MIDParallel/data/qfm/", qfm_surfs="/home/149/mt3516/island_damping/MIDParallel/data/qfm/qfm_benchmark_surfaces.jld2", target_freq=0.3, nev=200)'
>>mpiexec -n 2 julia -e 'using MIDParallel; qfm_spectrum_from_file(dir="/scratch/y08/mt3516/test/", qfm_surfs="/scratch/y08/mt3516/qfm_data/surfaces/total_chaos_surfs.jld2")'
"""


#now we can actually have a look at this and see what is going on.
#first process the output
#%%
par_post_process(dir_base)
#%%

evals = evals_from_file(dir=dir_base);

continuum_plot(evals)#, ymax=40.5)

ind = find_ind(evals, 0.278)
ind = find_ind(evals, 0.281)
ind = find_ind(evals, 0.2818)

ϕft = efunc_from_file(dir = dir_base, ind=ind);

potential_plot(ϕft, grids)
