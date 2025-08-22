
using MID
using MIDParallel

#this is probably a bit useless now!

#%%
#define the problem to solve

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


#%%

#this is not being split properly!
qlist, plist = farey_tree(4, 2, 1, 3, 1)

#%%

#par_construct_surfaces(plist, qlist, 0.5 .* ones(length(qlist)), prob, "/Users/matt/phd/MIDParallel/data/qfm/")

"""
Now run 
>>mpiexec -n 2 julia examples/qfm_surfaces.jl

Then the data must be combined in serial, by running
>>julia -e 'using MIDParallel; gather_surfs("/home/149/mt3516/island_damping/MIDParallel/data/qfm/", NPROCS)'
>>julia -e 'using MIDParallel; gather_surfs("/Users/matt/phd/MIDParallel/data/qfm/", 2)'
"""



par_construct_surfaces(plist, qlist, 0.5 .* ones(length(qlist)), prob, "/home/149/mt3516/island_damping/MIDParallel/data/qfm/")
#par_construct_surfaces(plist, qlist, 0.5 .* ones(length(qlist)), prob, "/Users/matt/phd/MIDParallel/data/qfm/")
