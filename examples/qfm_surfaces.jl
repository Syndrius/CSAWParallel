
using MID
using MIDParallel




#= This is defined inside MID now.
function qfm_benchmark_q(r::Float64)
    #a = 0.954545
    #b = 1.515151
    
    a = 1.93333
    b = 1.66666

    return a + b*r^2, 2 * b * r
end
=#

#%%
#define the problem to solve

R0=4.0

#amp needs further thought!
#define the non-resonant island
isl = init_island(m0=3, n0=2, A=0.005)

geo = init_geo(R0=R0)

prob = init_problem(q=qfm_benchmark_q, geo=geo, isl=isl)

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



#par_construct_surfaces(plist, qlist, 0.5 .* ones(length(qlist)), prob, "/home/149/mt3516/island_damping/MIDParallel/data/qfm/")
par_construct_surfaces(plist, qlist, 0.5 .* ones(length(qlist)), prob, "/Users/matt/phd/MIDParallel/data/qfm/")
