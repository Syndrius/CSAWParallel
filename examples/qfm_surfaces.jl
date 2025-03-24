
using MID
using MIDParallel




function qfm_benchmark_q(r::Float64)
    #a = 0.954545
    #b = 1.515151
    
    a = 1.93333
    b = 1.66666

    return a + b*r^2, 2 * b * r
end

#%%
#define the problem to solve

R0=4.0

#amp needs further thought!
#define the non-resonant island
isl = init_island(m0=3, n0=2, A=0.005)

geo = init_geo(R0=R0)

prob = init_problem(q=qfm_benchmark_q, geo=geo, isl=isl)


qlist, plist = farey_tree(5, 2, 1, 3, 1)

par_construct_surfaces(plist, qlist, 0.5 .* ones(length(qlist)), prob, "/Users/matt/phd/MIDParallel/data/qfm/")
