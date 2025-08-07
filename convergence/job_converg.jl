
using MID
using FFTW
using JLD2

include("/home/149/mt3516/island_damping/MID/Convergence/FEMConvergence/solution.jl") #going to be annoying af
include("/home/149/mt3516/island_damping/MIDParallel/convergence/fff.jl") #going to be annoying af

Nfl = [8, 10, 12]
#Nsl = [2, 3]
ml = [1, 2, 3]
nl = [1, 2, 3]
zl = [1, 2, 3]
#dir_base = "/Users/matt/phd/MIDParallel/data/convergence"
dir_base = "/scratch/y08/mt3516/Helmholtz/test"
dvals, dfuncs = fff_convergence(Nfl, ml, nl, zl, dir_base)
save_object(dir_base*"dvals.jld2", dvals)
save_object(dir_base*"dfuncs.jld2", dfuncs)
