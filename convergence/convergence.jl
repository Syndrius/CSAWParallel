#other than MIDParallel being broken, this seems to be working ok
#we will need to actually solve in parallel with the bash script, and post process
#but otherwise should be good to do some serious testing!

#note that this only works on the Helmholtz branch of MID.

using MID
using MIDParallel
#include("/Users/matt/phd/MID/convergence/FEMConvergence/solution.jl") #going to be annoying af
include("/home/149/mt3516/island_damping/MID/Convergence/FEMConvergence/solution.jl") #going to be annoying af
using Plots
using FFTW
#include("/Users/matt/phd/MIDParallel/convergence/fss.jl") #going to be annoying af
#include("/Users/matt/phd/MIDParallel/convergence/ffs.jl") #going to be annoying af
#include("/Users/matt/phd/MIDParallel/convergence/fff.jl") #going to be annoying af
include("/home/149/mt3516/island_damping/MIDParallel/convergence/fss.jl") #going to be annoying af
include("/home/149/mt3516/island_damping/MIDParallel/convergence/ffs.jl") #going to be annoying af
include("/home/149/mt3516/island_damping/MIDParallel/convergence/fff.jl") #going to be annoying af
#%%
#this case gives errors with 48 cores, may signify a srs problemo.
#this was with n+1 for y and n-1 for z. also with clustered grid, sep1=0.4, sep2=0.6, frac=0.5.
#this is a big cause for concern, unsure why this is happening.
#Nfl = [8, 10, 12]
#Nsl = [2, 3]
#ml = [1, 2, 3]
#nl = [1, 2, 3]
#zl = [1, 2, 3]

Nfl = [6, 8]
#Nsl = [2, 3]
ml = [1, 2]
nl = [1, 2]
zl = [1, 2]


#dir_base = "/Users/matt/phd/MIDParallel/data/convergence"
dir_base = "/scratch/y08/mt3516/Helmholtz/test"
#fss_convergence_inputs([10, 15], [2, 3], [1, 2], [1, 2], [1, 2], dir_base)
#ffs_convergence_inputs(Nfl, Nsl, ml, nl, zl, dir_base)
fff_convergence_inputs(Nfl, ml, nl, zl, dir_base)
#%%
#assuming this is fine, the rest is a ok.
#MIDParallel is cooked, think the wrong MPI version is loaded,
#so we will just do this in serial as if it is in parallel for now!
#imagine this cell being run by run.sh and done in parallel...
#this also assumes post-processing has been completed!
for i in 1:2
    dir = dir_base * "/fss"*string(i) * "/"
    #dir = dir_base * "/ffs"*string(i) * "/"
    #dir = dir_base * "/fff"*string(i) * "/"
    MID.Solve.spectrum_from_file(dir, true)
end
#%%
#dvals, dfuncs = fss_convergence([10, 15], [2, 3], [1, 2], [1, 2], [1, 2], dir_base)
#dvals, dfuncs = ffs_convergence(Nfl, Nsl, ml, nl, zl, dir_base)
dvals, dfuncs = fff_convergence(Nfl, ml, nl, zl, dir_base)

#ϕ = efunc_from_file(dir=dir_base*"/fff1/", ind=10)
#%%
#ok so this does work, but the lack of convergence for all eigenvalues is causing issues with comparison
#to the analytical case.
p = plot()
for i in 1:2
    scatter!(dvals[i, :])
end
display(p)
#%%
p = plot()
for i in 1:2
    scatter!(dfuncs[i, :])
end
display(p)
