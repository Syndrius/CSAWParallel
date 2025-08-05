#other than MIDParallel being broken, this seems to be working ok
#we will need to actually solve in parallel with the bash script, and post process
#but otherwise should be good to do some serious testing!
#note that this only works on the test branch of MID.

using MID
using MIDParallel
include("/Users/matt/phd/MID/convergence/FEMConvergence/solution.jl") #going to be annoying af
using Plots
using FFTW
include("/Users/matt/phd/MIDParallel/convergence/fss.jl") #going to be annoying af
include("/Users/matt/phd/MIDParallel/convergence/ffs.jl") #going to be annoying af
include("/Users/matt/phd/MIDParallel/convergence/fff.jl") #going to be annoying af
#%%
Nfl = [6, 8]
Nsl = [2, 3]
ml = [1, 2]
nl = [1, 2]
zl = [1, 2]
dir_base = "/Users/matt/phd/MIDParallel/data/convergence"
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
    #dir = dir_base * "/fss"*string(i) * "/"
    #dir = dir_base * "/ffs"*string(i) * "/"
    dir = dir_base * "/fff"*string(i) * "/"
    MID.Solve.spectrum_from_file(dir, true)
end
#%%
#dvals, dfuncs = fss_convergence([10, 15], [2, 3], [1, 2], [1, 2], [1, 2], dir_base)
#dvals, dfuncs = ffs_convergence(Nfl, Nsl, ml, nl, zl, dir_base)
dvals, dfuncs = fff_convergence(Nfl, ml, nl, zl, dir_base)
#%%
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
