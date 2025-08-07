
#MIDParallel is not working for large evals in Helmholtz case that are easily found in serial
#this is true with a single proc, implying it is not the grid division that is wrong.
using MID
using MIDParallel
using FFTW

include("/home/149/mt3516/island_damping/MID/Convergence/FEMConvergence/solution.jl") #going to be annoying af

#%%
dir = "/scratch/y08/mt3516/Helmholtz/test"

Nfs = [6]
mtarg = [1, 2]
ntarg = [1, 2]
zrtarg = [1, 2]
i = 1

flr = MID.Structures.FLRT(δ=1e-8)
geo = init_geo(R0=1.0)
prob = MID.Helmholtz.TestProblemT(flr=flr, geo=geo)
an_ev, code = anal_evals(mtarg, ntarg, zrtarg)
display(an_ev)
display(sqrt.(an_ev))
solver = init_solver(nev=5, targets=sqrt.(an_ev), prob=prob)
#solver = init_solver(nev=5, targets=[7.29], prob=prob)

xgrid = init_grid(type=:rf, N = Nfs[i])#, sep1=0.45, sep2=0.55, frac=0.4)
ygrid = init_grid(type=:af, N=Nfs[i]+1)
zgrid = init_grid(type=:af, N=Nfs[i]-1)


grids = init_grids(xgrid, ygrid, zgrid)
mkpath(dir * "/fff")
inputs_to_file(prob=prob, grids=grids, solver=solver, dir=dir*"/fff/")
#%%
par_post_process(dir*"/fff/")
