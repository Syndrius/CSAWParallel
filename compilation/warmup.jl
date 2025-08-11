
using MIDParallel
using MID
using MPI
using PetscWrap
using SlepcWrap

dir = ARGS[1]
#dir = "/Users/matt/phd/MIDParallel/compilation/"

#going to have to do this completly differently!
#consiler 2 cases, base fff and qfm fff
#this will still miss diffferent q-profiles etc.
#but should hopefully cover most things.
Nr=4;
Nθ=2;
Nζ=2;
geo = init_geo(R0=4.0);
rgrid = init_grid(type=:rf, N=Nr, start=0.4, stop=0.6)
θgrid = init_grid(type=:af, N=Nθ, pf=1)
ζgrid = init_grid(type=:af, N=Nζ, pf=-1)

grids = init_grids(rgrid, θgrid, ζgrid);
isl = init_island(m0=3, n0=-2, A=0.01/3, flux=true, r0=0.5, qp=2.0)

prob = init_problem(q=cantori_q, geo=geo, type=:flux, isls=MID.IslandT[isl], met=:cylinder); 

#isl = init_island(m0=3, n0=-2, A=0.01/3, r0=0.5, qp=2.0)
#this tries to capture the other possibilities
prob2 = init_problem(q=cantori_q, geo=geo, isls=MID.IslandT[isl], type=:flux); 
surfs = construct_surfaces([(4, 3), (3, 2), (7, 5), (2, 1), (5, 4)], [0.5, 0.6, 0.58, 0.5, 0.5], prob, M=4, N=2)

solver = init_solver(nev=2, targets=[0.2, 0.5], prob=prob)

MPI.Init()

SlepcInitialize("-eps_nev 2 -st_type sinvert -eps_gen_hermitian")


W, I = MIDParallel.ParMatrix.preallocate_matrix(grids)

MIDParallel.par_construct(W, I, prob, grids)
MatAssemblyBegin(W, MAT_FINAL_ASSEMBLY)
MatAssemblyBegin(I, MAT_FINAL_ASSEMBLY)
MatAssemblyEnd(W, MAT_FINAL_ASSEMBLY)
MatAssemblyEnd(I, MAT_FINAL_ASSEMBLY)
destroy!(W)
destroy!(I)
W, I = MIDParallel.ParMatrix.preallocate_matrix(grids)
MIDParallel.par_construct(W, I, prob2, grids, surfs)

MatAssemblyBegin(W, MAT_FINAL_ASSEMBLY)
MatAssemblyBegin(I, MAT_FINAL_ASSEMBLY)
MatAssemblyEnd(W, MAT_FINAL_ASSEMBLY)
MatAssemblyEnd(I, MAT_FINAL_ASSEMBLY)

vecr, veci = MatCreateVecs(W)

nconv = MIDParallel.par_solve(W, I, solver, dir, vecr, veci)

destroy!(vecr)
destroy!(veci)
destroy!(W)
destroy!(I)
SlepcFinalize()
MPI.Finalize()
