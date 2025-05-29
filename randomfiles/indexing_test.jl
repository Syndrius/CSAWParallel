#something may be wrong with the ol indexing!
using MIDParallel
using MID
using MIDViz
using Plots; plotlyjs()
using Plots; gr()
using LinearAlgebra
#%%

geo = init_geo(R0=10.0)
isl = init_island(m0=2, n0=-1, w=0.03, r0=00.5, qp=2.0)

prob = MID.Structures.init_isl_problem(geo=geo, isl=isl)

prob = init_problem(geo=geo, q=fu_dam_q)

#%%
#10x3x1 reproduces error
#this is not present in the non-island case, implying that the problem may not be indexing?
#so error only occurs when pf=0 for both cases
#looks like this may just be a problem in the slepc alg
#which is also happens to occur in qfm case.
rgrid = init_grid(type=:rf, N=20, start=0.0, stop=0.999)#, left_bc=false)
θgrid = init_grid(type=:af, N=5, pf=0)
ζgrid = init_grid(type=:af, N=2, pf=0)

grids = init_grids(rgrid, θgrid, ζgrid)
#%%
solver = init_solver(prob=prob, target=0.00000, nev=300)
#%%
dir = "/Users/matt/phd/MIDParallel/data/example/"
inputs_to_file(dir=dir, prob=prob, grids=grids, solver=solver)

#%%
par_post_process(dir)
#%%
evals = evals_from_file(dir=dir);

continuum_plot(evals, ymax=0.1)
display(length(evals.ω))
#%%
evals, _, _ = compute_spectrum(prob=prob, grids=grids, solver=solver)

W, I = MID.Construct.construct(prob, grids)

display(W[4268, 4268])
display(I[4268, 4268])

#unsure if they are supposed to be positive def, but they are not. Could be aproblemo
isposdef(Matrix(W))
isposdef(I)
#defs supposed to be hermitian...
ishermitian(Matrix(W))
ishermitian(Matrix(I))
ishermitian(I)

for i in 1:size(W)[1]
    display(I[i, i])
end
println(real.(Matrix(W)[43, 43]))
println(real.(Matrix(I)[43, 43]))

res = Matrix(W) / Matrix(I)

eigvals(res)

println(res[43, 43])

spy(real.(W))
display(W[32:36, 40])

display(W[14, 14])
display(I[14, 15])
display(W[22, 46])
display(W[39, 48])
println(W[121, :])
println(I[121, :])

size(W)
MID.Indexing.index_to_grid(13, grids)
MID.Indexing.index_to_grid(48, grids)

#%%
function read_petsc_mat(dir)

    f = open(dir, "r")

    n = countlines(f) - 2 #probably only works for 1 proc.
    close(f)

    f = open(dir, "r")
    mat = zeros(n, n) #idealy this will actually be sparse

    for (i, line) in enumerate(readlines(f)[3:end])#first two has mpi info

        substr = split(line, ":")
        row = parse(Int64, split(substr[1], " ")[2])
        #display(substr[2])
        #display(strip(substr[2]))
        tups = filter(x -> !isspace(x), substr[2])

        for tup in split(tups, ")")[1:end-1]
            #display(tup)
            data = split(tup, ",")
            col = parse(Int64, data[1][2:end])
            val = parse(Float64, data[2][1:end])
            #display(col)
            #display(val)
            mat[row+1, col+1] = val
        end

        #data = (Array{Tuple}, substr[2])

        #display(line)
    end
    close(f)
    return mat

end
Ip = read_petsc_mat("/Users/matt/phd/I.dat")
Wp = read_petsc_mat("/Users/matt/phd/W.dat")

#%%
display(Wp[4266:4270, 4266:4270])
display(W[4266:4270, 4266:4270])
display(Ip[4266:4270, 4266:4270])
display(I[4266:4270, 4266:4270])

#these matrices look identical, we may need to actually compare them more specifically
#otherwise, we are going to have to look into the lu method to understand why row 4267 is causing problemos, doesn't really look like a prpblem with our code.
#we probably also need to check if our matrices are supposed to be positive definite, because thay currently are not.
#same as checking for Hermitian.


