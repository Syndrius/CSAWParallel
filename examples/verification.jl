
#verification with Bowden singular

using MID
using MIDParallel

using Plots; plotlyjs()


#first we define the problem and write to file.
#this is identical to MID.
#see what kind of mem is required for this!
#think the matrix is v sparse now, so hopefully memory is not as big a problemo.
Nr=4000;
Nθ=100
#rgrid = collect(LinRange(0, 1, N));
geo = GeoParamsT(R0=10.0);
prob = init_problem(q=bowden_singular_q, geo=geo, dens=bowden_singular_dens, δ=-4e-9, met=diagonal_toroidal_metric!); 
rgrid = init_fem_grid(N=Nr)
θgrid = init_fem_grid(N=Nθ, pf=1)
ζgrid = init_sm_grid(start=-1, count = 1)
grids = init_grids(rgrid, θgrid, ζgrid);


#looks like full path is needed... a bit annoying tbh.
#dir_base = "/home/149/mt3516/island_damping/MIDParallel/data/example/"
dir_base = "/scratch/y08/mt3516/ffs_verification/"

mkpath(dir_base)

dir = dir_base * "4000x100/"

mkdir(dir)

inputs_to_file(prob=prob, grids=grids, dir=dir);

ω = par_vals_from_file(dir*"vals.dat", geo.R0);

display(ω[1])

display(imag(ω[1])/real(ω[1]))

#then the efuncs
ϕ = par_func_from_file(dir*"funcs.dat", 1, grids);

ϕms = mode_structure(ϕ, grids);

#maybe we want an imag flag for this? i.e to match fig 5 of Bowden.
plot_potential(ϕms, grids, 1, 1)


#now we consider a proper convergence test by running through N and δ
#we will do this twice, once for Nθ=50 and once for Nθ=100 (assuming that is not too much to ask)
#may want to repeat this with off-diagonal just as a comparison.

#Nθ case

Nθ = 100

θgrid = init_fem_grid(N=Nθ, pf=1)
ζgrid = init_sm_grid(start=-1, count = 1)
geo = GeoParamsT(R0=10.0);

dir_base = "/scratch/y08/mt3516/ffs_verification/theta100/"

mkpath(dir_base)

Nlist = [500, 1000, 1500, 2000, 3000, 4000]
δlist = [-4.0e-7, -4.0e-8, -4.0e-9, -4.0e-10, -4.0e-11]

for N in Nlist

    rgrid = init_fem_grid(N=N)

    grids = init_grids(rgrid, θgrid, ζgrid)

    for δ in δlist

        prob = init_problem(q=bowden_singular_q, geo=geo, dens=bowden_singular_dens, δ=δ, met=diagonal_toroidal_metric!);

        dlab = parse(Int64, split(string(δ), "-")[end])

        dir = dir_base * "N_" * string(N) * "_" * "delta" * string(dlab) * "/"
        mkdir(dir)

        inputs_to_file(prob=prob, grids=grids, dir=dir)
    end
end


#now lets see if we can make a verification plot!

Nθ = 50

θgrid = init_fem_grid(N=Nθ, pf=1)
ζgrid = init_sm_grid(start=-1, count = 1)
geo = GeoParamsT(R0=10.0);

dir_base = "/scratch/y08/mt3516/ffs_verification/theta50/"


Nlist = [500, 1000, 1500, 2000, 3000, 4000]
δlist = [-4.0e-7, -4.0e-8, -4.0e-9, -4.0e-10, -4.0e-11]

ωlist = zeros(ComplexF64, length(Nlist), length(δlist));

for (i, N) in enumerate(Nlist)

    rgrid = init_fem_grid(N=N)

    grids = init_grids(rgrid, θgrid, ζgrid)

    for (j, δ) in enumerate(δlist)

        prob = init_problem(q=bowden_singular_q, geo=geo, dens=bowden_singular_dens, δ=δ, met=diagonal_toroidal_metric!);

        dlab = parse(Int64, split(string(δ), "-")[end])

        dir = dir_base * "N_" * string(N) * "_" * "delta" * string(dlab) * "/"
        
        #just trust the ev's to be number 1?

        ω = par_vals_from_file(dir * "vals.dat", geo.R0)

        ωlist[i, j] = ω[1]


    end
end


p = scatter()
#pretty close, but still not perf :(
#may need to check if the tae is actually the first one? perhap??
#will still need to compare with Axel's case, and check example with the off-diagonal terms.
#this converges to 0.018345, %diff of ~4.17, not so bad.
for i in 1:length(δlist)
    scatter!(Nlist, imag.(ωlist[:, i]) ./ real.(ωlist[:, 1]))
end

display(p)
display(ωlist[4, 3])