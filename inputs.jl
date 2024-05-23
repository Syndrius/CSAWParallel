

#inputs are cooking it from file, try to find out why...


using MID

N = 100;
#the collect is a bit annoying, but ok because we will typically use a clustered grid.
rgrid = collect(LinRange(0, 1, N));

geo = GeoParamsT(R0=10.0)

prob = init_problem(q=island_damping_q, geo=geo); #probbaly should use geo if it is part of prob,
#even if it is not really used.
grids = init_grids(rgrid=rgrid, mstart=2, mcount=2, nstart=-2, ncount=1);

inputs_to_file(prob=prob, grids=grids, dir="data/")

prob3, grids3 = inputs_from_file(dir="data/small_convergence/4e-4/")

problem_to_file(prob=prob, filename="data/prob_test.txt")
prob2 = problem_from_file(filename="data/prob_test.txt")


grids_to_file(grids=grids, filename="data/grids_test.txt")

grids2 = grids_from_file(filename="data/grids_test.txt")

ω, ϕ = construct_and_solve(prob=prob3, grids=grids3, full_spectrum=true);


#what are the inputs for this!!!, everything needs kwargs!!!
rcont, omcont, col = reconstruct_cont(ω = ω, ϕ = ϕ, rgrid = rgrid, pmd = grids.pmd, tmd = grids.tmd);

scatter(rcont, omcont, group=col, ylimits=(-0.05, 1.05))


tae_ind = find_ind(ω, 0.381)

plot_potential(r=rgrid, ϕ=ϕ, ind=tae_ind, pmd=grids.pmd, n=1)

grids_cont = init_grids(rgrid=rgrid[2:end], mstart=2, mcount=2, nstart=-2, ncount=1);

ω_cont = continuum(prob=prob, grids=grids_cont)

scatter(rgrid[2:end], ω_cont)