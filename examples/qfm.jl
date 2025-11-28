#this example with the chaos system seems to be working ok.
using MID
using MIDParallel
using Plots
#%%
geo = init_geometry(:tor, R0=4.0)

#probably cant just use w huh.
#this would require the root solve shite.
isl = init_island(m0=3, n0=1, A=0.1)

fields = init_fields(:ψ, q=cantori_q, isl=isl)

prob = init_problem(geometry=geo, fields=fields)

#%%

sgrid = init_grid(:s, 30, start=0.15, stop=0.9)
ϑgrid = init_grid(:ϑ, 6, pf=1)
ζgrid = init_grid(:ζ, 1, pf=-1)
#sgrid = init_grid(:s, 100, start=0.15, stop=0.9)
#ϑgrid = init_grid(:sm, 2, start=1)
#ζgrid = init_grid(:sm, 1, start=-1)
#sgrid = init_grid(:s, 40, start=0.15, stop=0.9)
#ϑgrid = init_grid(:ϑ, 10, pf=1)
#ζgrid = init_grid(:sm, 1, start=-1)

grids = init_grids(sgrid, ϑgrid, ζgrid)
#%%

#maybe need a guard against this!
solver = init_solver(prob=prob, full_spectrum=true)
solver = init_solver(prob=prob, targets=[0.2, 0.3, 0.4], nev=100)
#%%


dir_base = "/Users/matt/phd/MIDParallel/data/example/"
inputs_to_file(prob=prob, grids=grids, solver=solver, dir=dir_base);

par_post_process(dir_base) #unfort have we have to do this!
#now we can read the data in. first the eigenvalues,
evals = evals_from_file(dir_base);


scatter(evals.x1, real.(evals.ω))
ind = find_ind(evals, 0.3)
ϕft = efunc_from_file(dir_base, ind);

pgrid = MID.inst_grid(sgrid)
plot(pgrid, real.(ϕft[:, 1, 1]))
plot!(pgrid, real.(ϕft[:, 2, 1]))
plot!(pgrid, real.(ϕft[:, 3, 1]))
#############
ψgrid = init_grid(:ψ, 80, start=0.25, stop=0.8)
θgrid = init_grid(:θ, 30)
φgrid = init_grid(:φ, 10)
tor_grids = init_grids(ψgrid, θgrid, φgrid)
#%%
#pretty sure this only works for fff.
qfm_spectrum_to_tor(dir_base, tor_grids, "/Users/matt/phd/MID/test/data/benchmark_surfaces.jld2");
size(ϕ)

evals_tor = evals_from_file(joinpath(dir_base, "tor_map/"));

ϕ_tor = efunc_from_file(joinpath(dir_base, "tor_map/"), 33, ft=false);

p1grid = MID.inst_grid(ψgrid)
p2grid = MID.inst_grid(θgrid)

#looks kind of shite
#but shows that the mapping is still working!
contourf(p2grid, p1grid, real.(ϕ_tor[:, :, 1]))


continuum_plot(tor_evals)
tor_ind = find_ind(tor_evals, 0.25)
display(ϕft_tor)
#think these harmonic plots are fkn stoopid.
harmonic_plot(ϕft_tor, tor_grids, tae_ind, label_max=0.5)
harmonic_plot(ϕft_tor, tor_grids, cont_ind, label_max=0.5)
#this is kind of cool
#maybe shows us that our qfm choice was a bit cooked.
#i.e. the perturbation was perhaps a bit large.
contour_plot(ϕ, grids, cont_ind)
contour_plot(ϕ_tor, tor_grids, cont_ind)
