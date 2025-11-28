
using MID
using MIDParallel
using Plots; gr()
using Plots; plotlyjs()
#using MIDViz
#%%
#ideally we want to do the qfm and non-qfm case here.
#actually, probably won't even have this example here
#midparallel should be close to identical
#just need to read and write from grids.
#so we can probably just have a qfm and non-qfm example.

#so this file is just for getting the island mapping to work.

#this seems like an adequate choice.
isl = init_island(m0=1, n0=-1, w=0.1, qp=1.0, ψ0=1/2)
geo = init_geometry(:cyl, R0=1.0)
fields = init_fields(q=island_q, isl=isl)
prob = init_problem(geometry=geo, fields=fields)

#%%

ψgrid = init_grid(:ψ, 30, sep1=0.45, sep2=0.55, frac=0.5)
θgrid = init_grid(:θ, 5, pf=0)
φgrid = init_grid(:φ, 3, pf=0)

grids = init_grids(ψgrid, θgrid, φgrid)
#%%

solver = init_solver(prob=prob, target=0.0, nev=100)
#%%

dir_base = "/Users/matt/phd/MIDParallel/data/example/"
inputs_to_file(prob=prob, grids=grids, solver=solver, dir=dir_base);

par_post_process(dir_base)
evals = evals_from_file(dir_base);

continuum_plot(evals, xlimits=(0.45, 0.55), ylimits=(0.0, 0.1))

scatter(evals.x1, real.(evals.ω))

#got a bit of structure tbh!
#dang daniel we can actually get an island mode. wowee.
ind = find_ind(evals, 0.03339)

ϕft = efunc_from_file(dir_base, ind);

p1grid, p2grid, _ = MID.inst_grids(grids)
#%%
plot(p1grid, real.(ϕft[:, 1, 1]))
plot!(p1grid, real.(ϕft[:, 1, 2]))
plot!(p1grid, real.(ϕft[:, 1, 3]))
plot!(p1grid, real.(ϕft[:, 2, 1]))
plot!(p1grid, real.(ϕft[:, 2, 2]))
plot!(p1grid, real.(ϕft[:, 2, 3]))
plot!(p1grid, real.(ϕft[:, 3, 1]))
plot!(p1grid, real.(ϕft[:, 3, 2]))
plot!(p1grid, real.(ϕft[:, 3, 3]))
plot!(p1grid, real.(ϕft[:, 4, 1]))
plot!(p1grid, real.(ϕft[:, 4, 2]))
plot!(p1grid, real.(ϕft[:, 4, 3]))
plot!(p1grid, real.(ϕft[:, 5, 1]))
plot!(p1grid, real.(ϕft[:, 5, 2]))
plot!(p1grid, real.(ϕft[:, 5, 3]))
#%%
ϕ = efunc_from_file(dir_base, ind, ft=false);
contourf(p2grid, p1grid, real.(ϕ[:, :, 1]))
#%%

κgrid = init_grid(:ψ, 40, stop=0.999)
ᾱgrid = init_grid(:θ, 20)
τgrid = init_grid(:φ, 5)
isl_grids = init_grids(κgrid, ᾱgrid, τgrid)
#%%
qfm_spectrum_to_isl(dir_base, isl_grids, dir_base*"w1_surfaces.jld2")
isl_evals = evals_from_file(dir_base * "isl_map/")
ind = find_ind(isl_evals, 0.03339)
ϕft = efunc_from_file(dir_base*"/isl_map/", ind);

harmonic_plot(ϕ_isl, isl_grids, ind)
#%%
p1grid, p2grid, _ = MID.inst_grids(isl_grids)
plot(p1grid, real.(ϕft[:, 1, 1]))
plot!(p1grid, real.(ϕft[:, 1, 2]))
plot!(p1grid, real.(ϕft[:, 1, 3]))
plot!(p1grid, real.(ϕft[:, 2, 1]))
plot!(p1grid, real.(ϕft[:, 2, 2]))
plot!(p1grid, real.(ϕft[:, 2, 3]))
plot!(p1grid, real.(ϕft[:, 3, 1]))
plot!(p1grid, real.(ϕft[:, 3, 2]))
plot!(p1grid, real.(ϕft[:, 3, 3]))
plot!(p1grid, real.(ϕft[:, 4, 1]))
plot!(p1grid, real.(ϕft[:, 4, 2]))
plot!(p1grid, real.(ϕft[:, 4, 3]))
plot!(p1grid, real.(ϕft[:, 5, 1]))
plot!(p1grid, real.(ϕft[:, 5, 2]))
plot!(p1grid, real.(ϕft[:, 5, 3]))

#%%
using MIDCantori
isl = init_island(m0=1, n0=-1, w=0.1, qp=1.0, ψ0=1/2)
geo = init_geometry(:cyl, R0=1.0)
fields = init_fields(q=island_q, isl=isl)
prob = init_problem(geometry=geo, fields=fields)

rats = lowest_rationals(8, island_q(0.0, isl)[1], island_q(1.0, isl)[1])
gl = surface_guess(rats, prob.fields.q)

surfs = construct_surfaces(rats, gl, prob);

surfaces_to_file(surfs, dir_base*"w1_surfaces.jld2")
