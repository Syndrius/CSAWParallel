
using MID
using MIDParallel
using MIDViz
#%%

#first see if we can get some damn island modes in real small caseo
geo = init_geo(R0=1.0)
isl = init_island(m0=1, n0=-1, w=0.1, ψ0=0.5, qp=1.0, flux=true)


#need to change this
#prob = MID.Structures.init_isl_problem(geo=geo, isl=isl)
prob = init_problem(geo=geo, q=island_q, isl=isl, met=:cylinder)
#%%

rgrid = init_grid(type=:rf, N=50, sep1=0.45, sep2=0.55, frac=0.5)
#θgrid = asm_grid(start=-4, N=9)
θgrid = init_grid(type=:af, N=15, pf=1)
ζgrid = init_grid(type=:af, N=4, pf=-1)

grids = init_grids(rgrid, θgrid, ζgrid)
#%%
solver = init_solver(prob=prob, target=0.0, nev=200)
#%%
#evals, ϕ, ϕft = compute_spectrum(prob=prob, grids=grids, solver=solver);

dir = "/Users/matt/phd/MIDParallel/data/mapping/"

inputs_to_file(dir=dir, grids=grids, solver=solver, prob=prob)

#%%
par_post_process(dir)
#%%
Nκ = 51
Nᾱ = 16
Nτ = 5

κgrid = init_grid(type=:rf, N=Nκ, stop=0.999)
ᾱgrid = init_grid(type=:af, N=Nᾱ)
τgrid = init_grid(type=:af, N=Nτ)

isl_grids = init_grids(κgrid, ᾱgrid, τgrid)

MID.Mapping.qfm_spectrum_to_isl(dir, isl_grids, "/Users/matt/phd/MID/data/surfaces/cantori_island/w1_surfs.jld2")




#%%
evals = evals_from_file(dir=dir);

isl_om = evals.ω[0.4 .<evals.x1 .<0.6]
continuum_plot(evals)
#%%
#33 is a bit islandy, 35 is ok
#the 30's are ok ish
isl_ind = 36

ind = find_ind(evals, isl_om[isl_ind])

ϕft = efunc_from_file(dir=dir, ind=ind);

potential_plot(ϕft, grids, label_max=0.5)
#%%
ϕ = efunc_from_file(dir=dir, ind=ind, ft=false);
MIDViz.Plotting.contour_plot(ϕ, grids, 1)
#%%
#the data is not really islandy, but probably good enough to test the mapping out.

κgrid = init_grid(type=:rf, N=200, stop=0.999)
ᾱgrid = init_grid(type=:af, N=80)
τgrid = init_grid(type=:af, N=30)
isl_grids = init_grids(κgrid, ᾱgrid, τgrid)

#%%
#naturally, this is slow af.
#this is perhaps a bit to slow, I guess it is the interpolation that is the problemo?
#we may need to try and make an interpolation object storing the polynomials or something, 
#ok so interpolation is defs the problemo. Need to construct an interpolation object!
#our method is unnacetably slow. Unsure how to fix tbh.
#perhaps we use the normal interpolation for testing, then just use our shit version
#once we have good results? And just cop the huge time sink
#periodic stuff will be annoying af though
#this is just slow af unfort. Ceebs fixing tbh.
#looks like interpolations and hermite are giving the same thing, at least for this small shite version
@time MID.PostProcessing.tor_spectrum_to_isl(dir, isl_grids)

isl_evals = evals_from_file(dir=dir*"isl_map/");

continuum_plot(isl_evals)
