
using MID
using Plots; #plotlyjs()
using DelimitedFiles
using Printf

#first we make sure ind=1 is the tae for each case
Alist=[9e-6, 1.5e-5, 4e-5, 6.5e-5, 9e-5, 1.5e-4, 4e-4, 6.5e-4, 9e-4, 1.5e-3, 4e-3]

#julia wants to add a .0 to floats, v annoying Alist may have to be a list of strings
file_ω_1 = @sprintf("data/small_convergence/9e-6/eigvals.dat")
file_ϕ_1 = @sprintf("data/small_convergence/9e-6/eigfuncs.dat")

N = 1000
rgrid = clustered_grid(N, 0.4, 0.6, 0.5)
grids = init_grids(rgrid=rgrid, mstart=-3, mcount=12, nstart=-6, ncount=3, nincr=4);


ω_1 = readdlm(file_ω_1, ',', ComplexF64)
ϕ_base_1 = readdlm(file_ϕ_1, ',', ComplexF64)
ϕ_1 = reconstruct_phi(ϕ_base_1, length(ω_1), N, grids.pmd.count, grids.tmd.count);


#fkn even 9e-6 island has distorted shit...
#this tae does seem to have no extra damping though!
tae_ind_1 = 1
plot_potential(r=rgrid, ϕ = ϕ_1, ind=tae_ind_1, pmd = grids.pmd, n=2)
display(ω_1[tae_ind_1])



file_ω_2 = @sprintf("data/small_convergence/1.5e-5/eigvals.dat")
file_ϕ_2 = @sprintf("data/small_convergence/1.5e-5/eigfuncs.dat")



ω_2 = readdlm(file_ω_2, ',', ComplexF64);
ϕ_base_2 = readdlm(file_ϕ_2, ',', ComplexF64);
ϕ_2 = reconstruct_phi(ϕ_base_2, length(ω_2), N, grids.pmd.count, grids.tmd.count);


#fkn even 9e-6 island has distorted shit...
#this tae does seem to have no extra damping though!
tae_ind_2 = 1
plot_potential(r=rgrid, ϕ = ϕ_2, ind=tae_ind_2, pmd = grids.pmd, n=2)
display(ω_2[tae_ind_2])


file_ω_3 = @sprintf("data/small_convergence/4e-5/eigvals.dat")
file_ϕ_3 = @sprintf("data/small_convergence/4e-5/eigfuncs.dat")

ω_3 = readdlm(file_ω_3, ',', ComplexF64);
ϕ_base_3 = readdlm(file_ϕ_3, ',', ComplexF64);
ϕ_3 = reconstruct_phi(ϕ_base_3, length(ω_3), N, grids.pmd.count, grids.tmd.count);


#fkn even 9e-6 island has distorted shit...
#this tae does seem to have no extra damping though!
tae_ind_3 = 1
plot_potential(r=rgrid, ϕ = ϕ_3, ind=tae_ind_3, pmd = grids.pmd, n=2)
display(ω_3[tae_ind_3])


file_ω_4 = @sprintf("data/small_convergence/6.5e-5/eigvals.dat")
file_ϕ_4 = @sprintf("data/small_convergence/6.5e-5/eigfuncs.dat")

ω_4 = readdlm(file_ω_4, ',', ComplexF64);
ϕ_base_4 = readdlm(file_ϕ_4, ',', ComplexF64);
ϕ_4 = reconstruct_phi(ϕ_base_4, length(ω_4), N, grids.pmd.count, grids.tmd.count);


#fkn even 9e-6 island has distorted shit...
#this tae does seem to have no extra damping though!
tae_ind_4 = 1
plot_potential(r=rgrid, ϕ = ϕ_4, ind=tae_ind_4, pmd = grids.pmd, n=2)
display(ω_4[tae_ind_4])


#this is the first instance of modification to the tae mode structure.
#although tae_ind=2 is not as modified, and has much less damping...
#think this will always be an issue.
file_ω_5 = @sprintf("data/small_convergence/9e-5/eigvals.dat")
file_ϕ_5 = @sprintf("data/small_convergence/9e-5/eigfuncs.dat")

ω_5 = readdlm(file_ω_5, ',', ComplexF64);
ϕ_base_5 = readdlm(file_ϕ_5, ',', ComplexF64);
ϕ_5 = reconstruct_phi(ϕ_base_5, length(ω_5), N, grids.pmd.count, grids.tmd.count);


#fkn even 9e-6 island has distorted shit...
#this tae does seem to have no extra damping though!
tae_ind_5 = 1
plot_potential(r=rgrid, ϕ = ϕ_5, ind=tae_ind_5, pmd = grids.pmd, n=2)
display(ω_5[tae_ind_5])


file_ω_6 = @sprintf("data/small_convergence/1.5e-4/eigvals.dat")
file_ϕ_6 = @sprintf("data/small_convergence/1.5e-4/eigfuncs.dat")

ω_6 = readdlm(file_ω_6, ',', ComplexF64);
ϕ_base_6 = readdlm(file_ϕ_6, ',', ComplexF64);
ϕ_6 = reconstruct_phi(ϕ_base_6, length(ω_6), N, grids.pmd.count, grids.tmd.count);


#fkn even 9e-6 island has distorted shit...
#this tae does seem to have no extra damping though!
tae_ind_6 = 1
plot_potential(r=rgrid, ϕ = ϕ_6, ind=tae_ind_6, pmd = grids.pmd, n=2)
display(ω_6[tae_ind_6])


#this one has significant distortion, interestingly, basically only the m=2,3 modes are relavant in this case?
#not certain higher mode numbers will actually help??
file_ω_7 = @sprintf("data/small_convergence/4e-4/eigvals.dat")
file_ϕ_7 = @sprintf("data/small_convergence/4e-4/eigfuncs.dat")

ω_7 = readdlm(file_ω_7, ',', ComplexF64);
ϕ_base_7 = readdlm(file_ϕ_7, ',', ComplexF64);
ϕ_7 = reconstruct_phi(ϕ_base_7, length(ω_7), N, grids.pmd.count, grids.tmd.count);


#fkn even 9e-6 island has distorted shit...
#this tae does seem to have no extra damping though!
tae_ind_7 = 1
plot_potential(r=rgrid, ϕ = ϕ_7, ind=tae_ind_7, pmd = grids.pmd, n=2)
display(ω_7[tae_ind_7])


#properly lost the plot now, tae is still there but it is cooked!
file_ω_8 = @sprintf("data/small_convergence/6.5e-4/eigvals.dat")
file_ϕ_8 = @sprintf("data/small_convergence/6.5e-4/eigfuncs.dat")

ω_8 = readdlm(file_ω_8, ',', ComplexF64);
ϕ_base_8 = readdlm(file_ϕ_8, ',', ComplexF64);
ϕ_8 = reconstruct_phi(ϕ_base_8, length(ω_8), N, grids.pmd.count, grids.tmd.count);


#fkn even 9e-6 island has distorted shit...
#this tae does seem to have no extra damping though!
tae_ind_8 = 5
plot_potential(r=rgrid, ϕ = ϕ_8, ind=tae_ind_8, pmd = grids.pmd, n=2)
display(ω_8[tae_ind_8])

#this is actually useless data.
file_ω_9 = @sprintf("data/small_convergence/9e-4/eigvals.dat")
file_ϕ_9 = @sprintf("data/small_convergence/9e-4/eigfuncs.dat")

ω_9 = readdlm(file_ω_9, ',', ComplexF64);
ϕ_base_9 = readdlm(file_ϕ_9, ',', ComplexF64);
ϕ_9 = reconstruct_phi(ϕ_base_9, length(ω_9), N, grids.pmd.count, grids.tmd.count);


#fkn even 9e-6 island has distorted shit...
#this tae does seem to have no extra damping though!
tae_ind_9 = 3
plot_potential(r=rgrid, ϕ = ϕ_9, ind=tae_ind_9, pmd = grids.pmd, n=2)
display(ω_9[tae_ind_9])

file_ω_10 = @sprintf("data/small_convergence/1.5e-3/eigvals.dat")
file_ϕ_10 = @sprintf("data/small_convergence/1.5e-3/eigfuncs.dat")

ω_10 = readdlm(file_ω_10, ',', ComplexF64);
ϕ_base_10 = readdlm(file_ϕ_10, ',', ComplexF64);
ϕ_10 = reconstruct_phi(ϕ_base_10, length(ω_10), N, grids.pmd.count, grids.tmd.count);


#fkn even 9e-6 island has distorted shit...
#this tae does seem to have no extra damping though!
tae_ind_10 = 1
plot_potential(r=rgrid, ϕ = ϕ_10, ind=tae_ind_10, pmd = grids.pmd, n=2)
display(ω_10[tae_ind_10])


ω_data = imag.([ω_1[tae_ind_1], ω_2[tae_ind_2], ω_3[tae_ind_3], ω_4[tae_ind_4], ω_5[tae_ind_5], ω_6[tae_ind_6], ω_7[tae_ind_7]])#, ω_8[tae_ind_8]])
Amps = Alist[1:7]
println(Amps)

scatter(Amps, ω_data, xaxis=:log, legend=false, ylabel="Im(ω)", xlabel="Island Amplitude", dpi=600)
savefig("data/small_convergence/damping_results.png")