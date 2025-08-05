

function fss_convergence_inputs(Nfs, Nss, mtarg, ntarg, zrtarg, dir)

    flr = MID.Structures.FLRT(δ=1e-8)
    geo = init_geo(R0=1.0)
    prob = MID.Helmholtz.TestProblemT(flr=flr, geo=geo)
    an_ev, code = anal_evals(mtarg, ntarg, zrtarg) #access to this will be annoying af.
    solver = init_solver(nev=20, targets=sqrt.(an_ev), prob=prob)

    for i in 1:length(Nfs)
        xgrid = init_grid(type=:rf, N = Nfs[i])#, gp=3)
        ygrid = init_grid(type=:as, start=1, N=Nss[i])
        zgrid = init_grid(type=:as, start=1, N=Nss[i])
        
        grids = init_grids(xgrid, ygrid, zgrid)
        mkpath(dir * "/fss"*string(i))
        inputs_to_file(prob=prob, grids=grids, solver=solver, dir=dir * "/fss" * string(i) * "/")
    end
end


function fss_convergence(Nfs, Nss, mtarg, ntarg, zrtarg, dir_base)

    an_ev, code = anal_evals(mtarg, ntarg, zrtarg)

    eval_diff = zeros(length(Nfs), length(an_ev))
    efunc_diff = zeros(length(Nfs), length(an_ev))

    Nint = 2 * Nfs[end]
    int_grid = LinRange(0, 1, Nint)
    an_sol = zeros(length(an_ev), Nint)
    intphi = zeros(ComplexF64, length(an_ev), Nint) #needs to be complex unfort
    #gets the radial part of each eigenvalue.
    #comparing the angular parts is more effort than worth.
    for i in 1:length(an_ev)
        an_sol[i, :] = rad_anal_sol(code[i][1], code[i][2], int_grid)
    end

    for i in 1:length(Nfs)
        dir = dir_base * "/fss" * string(i) * "/"
        xgrid = init_grid(type=:rf, N = Nfs[i])#, gp=3)
        ygrid = init_grid(type=:as, start=1, N=Nss[i])
        zgrid = init_grid(type=:as, start=1, N=Nss[i])
        
        #could just read from file, but recreating is fine.
        grids = init_grids(xgrid, ygrid, zgrid)
        evals = evals_from_file(dir=dir)
        #don't need this with slice solver
        #nbcs = length(unique(MID.Indexing.compute_boundary_inds(grids)))

        eval_true = []
        mn_true = []
        for zr in zrtarg, m in mtarg, n in ntarg
            mat = []
            mn = []
            for j in 1:length(evals.ω)
                if abs((abs(evals.ω[j]) - 1)) < 0.2 #bc solution, so we ignore.
                    continue
                end
                if (m, n) == (abs(evals.modelabs[j][1]), abs(evals.modelabs[j][2]))
                    push!(mat, j)
                    push!(mn, evals.modelabs[j])
                end
            end
            for (j, ind) in enumerate(mat)
                ϕft = efunc_from_file(dir=dir, ind=ind)
                #only works with m, n starting from 1.
                if zero_crossings(ϕft[:, m, n, 1], Nfs[i]) == zr

                    push!(eval_true, ind)
                    push!(mn_true, mn[j])
                    break
                end
            end
        end
        eval_diff[i, :] = abs.(sort(real.(evals.ω[eval_true])) .^2 .- sort(an_ev))
        #now lets see if we can compare the efuncs as well.
        #probably easiest to just consider the radial part.
        #for this we will interpolate the solution, as the non-interpolated parts should match exactly.

        xgp = MID.inst_grid(xgrid)
        for k in 1:length(eval_true), j in 1:Nint
            #will only work if m,n are positive incrementing from 1.
            ϕft = efunc_from_file(dir=dir, ind=eval_true[k])
            intphi[k, j] = MID.Mapping.hermite_interpolation(int_grid[j], ϕft[:, mn_true[k][1], mn_true[k][2], :], xgp)
        end

        for k in 1:length(eval_true)
            am1 = argmax(abs.(an_sol[k, :]))
            am2 = argmax(abs.(real.(intphi[k, :])))
            sf = an_sol[k, am1] / real(intphi[k, am2])
            efunc_diff[i, k] = sum(abs.(an_sol[k, :] .- sf .* real.(intphi[k, :])))
        end

    end
    return eval_diff, efunc_diff
end
