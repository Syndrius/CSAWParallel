function fff_convergence_inputs(Nfs, mtarg, ntarg, zrtarg, dir)

    flr = MID.Structures.FLRT(δ=1e-8)
    geo = init_geo(R0=1.0)
    prob = MID.Helmholtz.TestProblemT(flr=flr, geo=geo)
    an_ev, code = anal_evals(mtarg, ntarg, zrtarg)
    solver = init_solver(nev=10, targets=sqrt.(an_ev), prob=prob)

    for i in 1:length(Nfs)
        xgrid = init_grid(type=:rf, N = Nfs[i])#, sep1=0.45, sep2=0.55, frac=0.4)
        ygrid = init_grid(type=:af, N=Nfs[i]+1)
        zgrid = init_grid(type=:af, N=Nfs[i]-1)
        

        grids = init_grids(xgrid, ygrid, zgrid)
        mkpath(dir * "/fff"*string(i))
        inputs_to_file(prob=prob, grids=grids, solver=solver, dir=dir*"/fff"*string(i) *"/")
    end

end
function fff_convergence(Nfs, mtarg, ntarg, zrtarg, dir_base)

    #assume Nfs and Nsf are hte same length, just the size of each grid.

    an_ev, code = anal_evals(mtarg, ntarg, zrtarg)

    eval_diff = zeros(length(Nfs), length(an_ev))
    efunc_diff = zeros(length(Nfs), length(an_ev))

    Nint = 2 * Nfs[end]
    int_grid = LinRange(0, 1, Nint)

    an_sol = zeros(length(an_ev), Nint)
    intphi = zeros(ComplexF64, length(an_ev), Nint, Nint, Nint) #needs to be complex unfort
    
    #gets the radial part of each eigenvalue.
    #comparing the angular parts is more effort than worth.
    for i in 1:length(an_ev)
        an_sol[i, :] = rad_anal_sol(code[i][1], code[i][2], int_grid)
    end

    for i in 1:length(Nfs)
        dir = dir_base * "/fff" * string(i) * "/"
        #xgrid = init_grid(type=:rf, N = Nfs[i])#, gp=3)
        #ygrid = init_grid(type=:af, N=Nfs[i])
        #zgrid = init_grid(type=:af, N=Nfs[i])
        
        #grids = init_grids(xgrid, ygrid, zgrid)
        _, grids, _  = inputs_from_file(dir=dir)
        evals = evals_from_file(dir=dir)

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
                if zero_crossings(ϕft[:, m+1, n+1, 1], Nfs[i]) == zr

                    push!(eval_true, ind)
                    push!(mn_true, mn[j])
                    break
                end
            end
        end
        #display("or here")
        display(sort(an_ev))
        display(sort(real.(evals.ω[eval_true]) .^2))
        eval_diff[i, :] = abs.(sort(real.(evals.ω[eval_true])) .^2 .- sort(an_ev))

        #xgp = MID.inst_grid(xgrid)
        #ygp = MID.inst_grid(ygrid)
        #zgp = MID.inst_grid(zgrid)
        xgp, ygp, zgp = MID.inst_grids(grids)
        for k in 1:length(eval_true), j in 1:Nint, l in 1:Nint, p in 1:Nint
            #will only work if m,n are positive incrementing from 1.
            #now m, n will be extra cooked.
            #ah yes, we have -1, need some kind of ind_to_mode map.
            #no longer actually needed as we need to interpolate in 2d.
            #mtrue = mode_to_ind(mn_true[k][1], ygrid)
            #mn_true for the last bit doesn't make any sense!
            #Feels like this should be complex???
            #guess if there is no pf? -> could signify a problemo
            ϕ = efunc_from_file(dir=dir, ind=eval_true[k], ft=false) 
            intphi[k, j, l, p] = MID.Mapping.hermite_interpolation(int_grid[j], int_grid[l], int_grid[p], ϕ, xgp, ygp, zgp)
        end

        intphi = fft(intphi, [3, 4])

        for k in 1:length(eval_true)
            #ntrue = mn_true[k][2] #don't think this is needed, just because we are essentially picking a random part.
            #love this.
            mtrue, ntrue = mn_true[k]
            if mtrue < 0
                mind = Nint + mtrue + 1
            else
                mind = mtrue + 1
            end
            if ntrue < 0
                nind = Nint + ntrue + 1
            else
                nind = ntrue + 1
            end
            am1 = argmax(abs.(an_sol[k, :]))
            am2 = argmax(abs.(real.(intphi[k, :, mind, nind])))
            sf = an_sol[k, am1] / real(intphi[k, am2, mind, nind])
            efunc_diff[i, k] = sum(abs.(an_sol[k, :] .- sf .* real.(intphi[k, :, mind, nind])))
        end

    end
    return eval_diff, efunc_diff

end
