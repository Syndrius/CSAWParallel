
using MPI
using PetscWrap
using SlepcWrap
using Printf

#10 evals for n=21 case.
"""
vpr = 9.84932752388983 + 5.04483945186781e-16im
vpr = 39.154786963877065 + 5.38966773856915e-16im
vpr = 87.19478064930598 + 2.2256992758899416e-15im
vpr = 152.78640450004195 - 1.845578095746258e-15im
vpr = 234.31457505076125 - 3.7596016412481745e-14im
vpr = 329.77179816602114 + 4.402812081455945e-15im
vpr = 436.80760020836357 + 1.4753474753966e-15im
vpr = 552.7864045000414 - 7.36311607069819e-15im
vpr = 674.8524279678148 + 5.2708419014125524e-14im
vpr = 800.0000000000013 - 5.908524128911189e-14im
vpr = 925.1475720321853 - 3.790290844304628e-14im
vpr = 1047.2135954999567 + 3.170094058138009e-14im
vpr = 1163.192399791636 - 3.0868390117695046e-13im
vpr = 1270.2282018339783 + 1.1221582182953232e-13im
vpr = 1365.6854249492364 - 2.8985293219786704e-13im
vpr = 1447.2135954999558 - 1.8990843632643625e-13im
vpr = 1512.8052193506933 - 2.35883882988961e-13im
vpr = 1560.8452130361225 - 7.231451069447117e-14im
vpr = 1590.150672476107 - 1.9962677967026047e-12im
vpr = 2.9651701223551994e17 - 2.8175829413263468e16im
vpr = -3.718546636382165e25 - 4.861413453107018e24im
"""


MPI.Init()

n = 21
Δx = 1. / (n - 1)

#SlepcInitialize("-eps_target 0 -eps_nev 10 -st_pc_factor_shift_type NONZERO -st_type sinvert")

lowe = 2
highe = 95

#ok so this works perfectly here wot...
#so it seems like this algorithm is just a bit shit,
#or the quadrature rule is garbage, seems like the chebyshev quadrule is better, but then it cannot actually find any imaginary values...
#slepcargs = @sprintf("-st_pc_factor_shift_type NONZERO -eps_type ciss -rg_type interval -rg_interval_endpoints %s,%s,0,0 -eps_ciss_integration_points 148 -eps_ciss_quadrule chebyshev", lowe, highe)

#ellipse seems so much better tbh.
#so this will work in general, but we require a more specific integration region that we probably expected.
slepcargs = @sprintf("-st_pc_factor_shift_type NONZERO -eps_type ciss -rg_type ellipse -rg_ellipse_center 10 -rg_ellipse_radius 20 -rg_ellipse_vscale 1")# -eps_ciss_quadrule chebyshev")
#SlepcInitialize("-eps_type ciss -rg_type interval -rg_interval_endpoints 20,100,-1,1")

SlepcInitialize(slepcargs)

A = MatCreate()
B = MatCreate()
MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n)
MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, n, n)
MatSetFromOptions(A)
MatSetFromOptions(B)
MatSetUp(A)
MatSetUp(B)


A_rstart, A_rend = MatGetOwnershipRange(A)
B_rstart, B_rend = MatGetOwnershipRange(B)


for i in A_rstart:A_rend-1
    if (i == 0)
        MatSetValues(A, [0], [0, 1], [-2., 1] / Δx^2, INSERT_VALUES)
    elseif (i == n-1)
        MatSetValues(A, [n-1], [n-2, n-1], [1., -2.] / Δx^2, INSERT_VALUES)
    else
        MatSetValues(A, [i], i-1:i+1, [1., -2., 1.] / Δx^2, INSERT_VALUES)
    end
end


for i in B_rstart:B_rend-1
    MatSetValue(B, i, i, -1., INSERT_VALUES)
end

(A_rstart == 0) && MatSetValues(A, [0], [0,1], [1., 0.], INSERT_VALUES)
(B_rstart == 0) && MatSetValue(B, 0, 0, 0., INSERT_VALUES)

(A_rend == n) && MatSetValues(A, [n-1], [n-2,n-1], [0., 1.], INSERT_VALUES)
(B_rend == n) && MatSetValue(B, n-1, n-1, 0., INSERT_VALUES)


MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY)
MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY)
MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY)
MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY)

eps = EPSCreate()
EPSSetOperators(eps, A, B)
EPSSetFromOptions(eps)
EPSSetUp(eps)


EPSSolve(eps)


nconv = EPSGetConverged(eps)

if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    for ieig in 0:nconv - 1
        vpr, vpi = EPSGetEigenvalue(eps, ieig)
        @show (vpr)
    end
end

MatDestroy(A)
MatDestroy(B)
EPSDestroy(eps)

SlepcFinalize()

MPI.Finalize()
