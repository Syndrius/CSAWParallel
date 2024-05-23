#!/bin/bash

DIR="data/small_convergence/"
#start with test of the different island sizes
#then we will ramp up the mcount and ncount
#until we can fix slepc implementation, we need to start small, hopefully it will be enough to show damping rate increases with island size.
R0=10.0
Alist=(9e-6 1.5e-5 4e-5 6.5e-5 9e-5 1.5e-4 4e-4 6.5e-4 9e-4 1.5e-3 4e-3)
N=1000
delta=-4e-9
tae_freq=0.0014653584 #taken from new_q_test, probbaly not accurate enough!
#simple no island case.
#first we define the inputs and write the inputs to file. 
#this is a simple case
#this is done in serial.
for i in ${!Alist[@]}; do
    dir_str=$DIR${Alist[$i]}"/"
    eval mkdir $dir_str
    inputs_str='
    using MID;
    qprof = island_damping_q;
    R0 = '"${R0}"';
    isl = IslandT(A='"${Alist[$i]}"', m0=5, n0=4)
    geo = GeoParamsT(R0=R0);
    prob = init_problem(q=qprof, geo=geo, isl=isl);
    N = '"${N}"';
    rgrid = clustered_grid(N, 0.4, 0.6, 0.5)
    grids = init_grids(rgrid=rgrid, mstart=-3, mcount=12, nstart=-6, ncount=3, nincr=4);
    inputs_to_file(prob=prob, grids=grids, dir='\""${dir_str}"\"');
    '
    #echo $inputs_str
    #julia -e "${inputs_str}"
    julia amplitude_convergence.jl ${Alist[$i]} ${dir_str}

    #now construct the matrix, which is done in parallel
    mpiexecjl -n 4 julia construct_from_inputs.jl $dir_str

    #then solve the matrix, in serial 
    julia solve_from_inputs.jl $dir_str $tae_freq



done
#so we have a copy of the scipt used to generate the data.
eval cp convergence.sh $DIR
