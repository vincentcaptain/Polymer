# LAMMPS script on the mean squared displacement of 10-mer
# settings
# a 3D 60 box
variable      x equal 60
variable      y equal 30
variable      z equal 30
# temperature t = 1
variable      t equal 1
#variable      rho equal 0.7
neigh_modify  delay 0 check yes one 200000 page 2000000


# System setup

units         lj
dimension     3
atom_style    bond
bond_style    harmonic
pair_style    lj/cut 3
#processors    1 4 6
read_data     data.nowall

velocity      all create $t 97287

group         polymer type 1

fix           1 polymer nve
fix           2 polymer langevin $t $t 1 1209


timestep      0.01
# equilibration run
thermo        1000
run           10000

reset_timestep 0
compute       com polymer com
compute       msd polymer msd com yes
compute       g polymer gyration
run 5000
dump          pos polymer xyz 10 pol.xyz
thermo_style custom c_g c_msd[4] c_com[1] 
thermo        10
run 1000000
