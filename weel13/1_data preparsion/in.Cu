# Melting
units 		metal
boundary 	p p p
atom_style 	atomic

lattice 		fcc 3.61
region 		box block 0 4 0 4 0 4
create_box 	1 box
create_atoms 	1 box
pair_style 	eam
pair_coeff 	* * "C:\Program Files\LAMMPS 64-bit 16Aug2018\Potentials\Cu_u3.eam"
 
timestep 	0.01
thermo 		1000

variable		 N equal step
variable 		pote equal pe
variable 		T equal temp
variable 		Press equal press
variable 		V equal vol
velocity 		all create 10 825577 dist gaussian

fix 		extra all print 100 "${N} ${T} ${V} ${pote} ${Press}" file data
dump 		1 all cfg 10000 a*.cfg mass type xs ys zs
dump_modify 	1 element Cu pad 6
dump        2 all xyz 10000 b*.xyz
dump        3 all dcd 10000 c.dcd
fix 		1 all npt temp 10 2000 1 iso 0 0 10
run 		120000
