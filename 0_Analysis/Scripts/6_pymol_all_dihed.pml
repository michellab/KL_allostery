import numpy as np

##### first set KLs #####

cmd.load("../0_TRAJECTORIES/0_system_A/first_frame.pdb","KL_backbone")
cmd.load("../0_TRAJECTORIES/0_system_A/first_frame.pdb","KL_sidechain")
cmd.load("../0_TRAJECTORIES/0_system_A/first_frame.pdb","KL_psi")
cmd.load("../0_TRAJECTORIES/0_system_A/first_frame.pdb","KL_phi")
cmd.load("../0_TRAJECTORIES/0_system_A/first_frame.pdb","KL_chi1")
cmd.load("../0_TRAJECTORIES/0_system_A/first_frame.pdb","KL_chi2")

cmd.hide("lines")
cmd.show("cartoon")
cmd.bg_color("grey70")

infile1 = open("KL_backbone.dat")
infile2 = open("KL_sidechain.dat")
infile3 = open("KL_psi.dat")
infile4 = open("KL_phi.dat")
infile5 = open("KL_chi1.dat")
infile6 = open("KL_chi2.dat")

KL_backbone = []
KL_sidechain = []
KL_psi = []
KL_phi = []
KL_chi1 = []
KL_chi2 = []

for line in infile1.readlines():\
	KL_backbone.append(float(line))\	
infile1.close()

for line in infile2.readlines():\
	KL_sidechain.append(float(line))\	
infile1.close()

for line in infile3.readlines():\
	KL_psi.append(float(line))\	
infile3.close()

for line in infile4.readlines():\
	KL_phi.append(float(line))\	
infile4.close()

for line in infile5.readlines():\
	KL_chi1.append(float(line))\	
infile5.close()

for line in infile6.readlines():\
	KL_chi2.append(float(line))\	
infile6.close()

cmd.alter('KL_backbone','b=0.0')
cmd.alter('KL_backbone and n. CA','b=KL_backbone.pop(0)')
cmd.alter('KL_sidechain','b=0.0')
cmd.alter('KL_sidechain and n. CA','b=KL_sidechain.pop(0)')
cmd.alter('KL_psi','b=0.0')
cmd.alter('KL_psi and n. CA','b=KL_psi.pop(0)')
cmd.alter('KL_phi','b=0.0')
cmd.alter('KL_phi and n. CA','b=KL_phi.pop(0)')
cmd.alter('KL_chi1','b=0.0')
cmd.alter('KL_chi1 and n. CA','b=KL_chi1.pop(0)')
cmd.alter('KL_chi2','b=0.0')
cmd.alter('KL_chi2 and n. CA','b=KL_chi2.pop(0)')

# To colour sidechain and backbone separately (not same scale as KL values might vary a lot)
cmd.spectrum("b", "white_orange_firebrick", "KL_backbone")
cmd.spectrum("b", "white_orange_firebrick", "KL_sidechain")
cmd.spectrum("b", "white_orange_firebrick", "KL_psi")
cmd.spectrum("b", "white_orange_firebrick", "KL_phi")
cmd.spectrum("b", "white_orange_firebrick", "KL_chi1")
cmd.spectrum("b", "white_orange_firebrick", "KL_chi2")

# To colour all to same scale, comment out the cmd.spectrum commands above, and uncomment below:
# cmd.spectrum("b", "white_orange_firebrick")

cmd.save("Dihedral_KL.pse")
