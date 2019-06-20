import numpy as np
cmd.load("../../0_TRAJECTORIES/0_system_A/first_frame.pdb","PC1_atom_contribution")
cmd.load("../../0_TRAJECTORIES/0_system_A/first_frame.pdb","PC2_atom_contribution")
cmd.loadall("*.pdb")

cmd.hide("lines")
cmd.show("cartoon")
cmd.bg_color("grey70")

infile1 = open("PC1_atom_contribution.dat", "r")
infile2 = open("PC2_atom_contribution.dat", "r")

PC1_bfact = []
PC2_bfact = []

for line in infile1.readlines():\ 
	PC1_bfact.append(float(line))\
infile1.close()

for line in infile2.readlines():\ 
	PC2_bfact.append(float(line))\
infile2.close()

cmd.alter("PC1", 'b=0.0')
cmd.alter("PC1 and n. CA", "b=PC1_bfact.pop(0)")
cmd.alter("PC2", 'b=0.0')
cmd.alter("PC2 and n. CA", "b=PC2_bfact.pop(0)")

cmd.spectrum("b", "white_grey60_grey20_orange_firebrick", "PC1")
cmd.spectrum("b", "white_grey60_grey20_orange_firebrick", "PC2")

cmd.save("PCA.pse")
