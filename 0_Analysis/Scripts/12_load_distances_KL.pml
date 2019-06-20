import numpy as np

atoms_indices_KLs = "0A_1B_atoms_indices_KLs.dat"

cmd.load("../../0_TRAJECTORIES/0_system_A/first_frame.pdb")
infile = np.loadtxt(atoms_indices_KLs)
cmd.hide("lines")
cmd.show("cartoon")
cmd.show("dashes")
cmd.set("cartoon_color","grey40")

# cut out any with KL < 20% max ?? - done
# assign colours to range of KL values - done
# hide labels  -DONE
# KL values are 3rd column of the input file - done

# Could add to atom and assign a value to "most moving" atoms?

# Some colors for range of values: 
dash_colour_list = ["white" , "paleyellow" , "orange" , "red", "firebrick"]

#print infile
count = 0
for i in infile:\
	atom1 = int(i[0])  \
	atom2 = int(i[1]) \
	kl_value = float(i[2]) \
	if kl_value < 14.8: \
		continue \	
	if kl_value < 14.9: \
		colour_scale = dash_colour_list[0] \
	elif (kl_value >= 14.9) and (kl_value < 14.95): \
		colour_scale = dash_colour_list[1] \
	elif (kl_value >=15.0) and (kl_value < 15.05): \
		colour_scale = dash_colour_list[2] \
	elif (kl_value >= 15.05) and (kl_value < 15.1): \
		colour_scale = dash_colour_list[3] \
	elif (kl_value >=15.1): \
		colour_scale = dash_colour_list[4] \
	filename = "distance_%s" % count \
	cmd.distance(filename,"rank %d" % atom1, "rank %d" % atom2) \
	cmd.hide("labels",filename) \
	cmd.set("dash_color",colour_scale,filename) \
	count = count + 1

cmd.set("dash_transparency","0.6")
