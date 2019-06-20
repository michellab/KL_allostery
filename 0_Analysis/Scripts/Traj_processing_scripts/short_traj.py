# coding: utf-8
import mdtraj as md
traj1 = md.load_dcd("traj0001.dcd",top="topology.parm7")

traj_short_A = traj1[0:(len(traj)/2)]
traj_short_B = traj1[(len(traj)/2):(len(traj)+1)]
traj_short_A.save_dcd("short_traj_A.dcd")
traj_short_B.save_dcd("short_traj_B.dcd")
