import numpy as np

set retain_order,1
set pdb_retain_ids,1

cmd.load("first_frame.pdb")

cmd.alter('all','chain="A"')

cmd.save("first_frame_ALL_CHAIN_A.pdb")
