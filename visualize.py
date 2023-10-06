from RiboGraphViz import RGV
from arnie.bpps import bpps
from arnie.mfe import mfe
import numpy as np
import matplotlib.pyplot as plt
import sys

seq = sys.argv[1]
name = sys.argv[2]

mRNA = seq

struct = mfe(mRNA, package='vienna')

# vector of p(unpaired)
bpp_vec = 1 - np.sum(bpps(mRNA, package='vienna'),axis=0)

rg = RGV(struct)

plt.figure(figsize=(12,12))
rg.draw(c=bpp_vec, cmap='plasma')
plt.savefig(f'{name}.png')
