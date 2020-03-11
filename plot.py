import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import axes

filename=sys.argv[1] # make default
file=open(filename)
data=[float(i) for i in file.readlines()]
file.close()

label=sys.argv[2] # make default

#plots
params = {
    'axes.labelsize': 14,
    'legend.fontsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'text.usetex': True,
    'figure.figsize': [5.0, 5.0]
    }
rcParams.update(params)

# get units
if label == "Energy":
    units = "eV"
if label == "Pressure":
    units = "kb"
if label == "Temperature":
    units = "K"

#plt.style.use('dark_background')
fig,ax=plt.subplots()
ax.plot(range(len(data)),data)
ax.set_xlabel("Timestep")
ax.set_ylabel(label + " (" + units+")")
ax.legend(loc="best")
fig.tight_layout()
fig.savefig(label+".png",dpi=600)

