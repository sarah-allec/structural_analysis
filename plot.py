import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import axes

#plots
params = {
    'font.family': 'serif',
    'font.serif': 'Times',
    'axes.labelsize': 14,
    'axes.xmargin': 0,
    'axes.ymargin': 0,
    'legend.fontsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'text.usetex': True,
    'figure.figsize': [5.0, 5.0]
    }
rcParams.update(params)

fig,ax=plt.subplots()
ax.plot(r_crystal,g_crystal,"-",color="black",label=r'$c-\mathrm{In}_{2}\mathrm{O}_{3}$')
ax.plot(r_amorphous,g_amorphous,"-",color="dodgerblue",label=r'$a-\mathrm{In}_{2}\mathrm{O}_{3}$')

ax.set_xlabel(r'r ($\mathrm{\AA}$)')
ax.set_ylabel("g(r)")
ax.legend(loc="best")
fig.tight_layout()
fig.savefig("gofr_in2o3.png",dpi=600)

