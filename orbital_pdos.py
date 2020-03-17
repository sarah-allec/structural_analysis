import numpy
from pylab import *
from matplotlib import rcParams
from matplotlib import axes

class OrbitalPDOS:
    def __init__(self,dos_file='DOSCAR',proj_file='projection'):
        self.dos_filename=dos_file
        self.proj_filename=proj_file

        self.dos_data=[]
        self.energy=[]
        self.total_dos=[]
        self.dos_up=[]
        self.dos_down=[]
        self.n_atoms=0
        self.n_bands=0
        self.e_fermi=0
        self.spin=False

        self.proj_data=[]
        self.elements=[]
        self.atom_numbers=[]
        self.orbitals=[]

    def read_dos(self):
        file=open(self.dos_filename)
        self.dos_data=file.readlines()
        file.close()
        self.n_atoms=int(self.dos_data[0].split()[0])
        self.n_bands=int(self.dos_data[5].split()[2])
        self.e_fermi=float(self.dos_data[5].split()[3])

    def read_proj(self):
        file=open(self.proj_filename)
        self.proj_data=file.readlines()
        file.close()

    def parse_proj(self):
        # parse proj_file
        for line in self.proj_data:
            if 'elements' in line:
                self.elements=line.split('=')[1].split('#')[0].split()
            if 'atom_numbers' in line:
                if 'None' in line:
                    self.atom_numbers.append('None')
                else:
                    temp=line.split('=')[1].split('#')[0].split(',')
                    temp=[t.split() for t in temp]
                    temp=[i for sublist in temp for i in sublist]
                    ranges=[i.find('-') for i in temp]
                    ranges=[i for i, x in enumerate(ranges) if x == 1]
                    for i in range(len(temp)):
                        if i in ranges:
                            ssplit=temp[i].split('-')
                            newlist=range(int(ssplit[0]),int(ssplit[1])+1)
                            for n in newlist:
                                self.atom_numbers.append(n)
                        else:
                            self.atom_numbers.append(int(temp[i]))
            if 'orbitals' in line:
                self.orbitals=line.split('=')[1].split('#')[0].split()

    def get_total_dos(self):
        
        total=[[float(i) for i in d.split()] for d in self.dos_data[6:6+self.n_bands]]
        self.energy=[d[0]-self.e_fermi for d in total]

        self.total_dos=[d[1] for d in total]
        if len(total[0]) > 4:
            # spin-polarized
            self.spin=True
            self.dos_up=[d[1] for d in total]
            self.dos_down=[d[2] for d in total]

    #def.get_pdos(self):

    def plot_dos(self,energy=None,dos=None,labels=None):
        if energy == None:
            energy=self.energy
        if dos == None:
            dos=self.total_dos
        params = {
            'axes.labelsize': 14,
            'legend.fontsize': 12,
            'xtick.labelsize': 8,
            'ytick.labelsize': 8,
            'text.usetex': False,
            'figure.figsize': [5.0, 5.0]
            }
        rcParams.update(params)
        fig,ax=plt.subplots()

        if labels==None:
            ax.plot(energy,dos)
            ax.set_xlabel("Energy (eV)")
            ax.set_ylabel("Density of states")
            fig.tight_layout()
            fig.savefig("total_dos.png",dpi=600)
        else: 
            for l in range(len(labels)):
                ax.plot(energy,dos[l],label=labels[l])
            ax.legend(loc="best")
            ax.set_xlabel("Energy (eV)")
            ax.set_ylabel("Density of states")
            fig.tight_layout()
            fig.savefig("pdos.png",dpi=600)

my_pdos = OrbitalPDOS()
my_pdos.read_proj()
my_pdos.parse_proj()
#my_pdos.read_dos()
#my_pdos.get_total_dos()
#my_pdos.plot_dos()

#if elements[0] != 'All':
#    if orbitals == 'None':
#        el_dos=np.zeros((len(elements),n_bands))
#    else:
#        el_dos=np.zeros((len(elements)*len(orbitals),n_bands))
#
#if atom_numbers != 'None':
#    if orbitals == 'None':
#        atno_dos=np.zeros((len(atom_numbers),n_bands))
#    else:
#        atno_dos=np.zeros((len(atom_numbers)*len(orbitals),n_bands))
#
#for atom in range(1,2):
#    start=atom*n_bands+atom+6
#    curr_dos=data[start:start+n_bands]
#     
