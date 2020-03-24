import numpy
from pylab import *
from matplotlib import rcParams
from matplotlib import axes

class OrbitalPDOS:
    def __init__(self,dos_file='DOSCAR',proj_file='projection'):
        self.dos_filename=dos_file
        self.proj_filename=proj_file
        self.elements=[]

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
        self.species=[]
        self.atom_numbers=[]
        self.orbitals=[]

    def read_dos(self):
        file=open(self.dos_filename)
        self.dos_data=file.readlines()
        file.close()
        self.n_atoms=int(self.dos_data[0].split()[0])
        self.n_bands=int(self.dos_data[5].split()[2])
        self.e_fermi=float(self.dos_data[5].split()[3])

    def read_poscar(self):
        file=open('POSCAR')
        data=file.readlines()
        file.close()

        species=data[5].split()
        n_per=data[6].split()
        self.elements = [item for sublist in ([species[i]]*int(v) for i,v in enumerate(n_per)) for item in sublist]

    def read_proj(self):
        file=open(self.proj_filename)
        self.proj_data=file.readlines()
        file.close()

    def parse_proj(self):
        # parse proj_file
        for line in self.proj_data:
            if 'species' in line:
                self.species=line.split('=')[1].split('#')[0].split()
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

    def get_pdos(self):
        el_dos=[]
        atno_dos=[]
        if self.species[0] != 'All':
            if self.orbitals[0] == 'None':
                el_dos=np.zeros((len(self.species),self.n_bands))
            else:
                el_dos=np.zeros((len(self.species)*len(self.orbitals),self.n_bands))
        
        if self.atom_numbers[0] != 'None':
            if self.orbitals[0] == 'None':
                atno_dos=np.zeros((len(self.atom_numbers),self.n_bands))
            else:
                atno_dos=np.zeros((len(self.atom_numbers)*len(self.orbitals),self.n_bands))

        for atom in range(1,self.n_atoms+1):
        #for atom in range(25,27):
            start=atom*self.n_bands+atom+6
            curr_dos=self.dos_data[start:start+self.n_bands]
            if len(el_dos) != 0:
                el_index=self.species.index(self.elements[atom-1])
                if self.orbitals[0] == 'None':
                    for row in range(len(curr_dos)):
                        el_dos[el_index][row]=el_dos[el_index][row]+sum([float(r) for r in curr_dos[row].split()[1:]])

                else:
                    for row in range(len(curr_dos)):
                        if 's' in self.orbitals:
                            el_dos[3*el_index][row]=el_dos[3*el_index][row]+sum([float(r) for r in curr_dos[row].split()[1:3]])
                        if 'p' in self.orbitals:
                            el_dos[3*el_index+1][row]=el_dos[3*el_index+1][row]+sum([float(r) for r in curr_dos[row].split()[3:9]])                          
                        if 'd' in self.orbitals:
                            el_dos[3*el_index+2][row]=el_dos[3*el_index+2][row]+sum([float(r) for r in curr_dos[row].split()[9:]])                          
        return el_dos

    def plot_dos(self,energy=None,dos=None,labels=None):
        if energy == None:
            energy=self.energy
            print('using default energy')
        if dos == None:
            print('using total dos')
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

my_pdos = OrbitalPDOS(dos_file='DOSCAR_ismear0')
my_pdos.read_poscar()
my_pdos.read_proj()
my_pdos.parse_proj()
my_pdos.read_dos()
my_pdos.get_total_dos()
pdos=my_pdos.get_pdos()
my_pdos.plot_dos(dos=pdos,labels=['Co(s)','Co(p)','Co(d)','S(s)','S(p)'])
     
