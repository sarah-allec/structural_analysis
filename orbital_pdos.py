import numpy
from pylab import *
from matplotlib import rcParams
from matplotlib import axes

class OrbitalPDOS:
    def __init__(self,dos_file='DOSCAR',proj_file='projection'):
        self.dos_filename=dos_file
        self.proj_filename=proj_file

dos_file=open('DOSCAR')
data=dos_file.readlines()
dos_file.close()

proj_file=open('projection') #file specifying which atoms/orbitals to project onto
proj=proj_file.readlines()
proj_file.close()

# parse proj_file
for line in proj:
    if 'elements' in line:
        elements=line.split('=')[1].split('#')[0].split()
    if 'atom_numbers' in line:
        temp=line.split('=')[1].split('#')[0].split(',')
        temp=[t.split() for t in temp]
        temp=[i for sublist in temp for i in sublist]
        ranges=[i.find('-') for i in temp]
        ranges=[i for i, x in enumerate(ranges) if x == 1]
        atom_numbers=[]
        for i in range(len(temp)):
            if i in ranges:
                ssplit=temp[i].split('-')
                newlist=range(int(ssplit[0]),int(ssplit[1])+1)
                for n in newlist:
                    atom_numbers.append(n)
            else:
                atom_numbers.append(int(temp[i]))
    if 'orbitals' in line:
        orbitals=line.split('=')[1].split('#')[0].split()

n_atoms=int(data[0].split()[0])
n_bands=int(data[5].split()[2])
e_fermi=float(data[5].split()[3])

# first get total dos
total=[[float(i) for i in d.split()] for d in data[6:6+n_bands]]
energy=[d[0] for d in total]

spin=False
# check whether or not dos is spin-polarized
if len(total[0]) < 5:
    # not spin-polarized
    dos=[d[1] for d in total]
else:
    # spin-polarized
    spin=True
    dos_up=[d[1] for d in total]
    dos_down=[d[2] for d in total]

if elements[0] != 'All':
    if orbitals == 'None':
        el_dos=np.zeros((len(elements),n_bands))
    else:
        el_dos=np.zeros((len(elements)*len(orbitals),n_bands))

if atom_numbers != 'None':
    if orbitals == 'None':
        atno_dos=np.zeros((len(atom_numbers),n_bands))
    else:
        atno_dos=np.zeros((len(atom_numbers)*len(orbitals),n_bands))

for atom in range(1,2):
    start=atom*n_bands+atom+6
    curr_dos=data[start:start+n_bands]
     
