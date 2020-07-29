from ase.io import read
import numpy as np
import pandas as pd

flatten = lambda l: [item for sublist in l for item in sublist]

class Polyhedra():
    ''' This is a class for computing the polyhedra network of
    a material. The methods in this class provide the number of 
    each type of polyhedra and the average distance between them.

    args:
      traj = (ASE atoms object) trajectory/structure 
    '''
    def __init__(self, traj):
        self.traj = traj
        self.elements = self.traj[0].get_chemical_symbols()
        self.species = list(set(self.elements))
        self.spec_indices = {}
        for s in self.species:
            self.spec_indices[s] = [i for i,v in enumerate(self.elements) if v == s]
    def get_shell(self,step, el1, el2, cutoff1, cutoff2):
        indices1 = self.spec_indices[el1] #indices for the central element
        indices2 = self.spec_indices[el2] #indices for the coordinating element

        r = self.traj[step].get_positions()
        LCC = self.traj[step].get_cell()
        n_atoms = len(r)        # number of atoms
        LCC_inv = np.linalg.inv(LCC)    # inverse of lattice vec matrix
        ij = [ [i] + [j] for i in indices1 for j in indices2 ]
        v_ij = [ [np.abs(np.subtract(r[j], r[i]))] for i in indices1 for j in indices2 ]
        vp_ij = [ np.mod( ( np.dot(LCC_inv,np.transpose(v)) + 0.5), 1 ) - 0.5 for v in v_ij ]
        vp_ij = [ np.matmul( LCC, vp) for vp in vp_ij ]

        r_ij = np.linalg.norm(vp_ij,axis = 1) 
        sorted_indices=np.argsort(np.array(r_ij), axis=None)
        r_ij=np.sort(r_ij,axis=None)
        print(r_ij[0:10])     
        shell_indices = pd.DataFrame([ ij[sorted_indices[i]] for i,d in enumerate(r_ij) if d >= cutoff1 and d <= cutoff2 ],
                        columns = ['col1','col2'])
        shell_distances = [ d for d in r_ij if d > cutoff1 and d < cutoff2 ]
 
        shell_indices = shell_indices.groupby('col1')['col2'].apply(list)
        return [shell_indices, shell_distances]

my_traj = read('CONTCAR_big',index=':')
my_poly = Polyhedra(my_traj)
[i,r]=my_poly.get_shell(0,'In','In',0.0,2.4)
pd.DataFrame( i ).to_csv( "second_shell.dat",sep="\t",header=False,index=True )


