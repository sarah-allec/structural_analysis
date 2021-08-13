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
        shell_indices = pd.DataFrame([ ij[sorted_indices[i]] for i,d in enumerate(r_ij) if d >= cutoff1 and d <= cutoff2 ],
                        columns = ['col1','col2'])
        shell_distances = [ d for d in r_ij if d > cutoff1 and d < cutoff2 ]
 
        shell_indices = shell_indices.groupby('col1')['col2'].apply(list)
        return [shell_indices, shell_distances]

    def get_polyhedra(self,step,el1,el2,cutoff1,cutoff2,cutoff3):
        ''' Computes the number of each type of polyhedra (corner-, edge-, and face-sharing)
        for each atom of type el1. 

        args:
          step = (scalar) timestep to compute polyhedra for
          el1 = (string) element type (e.g. "In") for which to compute polyhedra
          el2 = (string) element type (e.g. "O") in first coord. shell of el1
          cutoff1 = (scalar) lower cutoff distance for first shell (usually 0.0)
          cutoff2 = (scalar) upper cutoff distance for first shell (can be obtained as first minimum
                    in el1-el2 g(r))
          cutoff3 = (scalar) upper cutoff for computing polyhedra - e.g. for In2O3, corresponds
                    to the largest distance in the third shell around In (can be obtained as second 
                    minimum in el1-el1 g(r))'''

        # first step is to get first shell of el
        first_shell = self.get_shell(step, el1, el2, cutoff1, cutoff2)[0]
        in_shells = self.get_shell(step, el1, el1, cutoff2, cutoff3)[0]

        n_el1 = len(self.spec_indices[el1])
        corner = np.zeros( n_el1 )
        edge = np.zeros( n_el1 )
        face = np.zeros( n_el1 )
        corner_indices = []
        edge_indices = []
        face_indices = []
        # loop through each In atom
        for i in range(n_el1):
            O_atoms = first_shell[i]
            c_indices_temp = []
            e_indices_temp = []
            f_indices_temp = []
            for j in in_shells[i]:
                shared_O_atoms = [ v for v in O_atoms if v in first_shell[j] ]
                if len(shared_O_atoms) == 1:
                    c_indices_temp.append(j) 
                    corner[i] += 1
                if len(shared_O_atoms) == 2:
                    e_indices_temp.append(j)
                    edge[i] += 1
                if len(shared_O_atoms) == 3:
                    f_indices_temp.append(j)
                    face[i] += 1
            corner_indices.append([i]+c_indices_temp)
            edge_indices.append([i]+e_indices_temp)
            face_indices.append([i]+f_indices_temp)
        return [ np.column_stack( (corner, edge, face) ), [ corner_indices, edge_indices, face_indices ] ]

