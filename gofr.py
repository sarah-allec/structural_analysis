import numpy as np
from math import floor, ceil
from ase import atoms
from ase.io import read

class GofR():
    ''' This is a class for computing the pair distribution function of a given
    structure.

    args:
      traj = (ASE atoms object) trajectory/structure
    '''
    def __init__( self, traj ):
        self.traj = traj
        self.elements = self.traj[0].get_chemical_symbols()
        self.species = list(set(self.elements))
        self.spec_indices = {}
        for s in self.species:
            self.spec_indices[s] = [i for i,v in enumerate(self.elements) if v == s]
        self.N = 0

    def volume( self, LCC ):
        ''' Computes the volume of the simulation box
        from the matrix of lattice vectors.

        args:
          LCC = (3x3 array) lattice vector matrix [a,b,c]

        returns:
          V = (scalar) the volume of the simulation box'''

        return np.dot(LCC[0],np.cross(LCC[1],LCC[2]))

    def rho_0( self, N, V ):
        ''' Computes the average density, rho_0,
        from number of atoms, N, and the volume, V,
        of the computational system.
    
        args:
          N = (scalar) number of atoms in simulation box
          V = (scalar) volume of simulation box
        returns:
          rho_0 = (scalar) the average density'''
    
        return N/V

    def r_ij( self, LCC, positions, indices1, indices2 ):
        ''' Computes r_ij for a given timestep 
    
        args:
          LCC = (3x3 array) lattice vector matrix [a,b,c]
          positions = (n_atoms x 3 array) array of atomic coordinates
    
        returns:
          r_ij = (array) array of the shortest distances (across 
                 periodic boundaries) of all pairs of atoms'''
    
        r = positions

        if indices1 != None and indices2 != None:
            v_ij = [ np.abs(np.subtract(r[j], r[i])) for i in indices1 for j in indices2 if j != i ]
        else:
            n_atoms = len(positions)        # number of atoms
            v_ij = [ np.abs(np.subtract(r[j], r[i])) for i in range(n_atoms) for j in range(i+1, n_atoms) ]

        LCC_inv = np.linalg.inv(LCC)    # inverse of lattice vec matrix
        vp_ij = [ np.mod( ( np.dot(LCC_inv,v) + 0.5), 1 ) - 0.5 for v in v_ij ] 
        vp_ij = [ np.matmul( LCC, vp) for vp in vp_ij ] 
    
        r_ij = np.linalg.norm(vp_ij,axis = 1)
    
        return r_ij 

    def f( self, r_grid, r_ij, sigma ):
        ''' Computes f, the atomic contributions to the radial density.
    
        args:
          r_grid = (1D array) the grid over which the Gaussian is computed
          r_ij = (array) the shortest distance(s) between the 
                 ith and jth atoms, accounting for PBC
          sigma = (scalar) width of Gaussian
    
        returns:
          f = (array) the atomic contribution(s) to the radial density'''
    
        pi_3_2 = np.power( np.sqrt(np.pi), 3)
        a = 1 / ( 4 * pi_3_2 * sigma * np.square(r_ij) )
        b = - 1 / (sigma * sigma)
        return [a[i]*np.exp(b*(r_grid-r)**2) for i,r in enumerate(r_ij)]

    def rho( self, cell, positions, ngrid, sigma, indices1=None, indices2=None ):
        ''' Computes rho, the radial density, as a function of r.
    
        args:
            N = (scalar) number of atoms in simulation box
            r_grid = (1D array) the grid over which the Gaussian is computed
            r_ij = (array) the shortest distance(s) between the
                   ith and jth atoms, accounting for PBC
            sigma = (scalar) width of Gaussian
    
        returns:
            rho = (array) the radial density as a function of r'''
        
        curr_r_ij = self.r_ij(cell,positions,indices1,indices2) 
        r_grid = np.linspace(0.5,8.0,ngrid)
        if indices1 == None or indices2 == None:
            return [r_grid, np.divide( 2*np.sum( self.f(r_grid, curr_r_ij, sigma), axis = 0 ), self.N )]
        else:
            N1 = len(indices1)
            N2 = len(indices2)
            return [r_grid, np.divide( self.N*np.sum( self.f(r_grid, curr_r_ij, sigma), axis = 0 ), (N1*N2) )]

    def little_g( self, step, sigma=0.15, ngrid = 100, el1=None, el2=None ):
        ''' Computes R(r), the radial distribution function, and 
        g(r), the pair distribution function.
    
        args:
            sigma = (scalar) width of Gaussian
    
        returns:
            [r, g], where
            r = (array) the grid over which the Gaussian is compute
            g = (array) the pair distribution function for all element types
            R = (array) the radial distribution function for all element types'''

        positions = self.traj[step].get_positions()
        cell = self.traj[step].get_cell()
        self.N = len( positions ) 
   
        V = self.volume(cell)
       
        if el1 != None and el2 != None:
            indices1 = self.spec_indices[el1]
            indices2 = self.spec_indices[el2]
            pair = el1+"-"+el2
            coeff = self.N / (len(indices1)*len(indices2))
            rho_r = self.rho(cell,positions,ngrid,sigma,indices1,indices2)

        else:
            coeff = self.N
            rho_r = self.rho(cell,positions,ngrid,sigma)
        r = rho_r[0] 
        g = rho_r[1]/self.rho_0(self.N,V)
        return np.transpose([r, g])
