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
        self.r_grid = 0
        self.N = 0
        self.V = 0 

        self.sigma = 0
        self.r_ij_vec = 0
        self.f_vec = 0
        self.g_vec = 0

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
    
        self.r_ij_vec = np.linalg.norm(vp_ij,axis = 1) 
        return self.r_ij_vec 

    def f( self, r_grid=None, r_ij=None, sigma=None ):
        ''' Computes f, the atomic contributions to the radial density.
    
        args:
          r_grid = (1D array) the grid over which the Gaussian is computed
          r_ij = (array) the shortest distance(s) between the 
                 ith and jth atoms, accounting for PBC
          sigma = (scalar) width of Gaussian
    
        returns:
          f = (array) the atomic contribution(s) to the radial density'''
    
        if r_grid is None:
            r_grid = self.r_grid
        if r_ij is None:
            r_ij = self.r_ij_vec
        if sigma is None:
            sigma = self.sigma

        pi_3_2 = np.power( np.sqrt(np.pi), 3)
        a = 1 / ( 4 * pi_3_2 * sigma * np.square(r_ij) )
        b = - 1 / (sigma * sigma)
        self.f_vec = [a[i]*np.exp(b*(r_grid-r)**2) for i,r in enumerate(r_ij)]
        return self.f_vec

    def df( self ):
        # this function should be used when little_g has been used with the current
        # instance of the class. 
        # To do: overload this function to be used without having ran little_g
        coeff = (2 / ( self.sigma * self.sigma ) )
        numerator = [ coeff * (r - self.r_grid) for r in self.r_ij_vec ]
        df_dr = np.sum( np.multiply( numerator, self.f_vec ), axis = 0 )
        return df_dr

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
        self.r_grid = np.linspace(0.5,8.0,ngrid)
        if indices1 == None or indices2 == None:
            return [self.r_grid, np.divide( 2*np.sum( self.f(), axis = 0 ), self.N )]
        else:
            N1 = len(indices1)
            N2 = len(indices2)
            return [self.r_grid, np.divide( self.N*np.sum( self.f(), axis = 0 ), (N1*N2) )]

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

        self.sigma = sigma
        positions = self.traj[step].get_positions()
        cell = self.traj[step].get_cell()
        self.N = len( positions ) 
   
        self.V = self.volume(cell)
       
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
        g = rho_r[1]/self.rho_0(self.N,self.V)
        self.g_vec =  np.transpose([r, g])
        return self.g_vec
