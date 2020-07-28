import multiprocessing as mp
import numpy as np
from math import floor, ceil
from ase import atoms
from ase.io import read

def volume( LCC ):
    ''' Computes the volume of the simulation box
    from the matrix of lattice vectors.

    args:
      LCC = (3x3 array) lattice vector matrix [a,b,c]

    returns:
      V = (scalar) the volume of the simulation box'''

    return np.dot(LCC[0],np.cross(LCC[1],LCC[2]))

def rho_0(N,V):
    ''' Computes the average density, rho_0,
    from number of atoms, N, and the volume, V,
    of the computational system.

    args:
      N = (scalar) number of atoms in simulation box
      V = (scalar) volume of simulation box
    returns:
      rho_0 = (scalar) the average density'''

    return N/V

def r_ij(LCC,positions):
    ''' Computes r_ij for a given timestep 

    args:
      LCC = (3x3 array) lattice vector matrix [a,b,c]
      positions = (n_atoms x 3 array) array of atomic coordinates

    returns:
      r_ij = (array) array of the shortest distances (across 
             periodic boundaries) of all pairs of atoms'''

    r = positions
    n_atoms = len(positions)        # number of atoms
    LCC_inv = np.linalg.inv(LCC)    # inverse of lattice vec matrix
    v_ij = [ np.abs(np.subtract(r[j], r[i])) for i in range(n_atoms) for j in range(i+1, n_atoms) ]
    
    vp_ij = [ np.mod( ( np.dot(LCC_inv,v) + 0.5), 1 ) - 0.5 for v in v_ij ] 
    vp_ij = [ np.matmul( LCC, vp) for vp in vp_ij ] 

    r_ij = np.linalg.norm(vp_ij,axis = 1)

    return r_ij 

def f(r_grid, r_ij, sigma):
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
    #a=1/(sigma*2*np.pi*r_ij)
    #b=-0.5/(sigma**2)
    return [a[i]*np.exp(b*(r_grid-r)**2) for i,r in enumerate(r_ij)]

def rho(cell,positions,ngrid,sigma):
    ''' Computes rho, the radial density, as a function of r.

    args:
        N = (scalar) number of atoms in simulation box
        r_grid = (1D array) the grid over which the Gaussian is computed
        r_ij = (array) the shortest distance(s) between the
               ith and jth atoms, accounting for PBC
        sigma = (scalar) width of Gaussian

    returns:
        rho = (array) the radial density as a function of r'''
    N = len(positions)
    curr_r_ij = r_ij(cell,positions) 
    r_grid = np.linspace(0.5,ceil(max(curr_r_ij)),ngrid)
    return [r_grid, np.divide( 2*np.sum( f(r_grid, curr_r_ij, sigma), axis = 0 ), N)]

def rdf(traj,sigma=0.15,ngrid = 100,n_cores = 8):
    ''' Computes R(r), the radial distribution function, and 
    g(r), the pair distribution function.

    args:
        sigma = (scalar) width of Gaussian

    returns:
        [r, g], where
        r = (array) the grid over which the Gaussian is compute
        g = (array) the pair distribution function for all element types
        R = (array) the radial distribution function for all element types'''

    n_steps = len(traj)
    N = len(traj[0].get_positions())
    g = np.zeros( (ngrid,) )
    R = np.zeros( (ngrid,) )

    groups = range(0, n_steps, n_cores)
    print(groups)
    procs = []
    for g in groups:
        for c in range(n_cores):
            step = g+c
            if step == n_steps:
                break
            cell = traj[step].get_cell()
            positions = traj[step].get_positions()
            procs.append(mp.Process(target=rho(), args=(cell, positions, ngrid, sigma )))
    #rho_r = np.divide(np.sum([ rho(V[i,:,:],positions[i,:,:],ngrid,sigma) for i in range(n_steps) ],axis=0), n_steps)
    #r = rho_r[0]
    #g = rho_r[1]/rho_0(N,V)
    #R = np.dot( np.dot(r,r), rho_r[1] )
    #R = 4*np.pi*R
    #g = g
    #return np.transpose([r, g, R])
