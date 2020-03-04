import numpy as np
import sys

def angle(v1,v2):
    return np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))

def volume(a,b,c):
    return np.dot(a,np.cross(b,c))    

coords_filename=sys.argv[1] #should be XDATCAR; need to make a default for this
xyz_filename=sys.argv[2] #Custom name for xyz file; need to make a default for this

coords_file=open(coords_filename)
coords_lines=coords_file.readlines()
coords_file.close()

species = coords_lines[5].split()
n_per = coords_lines[6].split()
n_atoms = sum(int(i) for i in n_per)
n_steps = len(coords_lines)/(n_atoms+8)

#get list of elements (only needs to be done once)
species = coords_lines[5].split()
n_per = coords_lines[6].split()
elements = [int(n_per[i])*[s] for i,s in enumerate(species)]
elements = [item for sublist in elements for item in sublist]

with open(xyz_filename,'w') as new:
    for step in range(n_steps):
        #get lattice vectors
        a = [float(i) for i in coords_lines[step*(n_atoms+8)+2].split()]
        b = [float(i) for i in coords_lines[step*(n_atoms+8)+3].split()]
        c = [float(i) for i in coords_lines[step*(n_atoms+8)+4].split()]
        a_mag=np.linalg.norm(a)
        b_mag=np.linalg.norm(b)
        c_mag=np.linalg.norm(c)
        alpha=angle(b,c)
        beta=angle(a,c)
        gamma=angle(a,b)
        omega=volume(a,b,c)
   
        new.write(str(n_atoms)+'\n')
        new.write('Step '+str(step+1)+'\n')
        #get coordinates
        for n in range(n_atoms):
            coord = [float(i) for i in coords_lines[step*(n_atoms+8)+8+n].split()]
            u = coord[0]
            v = coord[1]
            w = coord[2]
            x=str(a_mag*u+b_mag*np.cos(gamma)*v+c_mag*np.cos(beta)*w)
            y=str(b_mag*np.sin(gamma)*v+c_mag*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma)*w)
            z=str(omega/(a_mag*b_mag*np.sin(gamma))*w)
            
            new.write("    ".join([elements[n],x,y,z])+'\n')
