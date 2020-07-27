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
        shell_distances = [ d for d in r_ij if d >= cutoff1 and d <= cutoff2 ]
 
        #print('Original DataFrame')
        #print(sorted_indices)
        shell_indices = shell_indices.groupby('col1')['col2'].apply(list)
        #print('Grouped by DataFrame')
        #print(newdf)
        return [shell_indices, shell_distances]
my_traj = read('CONTCAR_big',index=':')
my_poly = Polyhedra(my_traj)
[i,r]=my_poly.get_shell(0,'In','O',0.0,2.4)
pd.DataFrame( i ).to_csv( "first_shell.dat",sep="\t",header=False,index=True )

#def get_shell(el1, el2, cutoff1, cutoff2):
## first, for each In atom, find all O atoms within 2.36 Angstroms
#for n in range(1):
#    first_shell=[]
#    first_shell_distances=[]
#    in_shells=[]
#    print("Step: "+str(n+1))
#    start = n * ( n_atoms + 1 ) + 8
#    end = n * ( n_atoms + 1 ) + 7 + n_atoms
#    coords=[[float(i) for i in d.split()] for d in data[start:end+1]]
#    for i in range(n_per[0]):
#        coord_1=np.array(coords[i])
#        distances=[]
#        for j in range(n_per[0],n_per[0]+n_per[1]):
#            coord_2=np.array(coords[j])
#            pdiff=np.abs(coord_1-coord_2)
#            mod=((np.dot(LCC_inv,pdiff)+0.5)%1)-0.5
#            pdiff_corr=np.matmul([a,b,c],mod)
#            distances.append(np.linalg.norm(pdiff_corr))
#        sorted_indices=np.argsort(np.array(distances))
#        distances.sort()
#        first_shell.append([ sorted_indices[k]+256 for k in range(len(distances)) if distances[k] <= 2.36 ])
#        first_shell_distances.append([ d for d in distances if d <= 2.36 ])
#
#        distances=[]
#        j_indices=[] # to ensure that I map from index of distance to index of In atom
#        for j in range(n_per[0]): 
#            if i != j:
#                j_indices.append(j)
#                coord_2=np.array(coords[j])
#                pdiff=np.abs(coord_1-coord_2)
#                mod=((np.dot(LCC_inv,pdiff)+0.5)%1)-0.5
#                pdiff_corr=np.matmul([a,b,c],mod)
#                distances.append(np.linalg.norm(pdiff_corr))
#        sorted_indices=np.argsort(np.array(distances))
#        distances.sort() 
#        in_shells.append([ j_indices[sorted_indices[k]] for k in range(len(distances)) if distances[k] <= 4.0 ])
#
#    # now that we have the first shell of each In atom (i.e., the O atoms
#    # coordinating each In), we need to identify the different types of polyhedra.
#    # we can do this by determining the number of In neighbors that share 1, 2, or 3
#    # O atoms with the central In atom:
#    # 1 --> corner-sharing polyhedra
#    # 2 --> edge-sharing polyhedra
#    # 3 --> face-sharing polyhedra
#   
#    print(np.mean(flatten(first_shell_distances))) 
#    np.savetxt("first_shell_distances.txt",first_shell_distances,fmt='%s')
#    corner = np.zeros( n_per[0] )
#    edge = np.zeros( n_per[0] )
#    face = np.zeros( n_per[0] )
#    corner_indices = []
#    edge_indices = [] 
#    face_indices = [] 
#
#    # loop through each In atom
#    for i in range(n_per[0]):
#        O_atoms = first_shell[i]
#        corner_temp = [i]
#        edge_temp = [i]
#        for j in in_shells[i]:
#            shared_O_atoms = [ v for v in O_atoms if v in first_shell[j] ]
#            if len(shared_O_atoms) == 1:
#                corner[i] += 1
#                corner_temp.append(j)
#            if len(shared_O_atoms) == 2:
#                edge[i] += 1
#                edge_temp.append(j)
#            if len(shared_O_atoms) == 3:
#                face[i] += 1
#
#        corner_indices.append(corner_temp)
#        edge_indices.append(edge_temp)
#    for c in corner_indices:
#        print(c)
#    #edge_max=int(max(edge))
#    #edge_min=int(min(edge))
#    #corner_max=int(max(corner))
#    #corner_min=int(min(corner))
#    #print(corner_min)
#    #hist=np.zeros( (edge_max-edge_min+1, corner_max-corner_min+1) )
#    #for i in range(n_per[0]):
#    #    c = int(corner[i])
#    #    e = int(edge[i])
#    #    hist[e-edge_min,c-corner_min] += 1
#    #fig,ax=plt.subplots()
#    #ax.imshow(hist,cmap="Blues")
#    #ax.set_xlabel("Number of Corner Polyhedra")
#    #ax.set_ylabel("Number of Edge Polyhedra")
#    #ax.set_yticks(np.arange(0.5,edge_max-edge_min+1,1.0), minor='True')
#    #ax.set_xticks(np.arange(0.5,corner_max-corner_min+1,1.0), minor='True')
#    #ax.set_yticks(np.arange(0,edge_max-edge_min+1,1.0))
#    #ax.set_xticks(np.arange(0,corner_max-corner_min+1,1.0))
#    #ax.set_yticklabels(range(edge_min,edge_max+1))
#    #ax.set_xticklabels(range(corner_min,corner_max+1))
#    #ax.yaxis.grid(True, linestyle='--',linewidth=1.5, which='minor')
#    #ax.xaxis.grid(True, linestyle='--',linewidth=1.5, which='minor')
#    #fig.savefig("polyhedra_heatmap_a_In2O3.png",dpi=600)
#
#    #fig,ax=plt.subplots()
#    #ax.hist(corner,9,color='mediumblue',edgecolor='dimgray',alpha=0.5,label='Corner-Sharing')
#    #ax.hist(edge,6,color='limegreen',edgecolor='dimgray',alpha=0.5,label='Edge-Sharing')
#    #ax.set_xlabel("Number of Polyhedra")
#    #ax.set_ylabel("Frequency")
#    #ax.legend()
#    #fig.tight_layout()
#    #fig.savefig("polyhedra_a_In2O3.png",dpi=600)
