import numpy as np
import sys
import os
#parentdir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..' ) )
#sys.path.append(parentdir)
from class_data import *
from class_oneframe import *
from util_ff import *
import matplotlib.pyplot as plt

d1 = data(savemem=0)
#d1.read_all_lmp( 'rigid_slab.dump')
d1.read_all_gsd('../rigid_slab.gsd') #print(d1.allframes[-1].L_atom)
for frame in d1.allframes:
    #frame.L_atom['res'] = frame.L_atom['type']
    #frame.L_atom['type'] = ['C'] * frame.Natom
    # nm to A
    frame.L_atom.loc[:, ['x', 'y', 'z' ] ] *= 10.
    frame.deltaX *= 10
    frame.deltaY *= 10
    frame.deltaZ *= 10
    frame.xhi *= 10
    frame.yhi *= 10
    frame.zhi *= 10
    frame.xlo *= 10
    frame.ylo *= 10
    frame.zlo *= 10
    # recover molid 
    Nchain = 100
    Nperchain = int(frame.Natom / Nchain )
    col_mol = []
    for i in range(0, Nchain):
        #col_mol += [ i+1 ] * Nperchain
        col_mol += [ i+1 for j in range(0,  Nperchain ) ]
    frame.L_atom['mol' ] = col_mol

    # drop rigid body center of mass
    mask = frame.L_atom['type'] != 'rbody'
    frame.L_atom = frame.L_atom[mask]
    frame.L_atom = frame.L_atom.reset_index( drop= True )
    frame.L_atom['id'] = frame.L_atom.index + 1

d1.allframes = d1.allframes[:: 100 ]
print("number of frames: ", len ( d1.allframes ) )

d1.gen_dcut( hps = '../stats_module.dat' )
H, x, y = d1.contact_map()
# plot
fig, ax = plt.subplots()
mesh = ax.pcolormesh( x, y, H )
cbar = fig.colorbar(mesh, ax=ax)
cbar.set_label("Counts")
ax.set_xlabel("Res 1")
ax.set_ylabel("Res 2")
plt.savefig('inter_contact_map.png', transparent=True, dpi=300)
#print(H, x, y)
#print(H.shape, x.shape, y.shape)
#print( d1.allframes[-1].L_atom )

"""
# add sigma (A)
ff_atomwise = read_ff_file( '../stats_module.dat' )
type_to_sig = ff_atomwise.copy()
for key in type_to_sig:
    type_to_sig[key] =  type_to_sig[key][2]
print( 'sigma for each type: ', type_to_sig )

for frame in d1.allframes:
    frame.L_atom['dcut'] =  frame.L_atom['type'].map(type_to_sig)

frame = d1.allframes[-1]
print(frame.L_atom)
#d_2 = frame.pairwise_d_2( frame.L_atom)
#print(d_2)
idx1, idx2, d_2, d_cut_2 = frame.inter_cont_list()

L_contact_pair = np.vstack( [ idx1, idx2] ).T
print(L_contact_pair)
print(L_contact_pair.shape)
print(d_2)
print(np.max(d_2) )
print(d_cut_2)
"""
