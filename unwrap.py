import pandas as pd 
import numpy as np

def unwrap(atom_group, box ):
    atom_group['ux'] = atom_group['x'] + atom_group['ix'] * (box[0][1]-box[0][0] )
    atom_group['uy'] = atom_group['y'] + atom_group['iy'] * (box[1][1]-box[1][0] )
    atom_group['uz'] = atom_group['z'] + atom_group['iz'] * (box[2][1]-box[2][0] )
    return atom_group['ux'].values, atom_group['uy'].values, atom_group['uz'].values
    # no need to return values actually,
    # atom_group is an object/pointer to the object

