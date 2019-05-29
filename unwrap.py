import pandas as pd 
import numpy as np

def unwrap(atom_group, box ):
    ux = atom_group['x'] + atom_group['ix'] * (box[0][1]-box[0][0] )
    uy = atom_group['y'] + atom_group['iy'] * (box[1][1]-box[1][0] )
    uz = atom_group['z'] + atom_group['iz'] * (box[2][1]-box[2][0] )
    return ux.values, uy.values, uz.values

