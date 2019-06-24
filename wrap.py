import pandas as pd 
import numpy as np

def wrap(atom_group, box ):
    atom_group['ix'] = 0
    atom_group['iy'] = 0
    atom_group['iz'] = 0
    xlo = box[0][0]
    xhi = box[0][1]
    ylo = box[1][0]
    yhi = box[1][1]
    zlo = box[2][0]
    zhi = box[2][1]
    deltaX = abs(xhi -xlo)
    deltaY = abs(yhi -ylo)
    deltaZ = abs(zhi -zlo)
    for (idx, row) in atom_group.iterrows():
        if row['ux'] > xhi:
            atom_group.loc[idx, 'ix'] = int( (row['ux']-xlo) / deltaX)
        if row['ux'] < xlo:
            atom_group.loc[idx, 'ix'] = int( (row['ux']-xhi) / deltaX)
        if row['uy'] > yhi:
            atom_group.loc[idx, 'iy'] = int( (row['uy']-ylo) / deltaY)
        if row['uy'] < ylo:
            atom_group.loc[idx, 'iy'] = int( (row['uy']-yhi) / deltaY)
        if row['uz'] > zhi:
            atom_group.loc[idx, 'iz'] = int( (row['uz']-zlo) / deltaZ)
        if row['uz'] < zlo:
            atom_group.loc[idx, 'iz'] = int( (row['uz']-zhi) / deltaZ)

    atom_group['x'] = atom_group['ux'] - atom_group['ix']*self.deltaX
    atom_group['y'] = atom_group['uy'] - atom_group['iy']*self.deltaY
    atom_group['z'] = atom_group['uz'] - atom_group['iz']*self.deltaZ
    # all the above are pointer operatoins.
    # no need to return any value, actually.
    return atom_group['x'].values, atom_group['y'].values, atom_group['z'].values
