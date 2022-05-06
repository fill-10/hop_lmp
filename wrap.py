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
        if row['xu'] > xhi:
            atom_group.loc[idx, 'ix'] = int( (row['xu']-xlo) / deltaX)
        if row['xu'] < xlo:
            atom_group.loc[idx, 'ix'] = int( (row['xu']-xhi) / deltaX)
        if row['yu'] > yhi:
            atom_group.loc[idx, 'iy'] = int( (row['yu']-ylo) / deltaY)
        if row['yu'] < ylo:
            atom_group.loc[idx, 'iy'] = int( (row['yu']-yhi) / deltaY)
        if row['zu'] > zhi:
            atom_group.loc[idx, 'iz'] = int( (row['zu']-zlo) / deltaZ)
        if row['zu'] < zlo:
            atom_group.loc[idx, 'iz'] = int( (row['zu']-zhi) / deltaZ)

    atom_group['x'] = atom_group['xu'] - atom_group['ix']*deltaX
    atom_group['y'] = atom_group['yu'] - atom_group['iy']*deltaY
    atom_group['z'] = atom_group['zu'] - atom_group['iz']*deltaZ
    # all the above are pointer operatoins.
    # no need to return any value, actually.
    return atom_group['x'].values, atom_group['y'].values, atom_group['z'].values
