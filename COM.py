import numpy as np
import re
import pandas as pd

##--- center of mass ---
##--- for unwrapped data only ---

def COM(atom_group, box, mapping_col = 'element'):
    mass_table = {'H':1.008, 'Li':6.94, 'B':10.81, 'C':12.011, 'N':14.007, 'O':15.999, 'F':18.998, 'Na':22.990, 'Mg':4.305, 'Al':26.982, 'Si':28.085, 'P':30.974, 'S':32.06, 'Cl':35.45, 'K':39.098, 'Br':79.904 }
    
    # box info
    xlo = box[0][0]
    xhi = box[0][1]
    ylo = box[1][0]
    yhi = box[1][1]
    zlo = box[2][0]
    zhi = box[2][1]
    deltaX = xhi-xlo
    deltaY = yhi-ylo
    deltaZ = zhi-zlo
    #print(deltaX, deltaY, deltaZ)
    masses = atom_group[mapping_col].str.replace('\d+', '').map(mass_table)
    total_mass = masses.sum()

    # image flags for the COM
    ix = 0
    iy = 0
    iz = 0

    # save in local var
    group = atom_group.copy().reset_index(drop = True)

    # read in image flags
    if 'ix' in atom_group.columns:
        ix = atom_group.iloc[0]['ix'] # first atom as reference
    else:
        group.loc[:,'ix'] = 0 # if image flag not given, set them to 0 for all atoms
    if 'iy' in atom_group.columns:
        iy = atom_group.iloc[0]['iy']
    else:
        group.loc[:,'iy'] = 0
    if 'iy' in atom_group.columns:
        iz = atom_group.iloc[0]['iz']
    else:
        group.loc[:,'iz'] = 0

    """
    # fix the wrong image flags, #gmx
    # this is not good , have been rewritten below.
    for (idx, row) in group.iterrows():
        if idx ==0:
            continue  # skip the first
        if row['x'] - group.iloc[0]['x'] > deltaX/2:
            if row[idx]['ix'] == group.iloc[0]['ix']:
                group.loc[idx,'ix'] -=1 # must use loc[ , ], or cannot access original data
        elif row['x'] - group.iloc[0]['x'] < -deltaX/2:
            if row['ix'] == group.iloc[0]['ix']:
                group.loc[idx,'ix'] +=1
        if row['y'] - group.iloc[0]['y'] > deltaY/2:
            if row['iy'] == group.iloc[0]['iy']:
                group.loc[idx,'iy'] -=1
        elif row['y'] - group.iloc[0]['y'] < -deltaY/2:
            if row[idx]['iy'] == group.iloc[0]['iy']:
                group.loc[idx,'iy'] +=1
        if row['z'] - group.iloc[0]['z'] > deltaZ/2:
            if row['iz'] == group.iloc[0]['iz']:
                group.loc[idx,'iz'] -=1
        elif row['z'] - group.iloc[0]['z'] < -deltaZ/2:
            if row['iz'] == group.iloc[0]['iz']:
                group.loc[idx,'iz'] += 1
    print(group)
    """
    # unwrap
    if not 'ux' in atom_group.columns:
        group['ux'] = group['x'] + deltaX * group['ix']
    if not 'uy' in atom_group.columns:
        group['uy'] = group['y'] + deltaY * group['iy']
    if not 'uz' in atom_group.columns:
        group['uz'] = group['z'] + deltaZ * group['iz']
    #print(group)
    
    # fix the wrong image flags, #gmx
    # each pair of ions must not span further than half box
    # they should always stay in the same pbc image
    # the error comes from the gmx convert or image flag reset
    for (idx, row) in group.iterrows():
        if idx ==0:
            continue  # skip the first
        while group.iloc[idx]['ux'] - group.iloc[0]['ux'] > deltaX/2:
                group.loc[idx,'ux'] -= deltaX
                group.loc[idx,'ix'] -=1 # must use loc[ , ], or cannot access original data
        while group.iloc[idx]['ux'] - group.iloc[0]['ux'] < -deltaX/2:
                group.loc[idx,'ux'] += deltaX
                group.loc[idx,'ix'] +=1
        while group.iloc[idx]['uy'] - group.iloc[0]['uy'] > deltaY/2:
                group.loc[idx,'uy'] -= deltaY
                group.loc[idx,'iy'] -=1 # must use loc[ , ], or cannot access original data
        while group.iloc[idx]['uy'] - group.iloc[0]['uy'] < -deltaY/2:
                group.loc[idx,'uy'] += deltaY
                group.loc[idx,'iy'] +=1
        while group.iloc[idx]['uz'] - group.iloc[0]['uz'] > deltaZ/2:
                group.loc[idx,'uz'] -= deltaZ
                group.loc[idx,'iz'] -=1 # must use loc[ , ], or cannot access original data
        while group.iloc[idx]['uz'] - group.iloc[0]['uz'] < -deltaZ/2:
                group.loc[idx,'uz'] += deltaZ
                group.loc[idx,'iz'] +=1
    #print(group)

    # COM
    we_ux = np.dot(group['ux'] ,masses)/total_mass
    we_uy = np.dot(group['uy'] ,masses)/total_mass
    we_uz = np.dot(group['uz'] ,masses)/total_mass

    # wrap back
    comx =  we_ux.sum()
    comy =  we_uy.sum()
    comz =  we_uz.sum()
    #print(comx, comy, comz, ix, iy, iz)

    comx =  we_ux.sum() - ix*(xhi-xlo)
    comy =  we_uy.sum() - iy*(yhi-ylo)
    comz =  we_uz.sum() - iz*(zhi-zlo)
    #print(comx, comy, comz, ix, iy, iz)

    # check if out of boundary:
    if comx > xhi:
        comx -= (xhi -xlo)
        ix += 1
    elif comx < xlo:
        comx += (xhi-xlo)
        comx -= 1
    if comy > yhi:
        comy -= (yhi -ylo)
        iy += 1
    elif comy < ylo:
        comy += (yhi-ylo)
        iy -= 1
    if comz > zhi:
        comz -= (zhi -zlo)
        iz += 1
    elif comz < zlo:
        comz += (zhi-zlo)
        iz -= 1
    return comx, comy, comz, ix, iy, iz

if __name__ == '__main__':
    atoms = pd.DataFrame([  [1, '1H2', 1, 1., 0., 0. , 1, -1, 1], \
                            [2, '2H', 1, 1., 1., 2., 1, -1, 1], \
                            [3, 'H', 1, 1., 2., 1., 1, -1, 1]  ] , \
                            index = [19,20,21], \
                            columns = ['atom', 'element', 'mol', 'x', 'y','z', 'ix', 'iy', 'iz'] \
                        )

    com_x, com_y, com_z, ix, iy ,iz= COM(atoms, [[-0.1,4.9],[-0.1,3.9],[0.,5]])
    print(com_x, com_y,com_z, ix, iy ,iz)

