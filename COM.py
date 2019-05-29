import numpy as np
import re
import pandas as pd

##--- center of mass ---
##--- for unwrapped data only ---

def COM(atom_group, box, mapping_col = 'element', uxyz = False ):
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

    # initial image flags for the COM
    ix = 0
    iy = 0
    iz = 0

    # unwrap
    if uxyz:
        unwrap_x =  atom_group['ux']
        unwrap_y =  atom_group['uy']
        unwrap_z =  atom_group['uz']
    else:
        unwrap_x =  atom_group['x'] + deltaX * atom_group['ix']
        unwrap_y =  atom_group['y'] + deltaY * atom_group['iy']
        unwrap_z =  atom_group['z'] + deltaZ * atom_group['iz']
    # if no uxyz and no ixyz
    # this function would not work correctly
    # to speed up, here is no check

    # no need: group = atom_group.copy().reset_index(drop=True)
    # COM
    we_ux = np.dot(unwrap_x ,masses)/total_mass
    we_uy = np.dot(unwrap_y ,masses)/total_mass
    we_uz = np.dot(unwrap_z ,masses)/total_mass
    
    comx =  we_ux.sum()
    comy =  we_uy.sum()
    comz =  we_uz.sum()
    #print(comx, comy, comz, ix, iy, iz)

    # wrap back
    # check if out of boundary:
    if comx> xhi:
        ix =  int((comx - xlo)/deltaX)
        comx = comx - deltaX * ix
    elif comx < xhi:
        ix =  int((comx - xhi)/deltaX)
        comx = comx - deltaX * ix
    if comy> yhi:
        iy =  int((comy - ylo)/deltaY)
        comy = comy - deltaY * iy
    elif comy < yhi:
        iy =  int((comy - yhi)/deltaY)
        comy = comy - deltaY * iy
    if comz> zhi:
        iz =  int((comz - zlo)/deltaX)
        comz = comz - deltaZ * iz
    elif comz < zhi:
        iz =  int((comz - zhi)/deltaZ)
        comz = comz - deltaZ * iz





    return comx, comy, comz, ix, iy, iz

if __name__ == '__main__':
    atoms = pd.DataFrame([  [1, '1H2', 1, 4., 0., 0. , 0, -1, 1], \
                            [2, '2H', 1, 1., 1., 2., 0, -1, 1], \
                            [3, 'H', 1, 1., 2., 1., 0, -1, 1]  ] , \
                            index = [19,20,21], \
                            columns = ['atom', 'element', 'mol', 'x', 'y','z', 'ix', 'iy', 'iz'] \
                        )

    com_x, com_y, com_z, ix, iy ,iz= COM(atoms, [[-0.1,4.9],[-0.1,3.9],[0.,5]])
    print(com_x, com_y,com_z, ix, iy ,iz)

