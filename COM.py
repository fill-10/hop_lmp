import numpy as np
import re
import pandas as pd

##--- center of mass ---
##--- for unwrapped data only ---

def COM(atom_group, box, mapping_col = 'element'):
    mass_table = {'H':1.008, 'Li':6.94,'LI':6.94, 'B':10.81, 'C':12.011, 'N':14.007, 'O':15.999, 'F':18.998, 'Na':22.990, 'NA':22.990, 'Mg':4.305, 'MG':4.305,'Al':26.982, 'AL':26.982,'Si':28.085, 'SI':28.085, 'P':30.974, 'S':32.06, 'Cl':35.45, 'CL':35.45, 'K':39.098, 'Br':79.904, 'BR':79.904}
    
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
    if 'ux' in atom_group.columns:
        unwrap_x =  atom_group['ux']
    else:
        unwrap_x =  atom_group['x'] + deltaX * atom_group['ix']
    if 'uy' in atom_group.columns:
        unwrap_y =  atom_group['uy']
    else:
        unwrap_y =  atom_group['y'] + deltaY * atom_group['iy']
    if 'uz' in atom_group.columns:
        unwrap_z =  atom_group['uz']
    else:
        unwrap_z =  atom_group['z'] + deltaZ * atom_group['iz']
    
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
    atoms = pd.DataFrame([  [401, 'N', 1,  32.770,  39.280,   6.190], \
                            [402, 'C1', 1, 33.270,  38.880,   7.400], \
                            [403, 'C2', 1, 32.030,  38.310,   5.610],  \
                            [404, 'N1', 1, 31.850,  37.400,   6.640] , \
                            [405, 'C3', 1, 32.580,  37.740,   7.770] , \
                            [413, 'H',  1, 34.090,  39.300,   7.960]  , \
                            [414, 'H1', 1, 31.720,  38.160,   4.590] , \
                            [415, 'H2', 1, 32.470,  37.120,   8.660]   ],
                            index = [19,20,21,22,23,24,25,26], \
                            columns = ['atom', 'element', 'mol', 'ux', 'uy','uz'] \
                        )


    com_x, com_y, com_z, ix, iy ,iz= COM(atoms, [[0.,51.694],[0., 51.694],[0.,51.694]])
    print(com_x, com_y,com_z, ix, iy ,iz)

