import numpy as np
import re
import pandas as pd

##--- center of mass ---
##--- for unwrapped data only ---

def correct_image(atom_group, box):
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
    x =  atom_group['x'].values
    y =  atom_group['y'].values
    z =  atom_group['z'].values
    ix = np.zeros(len(x))
    iy = np.zeros(len(y))
    iz = np.zeros(len(z))
    # fix the wrong image flags, #gmx
    # Each pair of ions must not span further than half box.
    # They should always stay in the same pbc image.
    # The error comes from the gmx convert or image flag reset
    for ct in range(0, len(x)):
        if ct == 0:
            continue  # skip the first
        if x[ct]- x[0]> deltaX/2:
                ix[ct] -=1
        if x[ct]- x[0] < -deltaX/2:
                ix[ct] +=1
        
        if y[ct]- y[0]> deltaY/2:
                iy[ct] -=1
        if y[ct]- y[0] < -deltaY/2:
                iy[ct] +=1
        
        if z[ct]- z[0]> deltaZ/2:
                iz[ct] -=1
        if z[ct]- z[0] < -deltaZ/2:
                iz[ct] +=1
    #print(atom_group)
    return  ix.astype(int), iy.astype(int), iz.astype(int)

if __name__ == '__main__':
    # correct input flags
    from class_oneframe import oneframe
    from read_1_frame import read_1_lmp
    f = open('../Init.lammpstrj','r')
    time, Natom, box, cols, atomlist, position  = read_1_lmp(f)
    f.close()
    onef = oneframe()
    onef.load_snap(time, box, atomlist, cols)
    
    outfile = open('old.lammpstrj','w')
    onef.export_lmptrj(outfile, onef.L_atom)
    outfile.close()
    #
    #
    # Anion correction
    L_AN_mol = range(1,1+400)
    for mol in L_AN_mol:
        new_ix, new_iy, new_iz = correct_image( onef.L_atom.loc[onef.L_atom['mol'] == mol, :] , [[onef.xlo, onef.xhi],[onef.ylo, onef.yhi],[onef.zlo, onef.zhi]])
        onef.L_atom.loc[onef.L_atom['mol']==mol, 'ix'] = new_ix
        onef.L_atom.loc[onef.L_atom['mol']==mol, 'iy'] = new_iy
        onef.L_atom.loc[onef.L_atom['mol']==mol, 'iz'] = new_iz


    # cation correction
    # generate the list of all cations [[cation1], [cation2], ...] 

    # cation mols
    L_CT_mol = range(401,401+10)
    
    Deg_poly=40
    ions_per_mono = 1
    Natom_ION = 8
    all_ions_ids = []
    for i in L_CT_mol:
        # in each molecule
        iloc_Ter = [27,1132]
        CT_atom_ids = onef.L_atom[ onef.L_atom['mol'] == i ]['id'].values  # save in numpy array
        ##--- drop the terminal atoms and save them separately
        CT_atom_ids = np.delete(CT_atom_ids, iloc_Ter)
        #
        reduced_ids = CT_atom_ids - CT_atom_ids[0] # reset index starting from 0
        # create an empty list for all ids of ions
        ions_ids = []
        for j in range(0, len(CT_atom_ids) ): # must use length
            if (j%29<=5 and j%29>=1) or (j%29>=12 and j%29<=14):
                ions_ids.append( CT_atom_ids[j])
        # splite ions ids into single ions
        for k in range(0, Deg_poly*ions_per_mono):  # number of ions in one polymer chain
            all_ions_ids.append( ions_ids[k*Natom_ION : (k+1)*Natom_ION] )
        del CT_atom_ids
        del reduced_ids
        del ions_ids

    # now update image flags and molid
    N_mol = max(onef.L_atom['mol'])
    mol_counter = 1
    print(N_mol)
    for ion in all_ions_ids:
        # update image flags
        new_ix, new_iy, new_iz = correct_image(  onef.L_atom.loc[ ( onef.L_atom['id'].isin(ion)  ) ,: ],  [ [onef.xlo, onef.xhi],[onef.ylo, onef.yhi],[onef.zlo, onef.zhi] ]  )
        
        onef.L_atom.loc[onef.L_atom['id'].isin(ion), 'ix'] = new_ix
        onef.L_atom.loc[onef.L_atom['id'].isin(ion), 'iy'] = new_iy
        onef.L_atom.loc[onef.L_atom['id'].isin(ion), 'iz'] = new_iz
        
        # update mol id 
        # give the cationic groups new mol ids
        # for lammps to group them
        onef.L_atom.loc[onef.L_atom['id'].isin(ion), 'mol'] = N_mol+mol_counter

        mol_counter += 1
    outfile = open('new.lammpstrj','w')
    onef.export_lmptrj(outfile, onef.L_atom)
    outfile.close()
























    """


    atoms = pd.DataFrame([  [1, '1H2', 1, 0., 4., 0., 0., 0,0,0], \
                            [2, '2H', 1, 0., 1., 1., 2.,0,0,0], \
                            [3, 'H', 1, 0., 1., 2., 1.,0,0,0]  ] , \
                            index = [19,20,21], \
                            columns = ['atom', 'element', 'mol', 'q','x','y','z', 'ix', 'iy','iz'] \
                        )

    print(atoms)
    box = [[-0.1,4.9],[-0.1,3.9],[0.,5]]
    print(box)

    newix, newiy, newiz= correct_image(atoms.loc[atoms['atom']<3],box)
    atoms.loc[atoms['atom']<3,'ix'] = newix
    atoms.loc[atoms['atom']<3,'iy'] = newiy
    atoms.loc[atoms['atom']<3,'iz'] = newiz

    print(atoms)


    """
