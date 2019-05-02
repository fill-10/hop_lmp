import sys
import numpy as np
import pandas as pd
from read_1_frame import *
from COM import COM
from pbc_dist import *

##--- define a single frame class
class oneframe():
    def __init__(self):
        # box info
        # update by hand 
        self.dim = 3
        self.time = 0.0
        self.Natom = 0
        self.xhi = 0.1
        self.xlo = 0.0
        self.yhi = 0.1
        self.ylo = 0.0
        self.zhi = 0.1
        self.zlo = 0.0
        self.deltaX = 0.1
        self.deltaY = 0.1
        self.deltaZ = 0.1
        self.alpha  = 90.
        self.beta   = 90.
        self.gama   = 90.
        self.pbc = {'x': 1, 'y':1,'z':1}

        # atom info
        self.L_atom = []
        self.L_CT = []
        self.L_AN = []
    
    def update_pbc(self,box):
        self.xlo = box[0][0]
        self.xhi = box[0][1]
        self.ylo = box[1][0]
        self.yhi = box[1][1]
        self.deltaX = abs(self.xhi-self.xlo)
        self.deltaY = abs(self.yhi-self.ylo)
        try:
            self.zlo = box[2][0]
            self.zhi = box[2][1]
            self.deltaZ = abs(self.zhi-self.zlo)
        except:
            print('No Z direction data')
            self.dim = 2
    #
    def load_atoms(self, atoms,  cols=['id','mol','type','q', 'x', 'y','z','ix', 'iy', 'iz','element' ]):
        self.L_atom = pd.DataFrame(atoms, columns = cols)
        self.Natom  =  self.L_atom.shape[0]
    def load_snap(self, time, box, atoms, cols = ['id','mol','type','q', 'x', 'y','z','ix', 'iy', 'iz','element' ]):
        self.time = time
        self.update_pbc(box)
        self.load_atoms(atoms, cols)
    #
    def ion_gen(self, L_molid = range(401, 401+10), iontype='101', Deg_poly=40, ions_per_mono = 1, ilocL_Ter = [18,781], Natom_ION = 8, pick_cr_per_mon= "sel[  (  (sel.index%20 <=5) & (sel.index%20>=1) ) | ( (sel.index%20>=9 ) & (sel.index%20<=11)) ]"):
        ##--- create an empty list for ion candidates
        L_ion = []
        ##--- Loop over all molecules
        for i in L_molid:
            ##--- select one chain
            sel = self.atoms[ self.atoms['mol'] == i ].reset_index(drop=True)
            ##--- drop the terminal atoms
            sel = sel.drop(ilocL_Ter).reset_index(drop = True)
            
            ##--- select each ion based on the input criterion(-ria)
            ##--- the eval() module runs the string and returns the results
            sel1 = eval(pick_cr_per_mon)

            ##--- delete the verlet (temporary) atom dataframe
            del sel
            
            ##--- now, only the atoms in ions are selected in sel1
            ##--- split them into independent ions
            for j in range(0,Deg_poly*ions_per_mono):
                ##--- compute center of mass for each ion
                ion_x, ion_y, ion_z, ion_ix, ion_iy, ion_iz = COM( sel1.iloc[ j*Natom_ION : (j+1)*Natom_ION ]  )
                ##--- append this ion to the final list
                L_ion += [  [ (  (i-L_molid[0] )*Deg_poly + j  ) * ions_per_mono + 1, iontype, i, ion_x, ion_y, ion_z, ion_ix, ion_iy, ion_iz]  ]

            ##--- delete the temp dataframe sel1
            del sel1
        
        ##--- convert the result into pandas dataframe
        L_ion = pd.DataFrame(L_ion, columns = ['atom', 'type', 'mol', 'x', 'y', 'z', 'ix', 'iy', 'iz'])
        return L_ion

    def find_asso(self, L_im, L_mo, r_cut, do_stat = 1, diff_L=1):
        # l_im is the mother list, bring l_mo elements into l1
        if not diff_L:
            print('same list not developed...')
            return 1
       # prepare two lists for the mobile ions: associate atoms and associated molecules.
        L_mo_asso_atom = []
        L_mo_asso_mol  = []
        #
        # iter over rows (all the mobile ions)
        for (idx,row_mo) in L_mo.iterrows():
            coor1 = row_mo[['x','y','z']].values
            # prepare a list(column) of flags for immobile ions
            # if each of them associate with the current mobile ion
            L_im_if_asso_2_mo = []
            #
            # Loop over all immobile ions, 
            # check if each of them associate with the current mobile ion
            # if yes, mark 1
            # if no, mark 0
            for coor2 in L_im[['x', 'y','z']].values:
                if ifconn(coor1, coor2, [self.deltaX, self.deltaY, self.deltaZ], r_cut):
                    L_im_if_asso_2_mo.append(1)
                else:
                    L_im_if_asso_2_mo.append(0)
            #
            # attach this temporary column to the original dataframe,
            # this column is updated every time when looping over all mobile ions
            L_im['if_asso'] = L_im_if_asso_2_mo
            c_mo_asso_atom = L_im[L_im['if_asso'] !=0 ]['atom'].values#.tolist()
            c_mo_asso_mol  = np.unique(L_im[L_im['if_asso'] !=0 ]['mol'].values)
            L_mo_asso_atom.append( c_mo_asso_atom  )
            L_mo_asso_mol.append(  c_mo_asso_mol   )
            #print(c_mo_asso_atom)
            #print(c_mo_asso_mol)
        L_mo['asso_atom'] = L_mo_asso_atom
        L_mo['asso_mol']  = L_mo_asso_mol
        
        # stat
        if do_stat:
            N_asso_atom = np.array([]).astype(int)
            N_asso_mol  = np.array([]).astype(int)
            for (index,row) in L_mo.iterrows():
                N_asso_atom = np.append(N_asso_atom,  len(row['asso_atom']) )
                N_asso_mol  = np.append(N_asso_mol ,  len(row['asso_mol' ]) )
            hist_asso_im_atom = np.histogram( N_asso_atom , bins=np.arange(0, max( N_asso_atom) +2 ) )
            hist_asso_im_mol  = np.histogram( N_asso_mol  , bins=np.arange(0, max( N_asso_mol ) +2 ) )
            self.hist_asso_im_atom = hist_asso_im_atom
            self.hist_asso_im_mol  = hist_asso_im_mol
            return  N_asso_atom, N_asso_mol

    def delete_full(self):  # release mem
        del self.L_atom

    def export_pdb(self, filename, L_ion, Nth_frame=1):
        f = open(filename, 'w')
        f.write('REMARK\n'*3)
        f.write('CRYST1'+ '%9.3f'%self.deltaX +'%9.3f'%self.deltaY +'%9.3f'%self.deltaZ + '%7.2f'%self.alpha +'%7.2f'%self.beta+'%7.2f'%self.gama+' P 1           1\n')
        f.write('MODEL ' + ' '*4 + '%4d' %Nth_frame +'\n')
        if self.pandas :
            for (idx, row) in L_ion.iterrows():
                if len(row['type']) == 1:
                    Tatom = ' '+str(row['type'])+' '*2
                elif len(row['type']) == 2:
                    Tatom = str(row['type']) +' '*2
                elif len(row['type']) == 3:
                    Tatom = str(row['type']) +' '
                elif len(row['type']) == 4:
                    Tatom = str(row['type'])
                else:
                    Tatom = str(row['type'])

                f.write('ATOM  '+ '%5d' %row['atom'] + ' ' \
                + '%4s' %Tatom + ' ' + 'ION' + '  ' \
                + '%4d' %row['mol']+ '    ' \
                + '%8.3f' %row['x']  + '%8.3f' %row['y'] + '%8.3f' %row['z'] \
                + ' '*24 +'\n')
            f.write('TER\n')
        f.write('ENDMDL\n')
        f.close()
        return 0

    def export_lmptrj(self, filename, L_sel):
        f = open(filename, 'a')
        # write head:
        f.write('ITEM: TIMESTEP\n')
        f.write('%d' %self.time +'\n')
        f.write('ITEM: NUMBER OF ATOMS\n')
        N_sel = L_sel.shape[0]
        f.write('%d' %N_sel +'\n')
        f.write('ITEM: BOX BOUNDS'+' pp'*3 +'\n') # only support pbc now
        #f.write()
        #f.write()
        #f.write()
        f.close()
            
##-----------------------------

if __name__ == '__main__':
    import time as timer
    ##--- read-in the whole file ---
    filename = '../first5.lammpstrj'
    f = open(filename, 'r')

    # read all
    ##--- extract frames ---
    read_pos = 0
    atom_dict = {}
    start = timer.perf_counter()
    atomlist, Natom, box, time, cols, position  = read_snap(f, read_pos)
    onef = oneframe()
    onef.load_snap(time, box, atomlist, cols)
    print(onef.L_atom.shape)
    onef.export_lmptrj('zz.txt', onef.L_atom)
    """
    al1, al2, read_pos, Natom1, Natom2, box, time = sel2l_bytype_1f(raw, read_pos, '*N*', '*S*')
    onef = oneframe(box=box, time=time, pandas=1) #initialize object
    onef.load_2_lists(al1, al2)
    #print(onef.atoms, '\n' ,  onef.Natom, '\n', onef.time,'\n', onef.box,'\n' ,onef.atom_list1[:5], '\n', onef.Natom1,'\n', onef.atom_list2[:5],'\n', onef.Natom2, '\n')
    """
    
    # select
    onef.L_CT = onef.ion_gen(range(401, 401+10),'CT', 40, 1, [18,781], 8, "sel[  (  (sel.index%20 <=5) & (sel.index%20>=1) ) | ( (sel.index%20>=9 ) & (sel.index%20<=11)) ]")
    #onef.L_AN = onef.ion_gen(range(1, 1+400), 'AN', 1, 1, [], 15, "sel[:]" )
    anion_args = (range(1, 1+400), 'AN', 1, 1, [], 15, 'sel[:]' )
    onef.L_AN = onef.ion_gen(*anion_args)
    del onef.atoms
    
    ##--- export ions to pdb ---
    onef.L_CT['atom'] += 2000
    onef.export_pdb('t'+str(onef.time)+'.pdb', pd.concat([onef.L_AN, onef.L_CT], ignore_index=True ) )
    ##--------------------------
    
    #hist_atom_NCT, hist_mol_NCT =
    onef.find_asso(onef.L_CT,onef.L_AN, 7.8)

    
    print(onef.hist_asso_im_atom, onef.hist_asso_im_mol)

    end   = timer.perf_counter()
    print(end-start)

    
    #print(onef.L_CT)
    #print(onef.L_AN)
    #print(onef.L_AN.columns.get_loc('type'))

