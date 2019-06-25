import sys
import numpy as np
import pandas as pd
from read_1_frame import *
from read_1_fix import read_1_fix
from COM import COM
from NGP import NGP
from VANHOVE import *
from pbc_dist import *
from unwrap import unwrap

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
    def ion_gen(self, L_molid = range(401, 401+10), iontype='101', Deg_poly=40, ions_per_mono = 1, ilocL_Ter = [18,781], Natom_ION = 8, pick_cr_per_mon= "sel[  (  (sel.index%20 <=5) & (sel.index%20>=1) ) | ( (sel.index%20>=9 ) & (sel.index%20<=11)) ]" , COM_map_col= 'element'):
        ##--- create an empty list for ion candidates
        L_ion = []
        ##--- Loop over all molecules
        for i in L_molid:
            ##--- select one chain
            sel = self.L_atom[ self.L_atom['mol'] == i ].reset_index(drop=True)
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
                ## COM in this version always returns x,y,z and ix,iy,iz
                ion_x, ion_y, ion_z, ion_ix, ion_iy, ion_iz = COM( sel1.iloc[ j*Natom_ION : (j+1)*Natom_ION ], [[self.xlo, self.xhi],[self.ylo, self.yhi],[self.zlo, self.zhi]], COM_map_col )
                ##--- append this ion to the final list
                L_ion += [  [ (  (i-L_molid[0] )*Deg_poly + j  ) * ions_per_mono + 1, i, iontype, ion_x, ion_y, ion_z, ion_ix, ion_iy, ion_iz]  ]
            ##--- delete the temp dataframe sel1
            del sel1
        
        ##--- convert the result into pandas dataframe
        L_ion = pd.DataFrame(L_ion, columns = ['id', 'mol', 'type', 'x', 'y', 'z', 'ix', 'iy', 'iz'])
        return L_ion
    
    def read_ion(self, atoms, iontype, ionspermol=1, dropuxyz = False):
        L_ion = pd.DataFrame(atoms, columns = ['id','ux','uy','uz']) 
        # The columns are from lammps 'fix' output.
        L_ion.loc[:,'ix'] = 0 
        L_ion.loc[:,'iy'] = 0 
        L_ion.loc[:,'iz'] = 0
        for (idx, row) in L_ion.iterrows():
            if row['ux'] > self.xhi:
                L_ion.loc[idx, 'ix'] = int( (row['ux']-self.xlo) / self.deltaX)
            if row['ux'] < self.xlo:
                L_ion.loc[idx, 'ix'] = int( (row['ux']-self.xhi) / self.deltaX)
            if row['uy'] > self.yhi:
                L_ion.loc[idx, 'iy'] = int( (row['uy']-self.ylo) / self.deltaY)
            if row['uy'] < self.ylo:
                L_ion.loc[idx, 'iy'] = int( (row['uy']-self.yhi) / self.deltaY)
            if row['uz'] > self.zhi:
                L_ion.loc[idx, 'iz'] = int( (row['uz']-self.zlo) / self.deltaZ)
            if row['uz'] < self.zlo:
                L_ion.loc[idx, 'iz'] = int( (row['uz']-self.zhi) / self.deltaZ)
        L_ion['type'] = iontype
        # write mol id
        col_mol = []
        for imol in range(1, int(L_ion.shape[0]/ionspermol)+1):
            for j in range(0, ionspermol):
                col_mol += [imol]
        L_ion['mol'] = col_mol
        ##
        L_ion['x'] = L_ion['ux'] - L_ion['ix']*self.deltaX
        L_ion['y'] = L_ion['uy'] - L_ion['iy']*self.deltaY
        L_ion['z'] = L_ion['uz'] - L_ion['iz']*self.deltaZ
        # L_ion.reset_index(drop=True)
        # L_ion.loc['id'] = L_ion.index+1
        if dropuxyz:
            return L_ion.drop(columns= ['ux', 'uy' ,'uz'])
        else:
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
            c_mo_asso_atom = L_im[L_im['if_asso'] !=0 ]['id'].values#.tolist()
            c_mo_asso_mol  = np.unique(L_im[L_im['if_asso'] !=0 ]['mol'].values)
            L_mo_asso_atom.append( c_mo_asso_atom  )
            L_mo_asso_mol.append(  c_mo_asso_mol   )
            #print(c_mo_asso_atom)
            #print(c_mo_asso_mol)
        L_mo['asso_atom'] = L_mo_asso_atom
        L_mo['asso_mol']  = L_mo_asso_mol
        
        # stat
        if do_stat:
            # create two arrays for number of asso. atoms and number of asso. chains for all the mobile ions
            N_asso_atom = np.array([]).astype(int)
            N_asso_mol  = np.array([]).astype(int)
            for (index,row) in L_mo.iterrows():
                N_asso_atom = np.append(N_asso_atom,  len(row['asso_atom']) ) # append the number of asso. atom for each mobile ion
                N_asso_mol  = np.append(N_asso_mol ,  len(row['asso_mol' ]) ) # append the number of asso. chain for each mobile ion
            hist_asso_im_atom = np.histogram( N_asso_atom , bins=np.arange(0, max( N_asso_atom) +2 ) ) # histo
            hist_asso_im_mol  = np.histogram( N_asso_mol  , bins=np.arange(0, max( N_asso_mol ) +2 ) ) # histo
            self.hist_asso_im_atom = hist_asso_im_atom
            self.hist_asso_im_mol  = hist_asso_im_mol
            return  N_asso_atom, N_asso_mol
        else:
            return 0

    def delete_full(self):  # release mem
        del self.L_atom

    def export_pdb(self, filename, L_ion, Nth_frame=1):
        f = open(filename, 'a')
        f.write('REMARK\n'*3)
        f.write('CRYST1'+ '%9.3f'%self.deltaX +'%9.3f'%self.deltaY +'%9.3f'%self.deltaZ + '%7.2f'%self.alpha +'%7.2f'%self.beta+'%7.2f'%self.gama+' P 1           1\n')
        f.write('MODEL ' + ' '*4 + '%4d' %Nth_frame +'\n')
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

    def export_lmptrj(self, f, L_sel):  # f is the file object, i.e., pointer
        # write head:
        f.write('ITEM: TIMESTEP\n')
        f.write('%d' %self.time +'\n')
        f.write('ITEM: NUMBER OF ATOMS\n')
        N_sel = L_sel.shape[0]
        f.write('%d' %N_sel +'\n')
        f.write('ITEM: BOX BOUNDS'+' pp'*3 +'\n') # only support pbc now
        f.write('{:.5f}    {:.5f}\n'.format(self.xlo,self.xhi) )
        f.write('{:.5f}    {:.5f}\n'.format(self.ylo,self.yhi) )
        f.write('{:.5f}    {:.5f}\n'.format(self.zlo,self.zhi) )
        f.write('ITEM: ATOMS '+ ' '.join(L_sel.columns.values.astype(str)) +'\n'  )
        # write atom body
        L_sel.to_csv(f, sep=' ', float_format='%.5f', header=False, index=False)
    
    def hoppingtype_AN(self, prev): #current and prev snap
        asso_type = []
        for ind, row in self.L_AN.iterrows():
            #print(ind, row['atom'],type(row['asso_atom']))
            if len(row['asso_atom']) == 0 or len(prev.L_AN.loc[ind,:]['asso_atom'] ) == 0:
                asso_type += [3]
            elif np.array_equal(  row['asso_atom'] , prev.L_AN.loc[ind,:]['asso_atom']  ):
                asso_type += [4]
            elif np.array_equal(  row['asso_mol']  , prev.L_AN.loc[ind,:]['asso_mol']   ):
                asso_type += [1]
            else:
                asso_type += [2]

        hist_hop_type = np.histogram(asso_type, bins=[1,2,3,4,5])

        return hist_hop_type # not normalized
    
    def unwrap(self, L_ion): # L_ion is a pointer
        L_ion['ux'], L_ion['uy'], L_ion['uz'] = unwrap(L_ion, [[self.xlo, self.xhi], [self.ylo, self.yhi], [self.zlo, self.zhi]])

    def nongauss(self,L_mobile_ions,ref): # non gaussian parameter data point
        return NGP(L_mobile_ions['ux'], L_mobile_ions['uy'], L_mobile_ions['uz'], ref['ux'], ref['uy'], ref['uz'])
    
    def vanhove_s(self, L_mobile_ions, ref, maxdist = 25.0, accuracy = 0.1): # VH_s data point
        return  VANHOVE_S( L_mobile_ions['ux'], L_mobile_ions['uy'], L_mobile_ions['uz'], ref['ux'], ref['uy'], ref['uz'], maxdist, accuracy ) 
    
    def vanhove_d(self, L_mobile_ions, ref, maxdist = 25.0, accuracy = 0.1): #VH_d data piont
        #return  VANHOVE_D( L_mobile_ions['ux'], L_mobile_ions['uy'], L_mobile_ions['uz'], ref['ux'], ref['uy'], ref['uz'], maxdist, accuracy ) 
        return  VANHOVE_D( L_mobile_ions['x'], L_mobile_ions['y'], L_mobile_ions['z'], ref['x'], ref['y'], ref['z'], maxdist, accuracy ) 
    
    def findfast(self, L_mobile_ions, ref, r_star=6.0):
        L_mobile_ions['fast'] = 0
        r_star_squared = r_star*r_star
        L_ds_squared = (L_mobile_ions['ux']-ref['ux'])**2 +(L_mobile_ions['uy']-ref['uy'])**2 + (L_mobile_ions['uz']-ref['uz'])**2
        for (idx, val) in L_ds_squared.iteritems():
            if val > r_star_squared:
                L_mobile_ions.loc[idx, 'fast'] = 1
        return 0 

    def findstring(self, L_mobile_ions, ref, cutoff=2.5):
        if 'fast' in L_mobile_ions.columns:
            L_mobile_ions['string'] = -1
            current_string = 1
            for (idx, row) in L_mobile_ions[L_mobile_ions['fast']>0].iterrows():
                L_mobile_ions[idx, 'string'] = 0
                # get the coor of mobile ion:
                coor_mob = np.array([row['x'], row['y'], row['z'])
                for (idx2, row2) in ref[ref.index>idx].iterrows():
                    if L_mobile_ions.loc[idx2, 'fast'] <=0:
                        continue # pass the slow ions
                    if L_mobile_ions.loc[idx2,'string'] == row['string']:
                        continue # pass previously combined
                    else:
                        # get the coor of ref ion:
                        coor_ref = np.array([ row2['x'], row2['y'], row['z'] )
                        if ifconn(coor_mob, coor_ref, np.array([self.deltaX, self.deltaY,self.deltaZ]), cutoff):
                            ref_string_id = L_mobile_ions.loc[idx2,'string']
                            if row['string'] <=0: # mob does not in string
                                if ref_string_id >0: # ref in string
                                    L_mobile_ions.loc[idx,'string'] = ref_string_id
                                else: # it's a new string
                                    L_mobile_ions.loc[idx, 'string'] = current_string
                                    L_mobile_ions.loc[idx2,'string'] = current_string
                                    current_string +=1
                            elif ref_string_id <=0: # mob is in string but ref not in string
                                L_mobile_ions.loc[idx2,'string'] = row['string']
                            else: # both two are in strings, need to combine two strings
                                larger_string_id = max(row.loc[idx2, 'string'], row['string'])
                                smaller_string_id = min(row.loc[idx2, 'string'], row['string'])
                                L_mobile.loc[ L_mobile_ions['string']==larger_string_id, 'string'] = smaller_string_id
            string_histo = np.histogram(L_mobile_ions[L_mobile_ions['fast']>0] , bins=np.arange(0, L_mobile_ions.shape[0] +2)
            return string_histo  # not normalized
        else:
            print(need to find 'fast' atoms!)
            return -1


##-----------------------------

if __name__ == '__main__':
    import time as timer
    ##--- read-in the whole file ---
    filename = '../400K_corrected_lmp_B/VImC2_every100ps_0-50ns.lammpstrj'
    ANfile = '../400K_corrected_lmp_B/com_AN_every10ps_0-10ns.dat'
    CTfile = '../400K_corrected_lmp_B/com_CT_every10ps_0-10ns.dat'
    f = open(filename, 'r')
    fa = open(ANfile, 'r')
    fc = open(CTfile, 'r')
    ##--- timer starts ---
    start = timer.perf_counter()

    ##--- initialize a data object ---
    onef = oneframe()
    ###--- write box dim manually ---
    onef.update_pbc([[0.0,55.8983],[0.02,55.9183],[0.01,55.9083]])
    ###--- read ions from lammps 'fix' --
    anions, Nanion, time, pos = read_1_fix(fa)
    cations, Ncations, time, pos = read_1_fix(fc)
    ###--- organize ions in to object instances ---
    onef.L_CT = onef.read_ion(cations, 2, 40)
    onef.L_AN = onef.read_ion(anions, 1, 1)
    onef.L_CT.loc[:,'mol'] += max(onef.L_AN['mol'])
    
    ##--- read second frame ---
    secf = oneframe()
    secf.update_pbc([[0.0,55.8983],[0.02,55.9183],[0.01,55.9083]])
    anions, Nanion, time, pos = read_1_fix(fa)
    cations, Ncations, time, pos = read_1_fix(fc)
    secf.L_CT = onef.read_ion(cations, 2, 40)
    secf.L_AN = onef.read_ion(anions, 1, 1)
    secf.L_CT.loc[:,'mol'] += max(secf.L_AN['mol'])
    
    ###--- read first frame again from lammps 'trj' ---
    thirdf = oneframe()
    time, Natom, box, cols, atoms, position  = read_1_lmp(f)
    thirdf.load_snap(time, box, atoms, cols)
    # select ions
    anion_args = (range(1, 1+400), 1, 1, 1, [], 15, 'sel[:]' )
    thirdf.L_AN = thirdf.ion_gen(*anion_args)
    thirdf.L_CT = thirdf.ion_gen(range(411, 411+400),2, 1, 1, [], 8, "sel[:]")
    # release mem
    del thirdf.L_atom
    # correct id
    thirdf.L_CT['id'] += 400
    # correct mol
    thirdf.L_CT['mol'] = 0
    NDeg_Poly_CT = 40
    NIon_permon = 1
    Nmol_CT = 10
    for j0 in range(0, Nmol_CT):
        for j1 in range(j0*NDeg_Poly_CT*NIon_permon, (j0+1)*NDeg_Poly_CT*NIon_permon):
            thirdf.L_CT.loc[j1, 'mol'] = j0+1+thirdf.L_AN.shape[0]
    
    ##--- export ions to lammpstrj ------
    onef.export_lmptrj(open('one_'+str(onef.time)+'.lammpstrj','w'), pd.concat([onef.L_AN, onef.L_CT], ignore_index=True) )
    thirdf.export_lmptrj(open('third_'+str(onef.time)+'.lammpstrj','w' ), pd.concat([thirdf.L_AN, thirdf.L_CT], ignore_index=True) )
    
    ##--- find associating ions and mols for 2 consequent frames ---
    hist_atom_NCT, hist_mol_NCT =onef.find_asso(onef.L_CT,onef.L_AN, 7.8)
    hist_atom_NCT, hist_mol_NCT =secf.find_asso(secf.L_CT,secf.L_AN, 7.8)
    print('acco_atoms: ', onef.hist_asso_im_atom, '\nacco_mols: ', onef.hist_asso_im_mol)
    
    ##--- calc hopping types for the 2 frames ---
    print('one2second hop:', secf.hoppingtype_AN(onef)[0] )
    
    ###--- non gauss ---
    print('non gauss data point,secf- onef: ',secf.nongauss(secf.L_AN, onef.L_AN))
    ###--- van hove ---
    print('van hove self, secf - onef: ', secf.vanhove_s(secf.L_AN, onef.L_AN, 25.0, 0.1))
    print('van hove distinct, secf - onef: ', secf.vanhove_d(secf.L_AN, onef.L_AN, 25.0, 0.1))
    ##--- timer ends ---
    end = timer.perf_counter()
    print('time used: ', end-start)
