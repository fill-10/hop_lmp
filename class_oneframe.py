import sys
import numpy as np
import pandas as pd
from read_1_frame import *
from COM import COM
from NGP import NGP
from MSD import MSD
from VANHOVE import *
from pbc_dist import *
from unwrap import unwrap
from wrap import wrap
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
    def load_atoms(self, atoms, cols):
        self.L_atom = pd.DataFrame(atoms, columns = cols)
        self.Natom  =  self.L_atom.shape[0]
    def load_snap(self, time, box, atoms, cols = ['id','mol','type', 'x', 'y','z','ix', 'iy', 'iz','element' ]):
        self.time = time
        self.update_pbc(box)
        self.load_atoms(atoms, cols)
    #
    def unwrap(self, L_ion): # L_ion is a pointer to an object
        unwrap(L_ion, [[self.xlo, self.xhi], [self.ylo, self.yhi], [self.zlo, self.zhi]])
    def wrap(self, L_ion):
        wrap(L_ion, [[self.xlo, self.xhi], [self.ylo, self.yhi], [self.zlo, self.zhi]])
    #
    def ion_gen(self, L_molid = range(401, 401+10), iontype='101', Deg_poly=40, ions_per_mono = 1, ilocL_Ter = [18,781], Natom_ION = 8,\
                pick_cr_per_mon= "sel[  (  (sel.index%20 <=5) & (sel.index%20>=1) ) | ( (sel.index%20>=9 ) & (sel.index%20<=11)) ]" , COM_map_col='type'):
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
            ##--- delete the verlet (temporary) atom DataFrame
            del sel
            ##--- now, only the atoms in ions are selected in sel1
            ##--- split them into independent ions
            for j in range(0,Deg_poly*ions_per_mono):
                ##--- compute center of mass for each ion
                ion_ux, ion_uy, ion_uz= COM( sel1.iloc[ j*Natom_ION : (j+1)*Natom_ION ], COM_map_col )
                ##--- append this ion to the final list
                L_ion += [  [ (  (i-L_molid[0] )*Deg_poly + j  ) * ions_per_mono + 1, i, iontype, ion_ux, ion_uy, ion_uz]  ]
            del sel1
        
        ##--- convert the result into pandas dataframe
        L_ion = pd.DataFrame(L_ion, columns = ['id', 'mol', 'type', 'xu', 'yu', 'zu'])
        return L_ion
    
    def read_ion(self, atoms, iontype, ionspermol=1):
        # read in lammps 'fix' file or any unwrapped COM data directly.
        L_ion = pd.DataFrame(atoms, columns = ['id','xu','yu','zu']) 
        L_ion['type'] = iontype
        # write mol id
        col_mol = []
        for imol in range(1, int(L_ion.shape[0]/ionspermol)+1):
            for j in range(0, ionspermol):
                col_mol += [imol]
        L_ion['mol'] = col_mol
        # L_ion.reset_index(drop=True)
        # L_ion.loc['id'] = L_ion.index+1
        return L_ion

    def delete_full(self):  # release mem
        del self.L_atom

    def export_pdb(self, f, L_sel, Nth_frame):
        f.write('REMARK\n'*3)
        f.write('CRYST1'+ '%9.3f'%self.deltaX +'%9.3f'%self.deltaY +'%9.3f'%self.deltaZ + '%7.2f'%self.alpha +'%7.2f'%self.beta+'%7.2f'%self.gama+' P 1           1\n')
        f.write('MODEL ' + ' '*4 + '%4d' %Nth_frame +'\n')
        for (idx, row) in L_sel.iterrows():
            Tatom = str(row['type'])
            if len( Tatom ) == 1:
                Tatom = Tatom +' '*3
            elif len( Tatom ) == 2:
                Tatom = Tatom +' '*2
            elif len( Tatom ) == 3:
                Tatom = Tatom +' '

            f.write('ATOM  '+ '%5d' %row['id'] + ' ' \
                + '%4s' %Tatom + ' ' + 'ION' + '  ' \
                + '%4d' %row['mol']+ '    ' \
                + '%8.3f' %row['xu']  + '%8.3f' %row['yu'] + '%8.3f' %row['zu'] \
                + ' '*24 +'\n')
                # pdb records unwrapped data
        f.write('TER\n')
        f.write('ENDMDL\n')

    def export_lmptrj(self, f, L_sel, col = \
                        ['id', 'mol', 'type', 'x', 'y', 'z', 'ix', 'iy', 'iz'] ):
        # f is the file object, i.e., pointer
        # the default columns are the wrapped data, vmd cannot rerad unwrapped-only lmp
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
        f.write('ITEM: ATOMS '+ ' '.join(col) +'\n'  )
        #f.write('ITEM: ATOMS '+ ' '.join(L_sel.columns.values.astype(str)) +'\n'  )
        # write atom body
        L_sel.loc[:,col].to_csv(f, sep=' ',  header=False, index=False)
        #float_format='%.6f',
    
    def export_gro(self, f,  L_sel):
        f.write(' t=' + str(self.time) +'\n' )
        f.write( '%d' %self.Natom +'\n')
        # write atom body
        for (idx, row) in L_sel.iterrows():
            f.write( '{:>5.5s}{:>5.5s}{:>5.5s}{:>5.5s}{:>8.3f}{:>8.3f}{:>8.3f}\n'.format(\
                str(row['mol']), row['res'], row['type'], str(row['id']), \
                row['xu'], row['yu'], row['zu'] )\
            )
            # use str() to convert the mol and id into str
            # and truncate them if too long.
            # lammps does not have a column-column trj.
            # if load lammps trj and save into gromacs,
            # it may destory the gro fromat if we have > 99999 atoms.
        f.write('  {:>9.5f}  {:>9.5f}  {:>9.5f} \n'.format(
            self.deltaX, self.deltaY, self.deltaZ      )  \
        )
    def find_asso(self, L_im, L_mo, r_cut, clean = 0):
        # input must have wrapped data x y z
        # prepare two lists for the mobile ions: associate atoms and associated molecules.
        L_mo_asso_atom = []
        L_mo_asso_mol  = []
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
        # create two arrays for number of asso. atoms and number of asso. chains for all the mobile ions
        N_asso_atom = np.array([]).astype(int)
        N_asso_mol  = np.array([]).astype(int)
        for (index,row) in L_mo.iterrows():
            N_asso_atom = np.append(N_asso_atom,  len(row['asso_atom']) ) # append the number of asso. atom for each mobile ion
            N_asso_mol  = np.append(N_asso_mol ,  len(row['asso_mol' ]) ) # append the number of asso. chain for each mobile ion
        hist_asso_im_atom = np.histogram( N_asso_atom , bins=np.arange(0, max( N_asso_atom) +2 ) ) # histo
        hist_asso_im_mol  = np.histogram( N_asso_mol  , bins=np.arange(0, max( N_asso_mol ) +2 ) ) # histo
        self.hist_asso_im_atom = hist_asso_im_atom # save stats as attributes
        self.hist_asso_im_mol  = hist_asso_im_mol
        if clean:
            L_mo = L_mo.drop(['asso_atom', 'asso_mol'], axis = 1)
            L_im = L_im.drop(['if_asso'], axis = 1)
        return  N_asso_atom, N_asso_mol # return numbers   
    def hoppingtype_AN(self, prev, selkey = ''): #current and prev snap
        asso_type = []
        try:
            looptable = self.L_AN[ self.L_AN[selkey] > 0 ]
        except:
            looptable = self.L_AN

        for ind, row in looptable.iterrows():
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
    
    def msd(self, L_mobile_ions, ref):
        return MSD(L_mobile_ions['xu'], L_mobile_ions['yu'], L_mobile_ions['zu'], ref['xu'], ref['yu'], ref['zu'])

    def nongauss(self,L_mobile_ions,ref): # non gaussian parameter data point
        return NGP(L_mobile_ions['xu'], L_mobile_ions['yu'], L_mobile_ions['zu'], ref['xu'], ref['yu'], ref['zu'])
    
    def vanhove_s(self, L_mobile_ions, ref, maxdist = 25.0, accuracy = 0.1): # VH_s data point
        return  VANHOVE_S( L_mobile_ions['xu'], L_mobile_ions['yu'], L_mobile_ions['zu'], ref['xu'], ref['yu'], ref['zu'], maxdist, accuracy ) 
    
    def vanhove_d(self, L_mobile_ions, ref, maxdist = 25.0, accuracy = 0.1): #VH_d data piont
        return  VANHOVE_D( L_mobile_ions['x'], L_mobile_ions['y'], L_mobile_ions['z'], ref['x'], ref['y'], ref['z'], [self.deltaX, self.deltaY, self.deltaZ] ,maxdist, accuracy ) 
    
    def fsqt(self, L_mobile_ions, ref, q):
        distances = np.sqrt( \
                    (  L_mobile_ions['xu']-ref['xu'])**2 \
                 +  (  L_mobile_ions['yu']-ref['yu'])**2 \
                 +  (  L_mobile_ions['zu']-ref['zu'])**2 \
                           )
        kr = distances*q
        return np.mean(np.sin(kr)/kr)

    def findfast(self, L_mobile_ions, ref, r_star=6.0):
        L_mobile_ions['fast'] = 0
        r_star_squared = r_star*r_star
        L_ds_squared = (L_mobile_ions['xu']-ref['xu'])**2 +(L_mobile_ions['yu']-ref['yu'])**2 + (L_mobile_ions['zu']-ref['zu'])**2
        for (idx, val) in L_ds_squared.iteritems():
            if val > r_star_squared:
                L_mobile_ions.loc[idx, 'fast'] = 1
        return L_mobile_ions[L_mobile_ions['fast']>0].shape[0]/L_mobile_ions.shape[0] 

    def findstring(self, L_mobile_ions, ref, cutoff=2.5, maxlength=20):
        if 'fast' in L_mobile_ions.columns:
            L_mobile_ions['string'] = -1
            L_mobile_ions.loc[L_mobile_ions['fast']>0, 'string'] = 0
            #print(L_mobile_ions)
            #print('ready to loop')
            current_string = 1
            for (idx, row) in L_mobile_ions[L_mobile_ions['fast']>0].iterrows():
                # get the coor of mobile ion:
                #print('current mob ion: \n', int(row['id'])   )
                coor_mob = np.array([row['x'], row['y'], row['z']])
                for (idx2, row2) in ref.iterrows():
                    if idx2 == idx : # skip same ion
                        continue
                    if L_mobile_ions.loc[idx2, 'fast'] <=0:
                        continue # skip the slow ion
                    if L_mobile_ions.loc[idx2,'string'] != 0 and ( L_mobile_ions.loc[idx2,'string'] == L_mobile_ions.loc[idx,'string'] ):
                        continue # skip previously combined
                    else:
                        # get the coor of ref ion:
                        coor_ref = np.array([ row2['x'], row2['y'], row2['z']] )
                        #print(coor_mob, coor_ref)
                        if ifconn(coor_mob, coor_ref, np.array([self.deltaX, self.deltaY,self.deltaZ]), cutoff):
                            ref_string_id = L_mobile_ions.loc[idx2,'string'] # read in string id of ref ion
                            #print('ions ids: ', idx+1, idx2+1, 'are in the same string!'  )
                            if L_mobile_ions.loc[idx,'string'] <=0: # mob does not in string
                                if ref_string_id >0: # ref in string
                                    L_mobile_ions.loc[idx,'string'] = ref_string_id
                                    #print('ref ion ', 'brings the string id:', ref_string_id)
                                else: # it's a new string
                                    L_mobile_ions.loc[idx, 'string'] = current_string
                                    L_mobile_ions.loc[idx2,'string'] = current_string
                                    #print('new string created: ', current_string)
                                    current_string +=1
                            elif ref_string_id <=0: # mob is in string but ref not in string
                                L_mobile_ions.loc[idx2,'string'] = L_mobile_ions.loc[idx,'string']
                                #print('mobile ion brings the string id: ', int(L_mobile_ions.loc[idx,'string'])    )
                            else: # both two are in strings, need to combine two strings
                                larger_string_id = max(L_mobile_ions.loc[idx2, 'string'], L_mobile_ions.loc[idx,'string'])
                                smaller_string_id = min(L_mobile_ions.loc[idx2, 'string'], L_mobile_ions.loc[idx,'string'])
                                L_mobile_ions.loc[ L_mobile_ions['string']==larger_string_id, 'string'] = smaller_string_id
                                #print('combine the two strings: ', larger_string_id, smaller_string_id)
                            #print( L_mobile_ions.loc[L_mobile_ions['fast']>0, ['id','x', 'y', 'z', 'fast','string']]   )
            
            stringid_histo = np.histogram(L_mobile_ions[L_mobile_ions['fast']>0]['string'] , bins=np.arange(0, L_mobile_ions.shape[0] +2) )
            ### command above: -1(slow) is excluded, stringid==0 is ns=1, stringid>0 is ns>1
            ### command below: exclude stringid==0 first, then add it back by hand.
            stringlength_histo = np.histogram(stringid_histo[0][1:], bins=np.arange(1,maxlength+2))
            stringlength_histo[0][0] =  stringid_histo[0][0]
            return stringlength_histo  # not normalized
        else:
            print('need to find fast atoms!')
            return -1

    def ht(self, r_cut):
        L_ht = []
        for (idx_an, anion)in self.L_AN.iterrows():
            for (idx_ct, cation) in self.L_CT.iterrows():
                L_ht.append(  int(  ifconn(   [ anion['x'], anion['y'], anion['z'] ], [cation['x'], cation['y'], cation['z'] ], np.array([self.deltaX, self.deltaY,self.deltaZ]), r_cut ))   )
        return L_ht

    def selectatom(self, sourcelist, L_molid, ilocL_Ter, selrule):
        L_select = []
        if len(L_molid):
            for i in L_molid:
                sel = sourcelist[ sourcelist['mol'] == i ].reset_index( drop=True )
                sel = sel.drop(ilocL_Ter).reset_index( drop = True ) 
                # Reset index in order to prepare for 
                # the selection rule based on the index(ices).
                sel1 = eval(selrule)
                if len(L_select):
                    L_select = L_select.append(sel1, ignore_index=True)
                else:
                    L_select = sel1
        else:
            L_select = sourcelist[eval(selrule)]
        return L_select.reset_index(drop=True)
        # reset_index before return
        # if not, there can be the bug, when:
        # only one molid. L_select would have wrong index 
        # orignially from sourcelist

    def bond_uw(self, sel1, sel2):
        L_b_2 = []
        for (idx1, row1) in sel1.iterrows():
            coor1 = row1[ ['xu', 'yu', 'zu'] ].values
            coor2 = sel2.iloc[idx1,:].loc[ ['xu', 'yu', 'zu'] ].values
            blength = np.linalg.norm(coor1-coor2)
            L_b_2.append(blength)
        return L_b_2
    def bond_w(self, sel1, sel2):
        L_b_2 = []
        for (idx1, row1) in sel1.iterrows():
            coor1 = row1[ ['x', 'y', 'z'] ].values
            coor2 = sel2.iloc[idx1,:].loc[ ['x', 'y', 'z'] ].values
            blength = pbc_dist(coor1, coor2, [self.deltaX, self.deltaY, self.deltaZ])
            L_b_2.append(blength)
        return L_b_2

    def sel(self, sourcelist, *args, **kwargs):
        if 'molid' in kwargs:
            pass
        else:
            pass
    
    def angle_uw(self, sel1, sel2, sel3):
        L_cos_a = []
        for (idx1, row1) in sel1.iterrows():
            coor1 = row1[ ['xu', 'yu', 'zu'] ].values
            coor2 = sel2.iloc[idx1,:].loc[ ['xu', 'yu', 'zu'] ].values
            coor3 = sel3.iloc[idx1,:].loc[ ['xu', 'yu', 'zu'] ].values
            cos_angle = np.dot( (coor1-coor2 ), ( coor3-coor2 ) ) \
                      / np.linalg.norm(coor1-coor2) /  np.linalg.norm( coor3-coor2)
            L_cos_a.append(cos_angle)
        return L_cos_a

    def vector_angle(self, sel1, sel2, sel3, sel4):
        # vector1: sel2-sel1
        # vector2: sel4-sel3
        vec1_col =   sel2[['xu', 'yu', 'zu']] \
                   - sel1[['xu', 'yu', 'zu']]
        vec2_col =   sel4[['xu', 'yu', 'zu']] \
                   - sel3[['xu', 'yu', 'zu']]
        cos_col  =   np.sum( vec1_col * vec2_col, axis = 1 ) \
                   / np.linalg.norm( vec1_col, axis = 1 )    \
                   / np.linalg.norm( vec2_col, axis = 1 )
        return cos_col

    def dihed_uw_old(self, sel1, sel2, sel3, sel4):
        L_cos_d = []
        for (idx1, row1) in sel1.iterrows():
            coor1 = row1[ ['xu', 'yu', 'zu'] ].values
            coor2 = sel2.iloc[idx1,:].loc[ ['xu', 'yu', 'zu'] ].values
            coor3 = sel3.iloc[idx1,:].loc[ ['xu', 'yu', 'zu'] ].values
            coor4 = sel4.iloc[idx1,:].loc[ ['xu', 'yu', 'zu'] ].values
            b0 = coor1 - coor2
            baxis = coor3 -coor2
            b1 = coor4 - coor3
            vaxis = baxis / np.linalg.norm( baxis )  # unit vector 
            #
            # For b0 and b1, get the component parallel to the axis.
            v0 = b0 - np.dot(b0, vaxis) * vaxis
            #  = b0 - component align with the axis
            #  = component perpendicular to the axis
            #
            # Do the same thing on b1
            v1 = b1 - np.dot(b1, vaxis) * vaxis
            #
            # v0 and v1 both perpendicular to the axis.
            # Thus, they are in the plane which is perpendicular to the axis.
            # Use dot product and arccos to determine the angle.
            #
            cos_dihedral = np.dot( v0, v1 )/ np.linalg.norm( v0 ) / np.linalg.norm( v1 )
            # 
            L_cos_d.append( cos_dihedral)
            #
        return L_cos_d
    def dihed_uw(self, sel1, sel2, sel3, sel4): # Vectorized
        # convert to np.array:
        coor1 = sel1[ ['xu', 'yu', 'zu'] ].values
        coor2 = sel2[ ['xu', 'yu', 'zu'] ].values
        coor3 = sel3[ ['xu', 'yu', 'zu'] ].values
        coor4 = sel4[ ['xu', 'yu', 'zu'] ].values
        # raw vector 1 and vector 2, axis
        b0 = coor1 - coor2
        baxis = coor3 -coor2
        b1 = coor4 - coor3
        # determine axis index: NO need to do, so commented out.
        # if only 1 dim, it will be [[xx,yy,zz]],
        # shape is (3,1)
        ind_ax = 1
        #if len(baxis.shape) == 1:
        #    ind_ax = 0
        #
        # Axis as unit vector:
        vaxis = (baxis.T /  np.linalg.norm( baxis, axis = ind_ax ) ).T
        # component of b0 and b1 perpendicular to the axis
        v0 = b0 - (vaxis.T * np.sum(b0 * vaxis, axis = ind_ax ) ).T
        v1 = b1 - (vaxis.T * np.sum(b1 * vaxis, axis = ind_ax ) ).T
        # dot product as in a row:
        cos_dih = np.sum( v0 * v1 , axis = ind_ax )
        # normalize:
        # must divide the norms of v0 and v1, respectively.
        # if v0 v1 are perpendicular, using norm of cos_dih will yield np.nan
        # np.arccos cannot get PI/2 from np.nan
        # cos_dih has n*1 dim:
        cos_dih =  (  cos_dih.T / np.linalg.norm(v0, axis = ind_ax) \
                    / np.linalg.norm(v1, axis = ind_ax) ).T
        # return cos vector as np.array
        return cos_dih

        # Test of speed for angle calc:
        # Per my test, np.dot is 40 X faster than np.cross, 
        # and 6 X faster than np.linalg.norm
        # Faster function is always preferred.
        # In numpy 1.14.2 or above, arccos is pretty good. 
        # No problem with pi/2 or pi angles. 
        # This implementation uses cos.
        # No need to do cross product or arctan.

if __name__ == '__main__':
    m = oneframe()
    a = pd.DataFrame([ [1,0,0], [2,0,0]   ], columns =['xu', 'yu', 'zu'])
    b = pd.DataFrame([ [1,1,0], [0,-5.5,0]], columns =['xu', 'yu', 'zu'])
    s1 = pd.DataFrame([[0,0,0],[0,0,0]    ], columns =['xu', 'yu', 'zu'])
    s2 = pd.DataFrame([[0,0,2],[0,0,4.2]  ], columns =['xu', 'yu', 'zu'])
    """
    a = pd.DataFrame([ [1,0,0]  ], columns =['xu', 'yu', 'zu'])
    b = pd.DataFrame([ [1,1,0]  ], columns =['xu', 'yu', 'zu'])
    s1 = pd.DataFrame([[0,0,0]  ], columns =['xu', 'yu', 'zu'])
    s2 = pd.DataFrame([[0,0,2]  ], columns =['xu', 'yu', 'zu'])
    """

    print(  m.dihed_uw(a,s1,s2, b) )
