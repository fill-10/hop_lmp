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
        f.write('REMARK\n')
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

            Resname = str(row['res'])
            if len( Resname ) == 1:
                Resname = Resname + ' '*3
            elif len( Resname ) == 2:
                Resname = Resname + ' '*2
            elif len( Resname ) == 3:
                Resname = Resname + ' '*1
            elif len( Resname ) == 4:
                pass
            else:
                Resname = 'ION '

            f.write('ATOM  '+ '%5d' %row['id'] + ' ' \
                + '%4s' %Tatom + ' ' + Resname  + ' ' \
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
        f.write(' t= ' + str(self.time) +'\n' )
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

    def selectatom(self, sourcelist, L_molid=[], ilocL_Ter=[], selrule=None):
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
        # vectorize
        vec = sel1[ [ 'xu', 'yu', 'zu' ] ].values \
            - sel2[ [ 'xu', 'yu', 'zu' ] ].values
        # convert dtype from object to float.
        # the orignal dataframe has the dtypes of object,
        # because different data types in the columns.
        # next step, object would not have sqrt attr. 
        vec = vec.astype(np.float32)
        return np.sqrt( np.sum( vec **2,  axis = 1 ) ) # ufunc

    def bond_w(self, sel1, sel2):
        # vectorize
        vec = sel1[ [ 'x', 'y', 'z' ] ].values \
            - sel2[ [ 'x', 'y', 'z' ] ].values
        # convert dtype from object to float.
        # the orignal dataframe has the dtypes of object,
        # because different data types in the columns.
        # next step, object would not have sqrt attr. 
        vec = vec.astype(np.float64)
        pbc_vec( vec, [self.deltaX, self.deltaY, self.deltaZ] )
        #return np.sqrt( np.sum( vec **2,  axis = 1 ).astype(np.float64) ) # ufunc, vec has been converted to float
        return np.sqrt( np.sum( vec **2,  axis = 1 ) ) # ufunc

    def zmat(self, sel1, is_wrapped = False, is_sqrt = True ): #Z-matrix gen
        sel1 =sel1.reset_index(drop = True)
        Nrow = sel1.shape[0]
        if is_wrapped:
            nm_col = ['x', 'y', 'z']
        else:
            nm_col = ['xu', 'yu', 'zu']
        #empty res
        res = np.array([])
        for (idx, row) in sel1.iterrows():
            if idx < min(3, Nrow - 1):
                #tile ref
                ref = sel1.iloc[ idx ].loc[nm_col].values
                # other atoms
                vec_rest = sel1.iloc[idx +1 : ].loc[:, nm_col ].values
                #print (  sel1.iloc[idx : ]  )
                d_vec =   vec_rest - ref 
                # calc dist
                if is_wrapped :
                    pbc_vec( d_vec, [self.deltaX, self.deltaY, self.deltaZ] )
                else:
                    pass
                dist_ref =  np.sum( d_vec **2,  axis = 1 ) 
                res = np.concatenate( (res, dist_ref), axis = 0 ) 
        if is_sqrt:
            res = np.sqrt( res.astype(np.float32) ) #ufunc cannot loop  np.float64 (??)

        self.L_zmat = res
        # result: 0vs1~n, 1 vs 2~n, 2 vs 3~n ...
        return res
    def zmat2xyz(self, zmat_in): # rev zmat to xyz
        # [0, N-1): dist to at0
        # [N-1, N-1+N-2): dist to at1
        # [2N-3, 3N-6): dist to at2
        Nat = int ( len(zmat_in)/3 +2 )
        def build_basis3( zmat ): # 1st 3 atoms by default
            at0 = [ 0. , 0. , 0. ]
            if len (zmat) < 2 :
                return at0
            at1 = [ zmat[ 0 ] , 0. , 0. ]
            if len (zmat) < 3 :
                return [ at0, at1 ]
            # (x) **2         + y **2  = zmat[1] **2
            # (x-zmat[0]) **2 + y **2  = zmat[Nat - 1] **2
            x_at2 = 0.5 * ( ( zmat[1]**2 - zmat[Nat-1]**2 ) / zmat[0] + zmat[0] )
            y_at2 = np.sqrt( zmat[1]**2 - x_at2**2)
            at2 = [ x_at2, y_at2, 0. ]
            return [ at0, at1, at2 ] # list
        
        basis = build_basis3( zmat_in )
        res = basis
        if Nat < 4:
            return res
        #else
        for  i in range ( 3 , Nat):
            d0 = zmat_in[i-1]
            d1 = zmat_in[ Nat-1 + i-2]
            d2 = zmat_in[ 2*Nat-3 + i-3 ]
            # ( x     )**2 + (y     )**2 + z**2 = d0 **2
            # ( x - x1)**2 + (y     )**2 + z**2 = d1 **2
            # ( x - x2)**2 + (y - y2)**2 + z**2 = d2 **2
            #
            #  x1*(2x-x1) = d0**2 - d1**2 
            #  where x1 = basis[1][0]
            x_i = 0.5 * ( (d0**2-d1**2) / basis[1][0] + basis[1][0] )
            # ( x2-x1)*(2*x_i -x1-x2) + y2 *(2*y-y2)  = d1**2 -d2**2
            # where x2 = basis[2][0], x1=basis[1][0], y2=basis[2][1]
            y_i = 0.5 * ( (d1**2-d2**2 -  (basis[2][0] - basis[1][0]) * (2*x_i - basis[2][0] - basis[1][0] ) ) \
                        / basis[2][1] + basis[2][1]  )
            z_i = np.sqrt( max (0, d0 **2 - x_i**2 - y_i **2 ) )
            # avoid truncation error when z_i ~= 0 that z_i<0
            res.append( [ x_i, y_i, z_i ] )
        return res  # list(float) in (Natom, 3)

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

    def pbcsolv(self, prot, sol, rcut, di, box_lo, box_hi, L_atom ):
    # This function is to prepare the traj for sasa analysis
    # vmd sasa does not count in pbc.
    # Therefore, the traj must be modified that protein
    # is fully surrounded by solvent molecules
    # when the protein is across the pbc.
    # This function is only for the solvent box of gmx trjconv -pbc whole,
    # i.e., the solvent is mostly in the box
    # prot: protein molecules, sol: solvent molecules, 
    # rcut: cut off of the 'shield', di: 'xu'/'yu'/'zu', etc. 
    # box_lo, box_hi: box boundary in 'di' direction,
    # L_atom: whole atom list for modification. 
    # Here must have L_atom, becasue a ref variable is needed so that
    # the data can be directly modified to move the solvent molecules.
    # This function returns the copied sol molecules, because 
    # pandas.DataFrame.append cannot change the referred DataFrame L_atom.

        # upper and lower bounds needed for solvation
        bup = prot[di].max() + rcut
        blo = prot[di].min() - rcut
        # empty DataFrame for copy:
        atom_copy = pd.DataFrame( columns=sol.columns )
        # extend above upper bound:
        if  bup > box_hi: # move or add sol mol to upper bound
            sel_mol = sol.loc[  sol[ di ] <=  bup - box_hi + box_lo ,  :  ]
            sel_molid = sel_mol['mol'].drop_duplicates( ).reset_index(drop=True)
            sel_row = L_atom['mol'].isin( sel_molid ) # select from L_atom!
            # copy or move:
            if blo - box_lo > bup - box_hi : # True: move
                L_atom.loc[ sel_row, di ] += box_hi - box_lo
                pass
            else : # False: copy
                more_atom = sol.loc[ sel_row, : ].copy( deep = True )
                more_atom.loc[:, di ]  +=  box_hi - box_lo
                atom_copy = atom_copy.append( more_atom, ignore_index = True )
        # extend below lower bound:
        if blo < box_lo: #  to lower
            sel_mol = sol.loc[  sol[ di ] >=  blo + box_hi - box_lo ,  :  ]
            sel_molid = sel_mol['mol'].drop_duplicates( ).reset_index(drop=True)
            sel_row = L_atom['mol'].isin( sel_molid ) # select from L_atom!

            # copy or move:
            if box_hi - bup > box_lo - blo  : # True: move
                L_atom.loc[ sel_row, di ] -= box_hi - box_lo
                pass
            else : # False: copy
                more_atom = sol.loc[ sel_row, : ].copy( deep = True )
                more_atom.loc[:, di ]  -=  box_hi - box_lo
                atom_copy = atom_copy.append( more_atom, ignore_index = True )
        return  atom_copy
    
    def prep_sasa(self, kw_prot, kw_sol, rcut):
        ## work on self.L_atom using self.pbcsolv().
        dim = 0
        rcut_shield = rcut
        for direction in ['xu', 'yu', 'zu']:
            if dim == 0:
                box_Lo = self.xlo
                box_Hi = self.xhi
            elif dim == 1:
                box_Lo = self.ylo
                box_Hi = self.yhi
            elif dim == 2:
                box_Lo = self.zlo
                box_Hi = self.zhi
            # 
            more_copy = self.pbcsolv(          \
                self.L_atom.loc[ eval( kw_prot ), : ], \
                self.L_atom.loc[ eval( kw_sol  ), : ], \
                rcut_shield,  direction,       \
                box_Lo, box_Hi, self.L_atom    )
            self.L_atom = self.L_atom.append( more_copy, ignore_index = True )
            # counter +1 
            dim += 1
        # update id and Natom
        self.L_atom.reset_index( drop = True )
        self.L_atom['id'] = self.L_atom.index + 1
        self.Natom = self.L_atom.shape[0]
        return 0
    
    def pairwise_d_2(self, sel1, sel2 ): # pairwise distance square
        xyz1 = sel1[['x', 'y', 'z' ]].values
        xyz2 = sel2[['x', 'y', 'z' ]].values
        # separately x y z to avoid memory issues on bridges2 gpu node
        dx = np.zeros( int (xyz1.shape[0] *  xyz2.shape[0] ) , dtype=np.float32 )
        dy = np.zeros( int (xyz1.shape[0] *  xyz2.shape[0] ) , dtype=np.float32 )
        dz = np.zeros( int (xyz1.shape[0] *  xyz2.shape[0] ) , dtype=np.float32 )
        pre_rows = 0
        for i in range(0, xyz1.shape[0] ):
            x1 = xyz1[i, 0]
            y1 = xyz1[i, 1]
            z1 = xyz1[i, 2]
            raw_dx = xyz2[:, 0] - x1 # broadcast
            raw_dy = xyz2[:, 1] - y1
            raw_dz = xyz2[:, 2] - z1
            #print(raw_vec.shape[0] )
            dx[ pre_rows :  pre_rows +  raw_dx.shape[0]] = raw_dx
            dy[ pre_rows :  pre_rows +  raw_dy.shape[0]] = raw_dy
            dz[ pre_rows :  pre_rows +  raw_dz.shape[0]] = raw_dz
            pre_rows += xyz2.shape[0]
        pbc_dx = pbc_x(dx, self.deltaX)
        pbc_dy = pbc_x(dy, self.deltaY)
        pbc_dz = pbc_x(dz, self.deltaZ)
        dist_2 = pbc_dx**2 + pbc_dy**2 + pbc_dz**2
        return dist_2

    def pairwise_d_2_s(self, sel1 ): # pairwise distance square self set
        dx = np.zeros( int (xyz1.shape[0] * ( xyz1.shape[0] - 1 ) / 2 ),  dtype=np.float32 )
        dy = np.zeros( int (xyz1.shape[0] * ( xyz1.shape[0] - 1 ) / 2 ),  dtype=np.float32 )
        dz = np.zeros( int (xyz1.shape[0] * ( xyz1.shape[0] - 1 ) / 2 ),  dtype=np.float32 )
        pre_rows = 0
        for i in range(0, xyz1.shape[0] - 1 ):
            x1 = xyz1[i, 0]
            y1 = xyz1[i, 1]
            z1 = xyz1[i, 2]
            raw_dx = xyz1[i+1 :, 0] - x1 # broadcast
            raw_dy = xyz1[i+1 :, 1] - y1
            raw_dz = xyz1[i+1 :, 2] - z1
            dx[ pre_rows :  pre_rows +  raw_dx.shape[0]] = raw_dx
            dy[ pre_rows :  pre_rows +  raw_dy.shape[0]] = raw_dy
            dz[ pre_rows :  pre_rows +  raw_dz.shape[0]] = raw_dz
            pre_rows += xyz1[i+1 :].shape[0]
        # convert to pbc vector
        pbc_dx = pbc_x(dx, self.deltaX)
        pbc_dy = pbc_x(dy, self.deltaY)
        pbc_dz = pbc_x(dz, self.deltaZ)
        dist_2 = pbc_dx**2 + pbc_dy**2 + pbc_dz**2
        return dist_2

    def inter_cont_list(self, factor = 1.12246204831, *args, **kwargs):
        is_diff = False # default: same mol
        if 'L_mol0' in kwargs :
            L_mol0 = np.array( kwargs['L_mol0'] )
        else :
            L_mol0 = self.L_atom['mol']
        L_mol0 = np.unique( L_mol0 )
        L_mol0 = np.sort(L_mol0).astype(np.int32)
        if 'L_mol1' in kwargs :
            L_mol1 = np.array( kwargs['L_mol1'] )
            L_mol1 = np.unique( L_mol1 )
            L_mol1 = np.sort(L_mol1).astype(np.int32)
            is_diff = True
        sel=self.L_atom
        N_chain0 = L_mol0.shape[0]
        mask0 = sel['mol'].isin( L_mol0 )
        sel0 = sel[mask0]
        N_AperCh0 = int( sel0.shape[0] / N_chain0 )
        if is_diff:
            N_chain1 = L_mol1.shape[0]
            mask1 = sel['mol'].isin( L_mol1 )
            sel1 = sel[mask1] 
            N_AperCh1 = int( sel1.shape[0] / N_chain1 )
            N_pair = int(  N_chain0 * N_AperCh0 *  N_chain1 * N_AperCh1  )
        else :
            N_pair = int(  N_chain0 * ( N_chain0 - 1 ) / 2  * (N_AperCh0 ** 2) )
        # organize result: 
        # astype np.int8 to avoid mem problems
        L_pair_idx0 = np.zeros( N_pair , dtype= np.int16 )
        L_pair_idx1 = np.zeros( N_pair , dtype= np.int16 )
        L_d_2       = np.zeros( N_pair , dtype= np.float32 )
        L_dcut_2    = np.zeros( N_pair , dtype= np.float32 )

        #
        if is_diff:
            mask0 = sel0['mol'] == L_mol0[0]
            mask1 = sel1['mol'] == L_mol1[0]
            in_sel0 = sel0[ mask0 ]
            in_sel1 = sel1[ mask1 ]
            in_idx0 = ( in_sel0.index - np.min( in_sel0.index) ).values\
                    . astype(np.int16)
            in_idx1 = ( in_sel1.index - np.min( in_sel1.index) ).values\
                    . astype(np.int16)
            til_idx0 = np.tile( in_idx0 , N_chain0 )
            til_idx1 = np.tile( in_idx1 , N_chain1 )
            L_idx0 = np.repeat( til_idx0 , N_chain1 * N_AperCh1 )
            L_idx1 = np.tile(   til_idx1 , N_chain0 * N_AperCh0 )
            # update this into res
            L_pair_idx0[:] = L_idx0
            L_pair_idx1[:] = L_idx1
            # calc d_cut ^ 2
            # selected atoms must have 'd_cut'
            # if not, set up 'd_cut' for each atom first
            loc_pre = 0
            for imol in range(0, N_chain0 ) :
                # calc d_2 
                mask0 = sel0['mol'] == L_mol0[imol]
                d_2 = self.pairwise_d_2( sel0[mask0] , sel1 ) 
                # N of d_2 =  N_AperCh0 * N_AperCh1 * N_Chain1
                L_d_2[ loc_pre : loc_pre + N_AperCh0 * N_AperCh1 * N_chain1] = d_2
                dcut_imol = sel0[mask0]['dcut'].values
                L_dcut_idx0 = np.repeat( dcut_imol , N_chain1 * N_AperCh1 )
                L_dcut_idx1 = np.tile( sel1['dcut'] , N_AperCh0 )
                #Test: if dim of pair is correct
                #if L_dcut_idx0.shape[0] != L_dcut_idx1.shape[0] :
                #    print( 'error dim')
                #    return -1
                #else : print ( L_dcut_idx1.shape[0]   )
                #---
                #mix rule, 2^(1/6) and square
                dcut_2 = ( 0.5 * ( L_dcut_idx0 + L_dcut_idx1 ) * factor ) ** 2
                N_cur = N_AperCh0 * N_AperCh1 * N_chain1
                # Test:
                #if N_cur != L_dcut_idx1.shape[0]:
                #    print('error dim N_cur' )
                L_dcut_2[ loc_pre : loc_pre + N_cur ]  =  dcut_2
                loc_pre += N_cur

        else : # is_diff = False
            loc_pre = 0
            for imol in range(0, N_chain0-1) :
                # select
                mask0 = sel0['mol'] == L_mol0[imol]
                mask1 = sel0['mol'] >  L_mol0[imol]
                in_sel0 = sel0[ mask0 ]
                in_sel1 = sel0[ mask1 ]
                # build internal index
                in_idx0 = ( in_sel0.index - np.min( in_sel0.index) ).values\
                        . astype(np.int16)
                in_idx1 = np.tile(  in_idx0, N_chain0-1-imol  ).astype(np.int16)
                # internal lindex list of all pairs
                L_idx0 = np.repeat( in_idx0 , in_idx1.shape[0] )
                L_idx1 = np.tile(   in_idx1 , in_idx0.shape[0] )
                # update this into res
                L_pair_idx0[ loc_pre : loc_pre + L_idx0.shape[0] ]  =  L_idx0
                L_pair_idx1[ loc_pre : loc_pre + L_idx1.shape[0] ]  =  L_idx1
                # calc d_2 
                d_2 = self.pairwise_d_2( in_sel0 , in_sel1 )
                L_d_2      [ loc_pre : loc_pre + L_idx1.shape[0] ]  =  d_2
                #

                # calc d_cut ^ 2
                # selected atoms must have 'd_cut'
                # if not, set up 'd_cut' for each atom first
                dcut_idx0 = in_sel0['dcut'].values
                dcut_idx1 = in_sel1['dcut'].values
                L_dcut_idx0 = np.repeat( dcut_idx0 , dcut_idx1.shape[0] )
                L_dcut_idx1 = np.tile(   dcut_idx1 , dcut_idx0.shape[0] )
                #Test: if dim of pair is correct
                #if L_dcut_idx0.shape[0] != L_dcut_idx1.shape[0] :
                #    print( 'error dim')
                #    return -1
                #else : print ( L_dcut_idx1.shape[0]   )
                #
                #
                #mix rule, 2^(1/6) and square
                dcut_2 = ( 0.5 * ( L_dcut_idx0 + L_dcut_idx1 ) * factor ) ** 2
                L_dcut_2   [ loc_pre : loc_pre + L_dcut_idx1.shape[0] ]  =  dcut_2
                # update previous loop idx locator in res
                loc_pre += L_idx1.shape[0]

        # contact numbers
        mask_contact = L_d_2 < L_dcut_2
        L_pair_idx0 = L_pair_idx0 [ mask_contact ]
        L_pair_idx1 = L_pair_idx1 [ mask_contact ]
        #print( len(L_d_2[L_d_2>0]))
        #print ('contact pairs: ',  L_pair_idx1. shape )
        return L_pair_idx0, L_pair_idx1, L_d_2 [ mask_contact ], L_dcut_2 [ mask_contact ]

    def sel_mol( self, L_mol):
        mask = self.L_atom['mol'].isin( L_mol)
        return self.L_atom[mask]

if __name__ == '__main__':
    m = oneframe()
    a = pd.DataFrame([ [1,0,0], [2,0,0]   ], columns =['xu', 'yu', 'zu'])
    b = pd.DataFrame([ [1,1,0], [0,-5.5,0]], columns =['xu', 'yu', 'zu'])
    s1 = pd.DataFrame([[0,0,0],[0,0,0]    ], columns =['xu', 'yu', 'zu'])
    s2 = pd.DataFrame([[0,0,2],[0,0,4.2]  ], columns =['xu', 'yu', 'zu'])
    print(  m.dihed_uw(a,s1,s2, b) )
