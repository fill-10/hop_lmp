import sys
import numpy as np
import pandas as pd
import copy
from class_oneframe import *

class data(object):
    def __init__(self, savemem = 1):
        self.allframes = []
        self.save_mem = savemem
        #print('data object intialized')
    def read_all_pdb(self, filename, pbc='nojump', *args, **kwargs):
        f = open(filename, 'r')
        read_pos = 0
        
        ##--- keywords interpretation ---
        ##--- this is to save memory ---
        if 'CT_gen' in kwargs and kwargs['CT_gen']:
            self.CT_gen = kwargs['CT_gen']
        else: self.CT_gen = 0
        if 'AN_gen' in kwargs and kwargs['AN_gen']:
            self.AN_gen = kwargs['AN_gen']
        else: self.AN_gen = 0
       
        ##--- read all time frames ---
        while 1:
            try:
                ##--- read each time frame into variables
                time, Natom, box, cols, atoms, pos = read_1_pdb(f)
                onef = oneframe()
                onef.load_snap(time, box, atoms, cols)
                if self.CT_gen:
                    onef.L_CT = onef.ion_gen(*self.CT_gen)
                if self.AN_gen:
                    onef.L_AN = onef.ion_gen(*self.AN_gen)
                if self.save_mem and (len(onef.L_CT) or len(onef.L_AN) ):
                    del onef.L_atom
                    onef.L_atom = []
                if len(onef.L_AN) and len(onef.L_CT):
                    onef.L_CT['id'] += max(onef.L_AN['id'])
 
                ##--- construct data structure.
                self.allframes += [ onef]  # pass 'pandas' var into.
            except Exception as error:
                print('reading finished.', 'error:', error)
                break
        
        return 0

    def read_all_lmp(self, filename, *args, **kwargs):
        f = open(filename, 'r')
        read_pos = 0
        ##--- keywords interpretation ---
        ##--- this is to save memory ---
        if 'CT_gen' in kwargs and kwargs['CT_gen']:
            self.CT_gen = kwargs['CT_gen']
        else: self.CT_gen = 0
        if 'AN_gen' in kwargs and kwargs['AN_gen']:
            self.AN_gen = kwargs['AN_gen']
        else: self.AN_gen = 0
        ##--- loop! ---
        while 1:
            try:
                ##--- read each time frame into variables
                time, Natom, box, cols, atoms, pos = read_1_lmp(f)
                onef = oneframe()
                onef.load_snap(time, box, atoms, cols) 
                #print(onef.time)
                if self.CT_gen:
                    onef.L_CT = onef.ion_gen(*self.CT_gen)
                if self.AN_gen:
                    onef.L_AN = onef.ion_gen(*self.AN_gen)
                if self.save_mem and (len(onef.L_CT) or len(onef.L_AN) ):
                    del onef.L_atom
                    onef.L_atom = []
                if len(onef.L_AN) and len(onef.L_CT):
                    onef.L_CT['id'] +=  len(onef.L_AN)
                #print('ion generated')
                ##--- construct data structure.
                self.allframes += [ onef ]
            except Exception as error:
                print('reading finished.', 'error:', error)
                break

    def AN_CT_gen(self, L_anion_kw, L_cation_kw):
        self.CT_gen = L_cation_kw
        self.AN_gen = L_anion_kw
        for frame in self.allframes:
            if self.CT_gen:
                frame.L_CT = frame.ion_gen(*self.CT_gen)
            if self.AN_gen:
                frame.L_AN = frame.ion_gen(*self.AN_gen)
            if self.save_mem:
                del frame.L_atom
                frame.L_atom = []
            if len(frame.L_AN) and len(frame.L_CT):
                frame.L_CT['id'] +=  len(frame.L_AN)

    def AN_CT_readfix(self, AN_fixf, AN_ionspermol, CT_fixf, CT_ionspermol, box, dropuxyz=False):
        fAN = open(AN_fixf, 'r')
        fCT = open(CT_fixf, 'r')
        # clean all frames
        self.allframes = []
        # read ion
        Nframe = 0
        try:
            while 1:
                # initialize one frame 
                onef = oneframe()
                onef.update_pbc( box  )
                #print(onef.deltaZ)
                # read AN
                anions, Nanion, aniontime, anionpos = read_1_fix(fAN)
                onef.time = aniontime
                onef.L_AN = onef.read_ion(anions, 1, AN_ionspermol, dropuxyz)
                cations, Ncation, cationtime, cationpos = read_1_fix(fCT)
                onef.L_CT = onef.read_ion(cations, 2, CT_ionspermol, dropuxyz)
                onef.L_CT.loc[:, 'mol'] += max(onef.L_AN['mol'])
                Nframe +=1
                self.allframes += [onef]
                #print(Nframe)
        except Exception as err:
            print(err)
            print('ion reading completed')
            #print(Nframe)
    
    def export_ions_lmptrj(self, fn, skip=0, col = \
            ['id', 'mol', 'type', 'x', 'y', 'z', 'ix', 'iy', 'iz'] ) :
        # vmd can only read unwrapped lmptrj 
        counter = 0
        f = open(fn, 'w')
        for frame in self.allframes:
            if counter%(skip+1):
                counter +=1
                continue
            if len(frame.L_AN) and not len(frame.L_CT):
                frame.export_lmptrj( f, frame.L_AN , col )
            if not len(frame.L_AN) and len(frame.L_CT):
                frame.export_lmptrj( f, frame.L_CT , col )
            if len(frame.L_AN) and len(frame.L_CT):
                frame.export_lmptrj( f, pd.concat( \
                    [frame.L_AN.loc[:, col], frame.L_CT.loc[:, col] ], \
                    ignore_index=True ) , col )
            counter += 1
        f.close()
    def export_all_lmptrj(self, fn, skip=0, col = \
            ['id', 'mol', 'type', 'x', 'y', 'z', 'ix', 'iy', 'iz'] ) :
        # vmd can only read unwrapped lmptrj 
        counter = 0
        f = open(fn, 'w')
        for frame in self.allframes:
            if counter%(skip+1):
                counter +=1
                continue
            if len(frame.L_AN) and not len(frame.L_CT):
                frame.export_lmptrj( f, \
                    pd.concat( [frame.L_AN.loc[:, col], frame.L_atom.loc[:,col]], \
                    ignore_index=True )    )
            if not len(frame.L_AN) and len(frame.L_CT):
                frame.export_lmptrj( f, \
                    pd.concat( [frame.L_CT.loc[:, col], frame.L_atom.loc[:, col]], \
                    ignore_index=True )    )
            if len(frame.L_AN) and len(frame.L_CT):
                frame.export_lmptrj( f, \
                    pd.concat( [frame.L_AN.loc[:, col], frame.L_CT.loc[:, col], \
                    frame.L_atom.loc[:, col]], ignore_index=True  )    )
            counter += 1
        f.close()
    
    def export_ions_pdb(self, fn, skip=0):
        # pdb records unwrapped data
        col = ['id', 'mol', 'type', 'ux', 'uy', 'uz']
        counter = 0
        f = open(fn, 'w')
        for frame in self.allframes:
            if counter%(skip+1):
                counter +=1
                continue
            if len(frame.L_AN) and not len(frame.L_CT):
                frame.export_pdb(  f, frame.L_AN , counter+1 )
            if not len(frame.L_AN) and len(frame.L_CT):
                frame.export_pdb(  f, frame.L_CT , counter+1 )
            if len(frame.L_AN) and len(frame.L_CT):
                frame.export_pdb(  f, \
                    pd.concat(  [frame.L_AN.loc[:,col], \
                    frame.L_CT.loc[:,col] ], ignore_index=True  ), counter +1 )
            counter += 1
    
    def export_all_pdb(self, fn, skip=0):
        # pdb records unwrapped data
        col = ['id', 'mol', 'type', 'ux', 'uy', 'uz']
        counter = 0
        f = open(fn, 'w')
        for frame in self.allframes:
            if counter%(skip+1):
                counter +=1
                continue
            if len(frame.L_AN) and not len(frame.L_CT):
                frame.export_pdb(  f, \
                    pd.concat(  [ frame.L_AN.loc[:,col], \
                    frame.L_atom.loc[:,col] ], ignore_index=True  ) , counter+1 )
            if not len(frame.L_AN) and len(frame.L_CT):
                frame.export_pdb(  f, \
                    pd.concat(  [ frame.L_CT.loc[:,col], \
                    frame.L_atom.loc[:,vol] ], ignore_index=True  ) , counter+1 )
            if len(frame.L_AN) and len(frame.L_CT):
                frame.export_pdb(  f, \
                    pd.concat(  [ frame.L_AN.loc[:,col], \
                    frame.L_CT.loc[:,col],\
                    frame.L_atom.loc[:,col] ],\
                    ignore_index=True  ) , counter+1 )
            counter += 1

    def wrapall_L(self, skip=0 ):
        Nframe = len( self.allframes )
        for i in range(0, Nframe, skip+1):
            try:
                self.allframes[i].wrap(self.allframes[i].L_AN)
            except:
                pass
            try:
                self.allframes[i].wrap(self.allframes[i].L_CT)
            except:
                pass
            try:
                self.allframes[i].wrap(self.allframes[i].L_atom)
            except:
                pass
    def unwrapall_L(self, skip=0 ):
        Nframe = len( self.allframes )
        for i in range(0, Nframe, skip+1):
            try:
                self.allframes[i].unwrap(self.allframes[i].L_AN)
            except:
                pass
            try:
                self.allframes[i].unwrap(self.allframes[i].L_CT)
            except:
                pass
            try:
                self.allframes[i].unwrap(self.allframes[i].L_atom)
            except:
                pass

    def unwrapall_AN(self, skip=0):
        Nframe = len( self.allframes )
        for i in range(0, Nframe, skip+1):
            self.allframes[i].unwrap(self.allframes[i].L_AN)
    def unwrapall_CT(self, skip=0):
        Nframe = len( self.allframes )
        for i in range(0, Nframe, skip+1):
            self.allframes[i].unwrap(self.allframes[i].L_CT)


    def find_asso_AN_CT(self, r_cut, skip=0, clean = 0 ):
        counter = 0
        # save all frames stats together
        Total_N_asso_atom = np.array([]).astype(int)
        Total_N_asso_mol  = np.array([]).astype(int)
        
        ###--- collect number of associated atoms/mols from each frame ---
        for frame in self.allframes:
            if counter%(skip+1):
                counter +=1
                continue
            if len(frame.L_AN) and len(frame.L_CT):
                N_asso_atom_1f, N_asso_mol_1f = frame.find_asso(frame.L_CT, frame.L_AN, r_cut, clean )
                ##--- mount to the total list ---
                Total_N_asso_atom = np.append(Total_N_asso_atom, N_asso_atom_1f)
                Total_N_asso_mol  = np.append(Total_N_asso_mol , N_asso_mol_1f )
            counter += 1
        
        # stat    
        self.hist_asso_atom = np.histogram(Total_N_asso_atom, bins = np.arange(0, np.amax(Total_N_asso_atom) + 2 ) )
        
        self.hist_asso_mol  = np.histogram(Total_N_asso_mol , bins = np.arange(0, np.amax(Total_N_asso_mol ) + 2 ) )

    
        norm_hist_ttl_atom = (  self.hist_asso_atom[0] / self.hist_asso_atom[0].sum() , self.hist_asso_atom[1]  )
        norm_hist_ttl_mol  = (  self.hist_asso_mol[0]  / self.hist_asso_mol[0].sum()  , self.hist_asso_mol[1]   )
        
        return norm_hist_ttl_atom, norm_hist_ttl_mol

    def hoppingtype_AN(self, skip=0):
        Nframe = len( self.allframes )
        counter = 0
        histosum = np.array([0,0,0,0])
        # loop over frames
        for i in range(1+skip, Nframe):
            # i from second to last
            #
            if i%(skip+1):
                continue
            # p from first to last but 1
            p = i - skip -1
            histo_hop_type = self.allframes[i].hoppingtype_AN( self.allframes[p])
            histosum += histo_hop_type[0]

        norm_hist_hop_type = ( histosum/histosum.sum() , np.array([1,2,3,4,5]) )
        return norm_hist_hop_type
   
    ##--- real time non gaussian (may have fluctuations) ---
    def nongauss_AN(self, skip=0):
        Nframe = len( self.allframes )
        # prepare output columns:
        time_column = []
        nongauss_data = []
        for i in range(1, Nframe, skip+1):
            time_column += [self.allframes[i].time- self.allframes[0].time]
            nongauss_data += [ self.allframes[i].nongauss(  self.allframes[i].L_AN , self.allframes[0].L_AN  ) ]
        return time_column, nongauss_data
    ##--- averaged non gaussian ---
    def nongauss_AN_avg(self, start_interval=1, fixsave_every=10, maxattemp = 500):
        Nframe = len( self.allframes )
        time_column = []
        nongauss_data = []
        for i in range(start_interval, Nframe):
            #print(' ',i)
            time_column.append( i*fixsave_every )
            nongauss_point = [] 
            for j in range(i, min(i+maxattemp,Nframe), 1):
                nongauss_point.append( self.allframes[j].nongauss(  self.allframes[j].L_AN , self.allframes[j-i].L_AN  ) )
            #print(len(nongauss_point))
            nongauss_data.append( np.average(nongauss_point))
        return time_column, nongauss_data
    
    ##--- averaged msd ---
    def msd_AN_avg(self, start_interval=1, fixsave_every=10, maxattemp=500):
        Nframe = len( self.allframes )
        time_column = []
        msd_data = []
        for i in range(start_interval, Nframe):
            time_column.append( i*fixsave_every )
            msd_point = [] 
            for j in range(i, min(i+maxattemp,Nframe), 1):
                msd_point.append( self.allframes[j].msd(  self.allframes[j].L_AN , self.allframes[j-i].L_AN  ) )
            msd_data.append( np.average(msd_point))
        return time_column, msd_data

    ##--- averaged van hove self ---
    def vanhove_s_AN_avg(self,interval_star=100, maxdist=25.0, accuracy =0.1, skip = 0):
        Nframe = len( self.allframes )
        vanhove_s_raw =[]
        for i in range(0, Nframe - interval_star, skip+1):
            vanhove_s_point = self.allframes[i+interval_star].vanhove_s(self.allframes[i+interval_star].L_AN, self.allframes[i].L_AN, maxdist, accuracy)
            vanhove_s_raw.append(vanhove_s_point[0])
        #print(vanhove_s_raw)
        #print(np.mean(vanhove_s_raw, axis =0) )
        #print(sum(np.mean(vanhove_s_raw, axis =0)) )
        return np.arange(0, maxdist, accuracy)[:-1]+accuracy/2 , np.mean(vanhove_s_raw, axis =0)
    def fpi_r2_vanhove_s_AN_avg(self,interval_star=100, maxdist=25.0, accuracy =0.1, skip = 0):
        dist_col, vanhove_s = self.vanhove_s_AN_avg(interval_star, maxdist, accuracy, skip)
        return dist_col, vanhove_s, vanhove_s*dist_col*dist_col*4*3.14159

    ##--- averaged van hove distinct ---
    def vanhove_d_AN_avg(self,interval_star=100, maxdist=25.0, accuracy =0.1):
        Nframe = len( self.allframes )
        vanhove_d_raw =[]
        for i in range(0, Nframe - interval_star):
            vanhove_d_point = self.allframes[i+interval_star].vanhove_d(self.allframes[i+interval_star].L_AN, self.allframes[i].L_AN, maxdist, accuracy)
            vanhove_d_raw.append([vanhove_d_point[0]])
        return np.arange(0, maxdist, accuracy)[:-1]+accuracy/2 , np.mean(vanhove_d_raw, axis =0)

    ##--- find fast ---
    def find_AN_fast(self, interval_star, rstar, skip=0):
        Nframe = len( self.allframes)
        mobile_percent = []
        for i in range(0, Nframe-interval_star, skip+1):
            mobile_percent_single_p = self.allframes[i+interval_star].findfast(self.allframes[i+interval_star].L_AN, self.allframes[i].L_AN, rstar)
            mobile_percent +=[mobile_percent_single_p]
        return np.mean(mobile_percent)

    ##--- find string ---
    def find_AN_string(self, interval_star, cutoff, maxlength=20, skip=0, include_rattle_ions = False):
        Nframe = len(self.allframes)
        pns_list = []
        for i in range(0, Nframe-interval_star, skip+1):
            pns_single = self.allframes[i+interval_star].findstring(self.allframes[i+interval_star].L_AN, self.allframes[i].L_AN, cutoff)
            
            weighted_pns_single =  pns_single[0] * pns_single[1][:-1]   #weigthed histo
            #
            if include_rattle_ions: # ns =1 includes rattling and fast, so correct weighted_pns_single[0]
                weighted_pns_single[0] =  self.allframes[i+interval_star].L_AN.shape[0] - np.sum(weighted_pns_single[1:])
                total_counted_ion = self.allframes[i+interval_star].L_AN.shape[0] 
            else:
                total_counted_ion = np.sum(weighted_pns_single) 
            if total_counted_ion>0: # to avoid 0 string when not include rattleing ions.
                pns_list += [ weighted_pns_single/total_counted_ion]
                print(pns_list)
        return np.mean(pns_list, axis=0), np.arange(1, maxlength+1)

    def CtSt(self,  r_cut, skip = 0, resol = 10, maxattemp= 500): # resolution = realtime/frame
        Nframe = len( self.allframes )
        L_all_ht = []
        for i in range(0, Nframe, skip+1):
            L_all_ht += [ self.allframes[i].ht(r_cut) ]

        L_all_ht = np.array(L_all_ht) # save in numpy to speed up

        # double loop average
        time_column = [0]
        Ct_column = [1]
        St_column = [1]
        for j in range(1, L_all_ht.shape[0]): # loop over time intervals
            time_column.append( j*resol* (skip+1)  )
            Ct_raw = []
            St_raw = []

            for k in range(j, min( L_all_ht.shape[0], j+maxattemp) ):  # loop over different start frames
                Ct_raw += [ np.sum(  np.all(L_all_ht[[k-j,k],:] > 0 , axis = 0).astype(int) )/np.sum( L_all_ht[k-j,:] ) ]
                St_raw += [ np.sum(  np.all(L_all_ht[k-j:k+1,:] > 0 , axis = 0).astype(int) )/np.sum( L_all_ht[k-j,:] ) ]

            Ct_column.append( np.mean( Ct_raw)  )
            St_column.append( np.mean( St_raw)  )
        # TODO caluclate tau_c and tau_s
        from scipy.optimize import curve_fit
        

        return time_column, Ct_column, St_column

