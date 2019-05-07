import sys
import numpy as np
import pandas as pd

from class_oneframe import *


class data(object):
    def __init__(self, pdbfilename, save_mem=1):
        self.allframes = []
        self.atom_dict = None
        self.filename = pdbfilename
        self.CT_gen = None
        self.AN_gen = None
        self.save_mem = save_mem
        #print('data object intialized')
    def readtotal(self, *args, **kwargs):
        f = open(self.filename, 'r')
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
                atomlist, Natom, box, time, cols, position  = read_snap(f)  # pass the 'translate' var into read func.
                onef = oneframe()
                onef.load_snap(time, box, atomlist, cols) 
                
                if self.CT_gen:
                    onef.L_CT = onef.ion_gen(*self.CT_gen)
                if self.AN_gen:
                    onef.L_AN = onef.ion_gen(*self.AN_gen)
                if self.save_mem and (len(onef.L_CT) or len(onef.L_AN) ):
                    del onef.atoms
                if len(onef.L_AN) and len(onef.L_CT):
                    onef.L_CT['atom'] += max(2000, len(onef.L_AN))

                ##--- construct data structure.
                self.allframes += [ onef]  # pass 'pandas' var into.

            except Exception as error:
                print('reading finished.', 'error:', error)
                break
        
        ##--- release mem ---
        return 0

    def AN_CT_gen(self, L_anion_kw, L_cation_kw):
        self.CT_gen = L_cation_kw
        self.AN_gen = L_anion_kw
        for frame in self.allframes:
            if self.CT_gen:
                frame.L_CT = frame.ion_gen(*self.CT_gen)
            if self.AN_gen:
                frame.L_AN = frame.ion_gen(*self.AN_gen)
            if self.save_mem:
                del frame.atoms
            if len(frame.L_AN) and len(frame.L_CT):
                frame.L_CT['atom'] += max(2000, len(frame.L_AN))
    def AN_CT_readfix(self, AN_fixf, AN_ionspermol, CT_fixf, CT_ionspermol, box):
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
                onef.L_AN = onef.read_ion(anions, 1, AN_ionspermol)
                cations, Ncation, cationtime, cationpos = read_1_fix(fCT)
                onef.L_CT = onef.read_ion(cations, 2, CT_ionspermol)
                onef.L_CT.loc[:, 'mol'] += max(onef.L_AN['mol'])
                Nframe +=1
                self.allframes += [onef]
                #print(Nframe)
        except Exception as err:
            print(err)
            print('ion reading completed')
            #print(Nframe)
    
    def export_ions_lmptrj(self, fn, skip=0):
        counter = 0
        f = open(fn, 'w')
        for frame in self.allframes:
            if counter%(skip+1):
                counter +=1
                continue
            if len(frame.L_AN) and not len(frame.L_CT):
                frame.export_lmptrj( f, frame.L_AN  )
            if not len(frame.L_AN) and len(frame.L_CT):
                frame.export_lmptrj( f, frame.L_CT )
            if len(frame.L_AN) and len(frame.L_CT):
                frame.export_lmptrj( f, pd.concat(  [frame.L_AN, frame.L_CT], ignore_index=True  ) )
            counter += 1
        f.close()



    def export_ions_pdb(self, fn_prefix, skip=0):
        counter = 0
        for frame in self.allframes:
            if counter%(skip+1):
                counter +=1
                continue
            if len(frame.L_AN) and not len(frame.L_CT):
                frame.export_pdb(  fn_prefix+str(frame.time)+'.pdb', frame.L_AN, counter  )
            if not len(frame.L_AN) and len(frame.L_CT):
                frame.export_pdb(  fn_prefix+str(frame.time)+'.pdb', frame.L_CT, counter  )
            if len(frame.L_AN) and len(frame.L_CT):
                frame.export_pdb(  fn_prefix+str(frame.time)+'.pdb', pd.concat(  [frame.L_AN, frame.L_CT], ignore_index=True  ), counter  )
            counter += 1

    def find_asso_AN_CT(self, r_cut, skip=0, calc_stat=True):
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
                N_asso_atom_1f, N_asso_mol_1f = frame.find_asso(frame.L_CT, frame.L_AN, r_cut, calc_stat )
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
            
        



if __name__ == '__main__':
    import time as timer
    ##--- read-in the whole file ---
    filename = ''
    anionfixfile = '../400K_corrected_lmp/com_AN_every10ps_0-10ns.dat'
    cationfixfile = '../400K_corrected_lmp/com_CT_every10ps_0-10ns.dat'
    anionfixfile = '../com_AN_every10ps_0-10ns.dat'
    cationfixfile = '../com_CT_every10ps_0-10ns.dat'
    d1 = data(filename)
    start = timer.perf_counter()
    
    d1.AN_CT_readfix(  anionfixfile, 1, cationfixfile, 40, [ [-0.02,58.5387],[-0.02,58.5387],[0.,58.5187]  ]  )
    d1.export_ions_lmptrj('ions.lammpstrj')

    ##--- generate the ions ---
    # Here are the most tricky lines as follows:
    # The terminal atoms to knock out, the %(mod) atoms to include(key words in the " "):
    # The numbers are the number of rows, i.e. the indices in pandas.DataFrame. 
    # They are equal to the corresponding atomid - 1 !!! 
    # e.g. knock out terminal [18, 781] means knock out atomid 19 and 782.
    # index%20 == 5 means atomid%20 == 6.
    #CT_gen_kw = [range(401, 401+10),'CT', 40, 1, [18,781], 8, "sel[  (  (sel.index%20 <=5) & (sel.index%20>=1) ) | ( (sel.index%20>=9 ) & (sel.index%20<=11)) ]"]
    #AN_gen_kw = [range(1, 1+400), 'AN', 1, 1, [], 15, 'sel[:]']
    
    ##--- calcuate histograms ---
    hist_atom, hist_mol = d1.find_asso_AN_CT(7.8)
    print(d1.allframes[-1].L_AN.loc[ d1.allframes[-1].L_CT['id'] ==782])
    
    print(hist_atom)
    print(hist_mol)

    norm_hist_hop_type = d1.hoppingtype_AN()
    print(norm_hist_hop_type)

    stop  = timer.perf_counter()
    
    ##--- save pdb files of ions for each frame ---
    
    ##--- output to screen ---
    #print( hist_atom, '\n'*2, hist_mol   )
    #print(norm_hist_hop_type )
    #print(d1.hist_hop_type)
    #print( d1.hist_asso_atom)
    #print( d1.hist_asso_mol )
    print(stop-start)
