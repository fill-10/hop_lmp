from class_data import data
if __name__ == '__main__':
    import time as timer
    ##--- read-in the whole file ---
    filename = '../nvtprodB_every100ps_500frames_0-50ns.pdb'
    d1 = data(filename)
    
    start = timer.perf_counter()
    
    ##--- generate the ions ---
    # Here are the most tricky lines as follows:
    # The terminal atoms to knock out, the %(mod) atoms to include(key words in the " "):
    # The numbers are the number of rows, i.e. the indices in pandas.DataFrame. 
    # They are equal to the corresponding atomid - 1 !!! 
    # e.g. knock out terminal [18, 781] means knock out atomid 19 and 782.
    # index%20 == 5 means atomid%20 == 6. 
    CT_gen_kw = [range(1, 1+10),'CT', 40, 1, [30,1171], 15, "sel[  (  (sel.index%30 <=9) & (sel.index%30>=5) ) | ( (sel.index%30>=15 ) & (sel.index%30<=24)) ]" ]
    AN_gen_kw = [range(11, 11+400), 'AN', 1, 1, [], 15, 'sel[:]']
    d1.readtotal( CT_gen=CT_gen_kw, AN_gen = AN_gen_kw ) 
    
    ##--- calcuate histograms ---
    #hist_atom, hist_mol = d1.find_asso_AN_CT(7.8)
    #norm_hist_hop_type = d1.hopping_types_AN()


    stop  = timer.perf_counter()
    
    ##--- save pdb files of ions for each frame ---
    d1.export_ions_pdb('AmC2_B_every100_0-50ns_')
    
    ##--- output to screen ---
    #print( hist_atom, '\n'*2, hist_mol   )
    #print(norm_hist_hop_type )
    #print(d1.hist_hop_type)
    #print( d1.hist_asso_atom)
    #print( d1.hist_asso_mol )
    print(stop-start)
