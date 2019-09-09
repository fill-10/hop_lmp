from class_data import data
if __name__ == '__main__':
    import time as timer
    ##--- read-in the whole file ---
    pdbfilename = '/work/04201/tg835006/stampede2/LiuSum/C5VIm/570/nvtprod_every100ps_00-50ns.pdb'
    ##--- timer start ---
    start = timer.perf_counter()
    
    ##--- initialize a data object ---
    d1 = data()
    d1.save_mem=0 # save atoms and ions

    ##--- read in pdb traj output ---
    AN_gen_kw =  [ range(1,1+400), 1, 1, 1, [], 15, "sel[:]", 'type']
    CT_gen_kw =  [ range(401,401+10), 2, 40, 1, [27, 1132], 8, "sel[ ((sel.index%29<=5 ) & (sel.index%29>=1)) | ( (sel.index%29>=12 ) & (sel.index%29<=14)) ]", 'type'   ]
    d1.read_all_pdb(pdbfilename, AN_gen = AN_gen_kw, CT_gen =CT_gen_kw )
    #print(d1.allframes[0].L_AN)
    print('last timestep read: ', d1.allframes[-1].time)
    for frame in d1.allframes: #adjust id and mol
        frame.L_AN['id'] += max(frame.L_atom['id'])
        frame.L_AN['mol'] += max(frame.L_atom['mol'])
        frame.L_CT['id'] += max(frame.L_atom['id'])
        frame.L_CT['mol'] += max(frame.L_atom['mol'])
    d1.wrapall_L()
    ##--- save lammpstrj of ions to calculate rdf---
    d1.export_all_lmptrj('/work/04201/tg835006/stampede2/LiuSum/C5VIm/570/ions_every100ps_0-50ns.lammpstrj')
    #d1.export_all_pdb('ions_every1ns_0-30ns.pdb')
    print('ions written in lammpstrj')

    ##--- timer stop ---
    stop = timer.perf_counter()
    print('time used in sec: ', stop-start)
