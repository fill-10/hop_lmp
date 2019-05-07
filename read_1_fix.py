def read_1_fix(f, start=0, dim=3):
    if start>0:
        f.seek(start)
    # skip head lines
    start = f.tell()
    for nheadline in range(start,3):
        f.readline()
    cline = f.readline().split()
    time = float( cline[0] )
    Nchunk = int( cline[1] )
    atoms = []
    if dim ==3:
        for l in range(0, Nchunk):
            cline = f.readline().split()
            # skip the empty lines
            if cline[1] == '0' and cline[2] == '0' and cline[3] == '0':
                continue
            else:
                atoms.append(  [ int(cline[0]), float(cline[1]), float(cline[2]), float(cline[3]) ]  )
    """
    elif dim==2: #2d is not implemented..
        return 0
    """
    return atoms, len(atoms), time, f.tell()


if __name__ == '__main__' :
    filename = './com_CT_every10ps_0-10ns.dat'
    f = open(filename, 'r')
    atoms, Nion, time, pos = read_1_fix(f)
    #print(time)
    #print(Nion)
    #print(atoms)
    atoms, Nion, time, pos = read_1_fix(f)
    print(time)
    #print(Nion)
    #print(atoms)
    f.close()
