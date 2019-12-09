
if __name__ == '__main__':
    import time as timer
    import numpy as np
    a = [[1,0,0,1,0,1]*1000, [1,0,1,1,1,1]*1000, [1,0,0,1,0,0]*1000]
    b = np.array(a)
    ##--- timer starts ---
    start = timer.perf_counter()
    for t in range(0,1000):
        n = []
        for i in range(0,6*1000):
            n.append( a[1][i]*a[2][i]*a[0][i]  )

    ##--- timer ends ---
    end = timer.perf_counter()
    print('time used: ', end-start)
    #print(n)

    ##--- timer starts ---
    start = timer.perf_counter()
    for t in range(0,1000):
        n1 = np.multiply(a[0],a[1])
        n1 = np.multiply(n1,a[2])
    ##--- timer ends ---
    end = timer.perf_counter()
    print('time used: ', end-start)



    ##--- timer starts ---
    start = timer.perf_counter()
    b = np.array(a)
    for t in range(0, 1000):
        n2 = (np.all(b!=0, axis =0)).astype(int) # axis = 0 : per column

    ##--- timer ends ---
    end = timer.perf_counter()
    print('time used: ', end-start)
    #print(n2)
    print(np.all(n1==n) )
    print(np.all(n2==n) )

