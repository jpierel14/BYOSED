from astropy.table import Table
import glob
import sys
import matplotlib.pyplot as plt
import numpy as np

keys=['high','low']
for key in keys:
    files=glob.glob('*%s*.txt'%key)
    dms=[]
    for f in files:
        dat=Table.read(f,format='ascii')
        dms.append(np.median(dat['Dm15']))
    print(np.median(dms))
    plt.hist(dms)
    plt.savefig(key,format='pdf')
    plt.close()
