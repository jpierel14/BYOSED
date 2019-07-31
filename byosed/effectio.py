import numpy as np
import pandas
from scipy.interpolate import interpn


__all__ = ['generate_ND_grids','read_ND_grids']

def _meshgrid2(*arrs):
    arrs = tuple(arrs)	#edit
    lens = list(map(len, arrs))
    dim = len(arrs)

    sz = 1
    for s in lens:
        sz*=s

    ans = []
    for i, arr in enumerate(arrs):
        slc = [1]*dim
        slc[i] = lens[i]
        arr2 = np.asarray(arr).reshape(slc)
        for j, sz in enumerate(lens):
            if j!=i:
                arr2 = arr2.repeat(sz, axis=j)
        ans.append(arr2)

    return tuple(ans)


def generate_ND_grids(func,filename=None,colnames=None,*arrs):
    g=_meshgrid2(*arrs)
    positions = np.vstack(list(map(np.ravel, g))).T
    res=func(*(positions[:,i] for i in range(positions.shape[1]))).reshape((positions.shape[0],1))
    gridded=np.hstack([positions,res])
    if filename is not None:
        if colnames is not None:
            header=' '.join(colnames)
        else:
            header=''
        np.savetxt(filename,gridded,fmt='%f',header=header)
    return(gridded)


def read_ND_grids(filename,scale_factor=1.):
    with open(filename,'r') as f:
        temp=f.readline()

        if temp[0]=='#':
            names=temp.strip('#').split()
            gridded=pandas.read_csv(filename,sep=' ',names=names,comment='#',header=None)
        else:
            gridded=pandas.read_csv(filename,sep=' ',comment='#',header=None)

    arrs=tuple(np.unique(gridded.values[:,i]) for i in range(len(gridded.columns)-1))

    dim=[len(x) for x in arrs]

    theta=np.array(gridded[gridded.columns[-1]]).reshape(dim)*scale_factor-1.


    return([x.upper() for x in gridded.columns][:-1],lambda interp_array:interpn(arrs,theta,xi=interp_array,method='linear',bounds_error=False,fill_value=0))