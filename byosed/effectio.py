import numpy as np
import pandas,glob,sncosmo,os,sys
from scipy.interpolate import interpn,interp1d
from sncosmo.constants import HC_ERG_AA, MODEL_BANDFLUX_SPACING
from sncosmo.utils import integration_grid
from astropy.io import ascii
from astropy.table import Table

__all__ = ['generate_ND_grids','read_ND_grids','kaepora_to_sed']

def kaepora_to_sed(data_folder,effect_keys,base_sed='hsiao',minWave=0,maxWave=np.inf,minPhase=-np.inf,
                   maxPhase=np.inf,waveStep=10,scale_band='bessellb'):
    band=sncosmo.get_bandpass(scale_band)

    wave, dwave = integration_grid(band.minwave(), band.maxwave(),
                               MODEL_BANDFLUX_SPACING)
    trans = band(wave)
    if isinstance(base_sed,str):
        base_sed=sncosmo.Model(base_sed)

    filelists=[]
    for k in effect_keys:
        filelists.append(glob.glob(os.path.join(data_folder,'*%s*'%k)))
    seds={}
    for i in range(len(filelists)):
        filelist=filelists[i]
        temp_key=effect_keys[i]
        temp_phase=[]
        temp_min_wave=[]
        temp_max_wave=[]
        temp_flux=[]
        for f in filelist:
            dat=ascii.read(f)
            #print(dat['Flux'][0],dat['1-Sigma Lower'][0],dat['1-Sigma Upper'][0])


            phase=f[f.find('phase=')+6:f.find('phase=')+6+(f[f.find('phase=')+6:]).find('_')]

            phase=-1*float(phase[1:]) if phase[0]=='m' else float(phase[1:])

            temp_phase.append(phase)

            dat=dat[dat['Wavelength']>=minWave]
            dat=dat[dat['Wavelength']<=maxWave]
            temp_min_wave.append(np.min(dat['Wavelength']))
            temp_max_wave.append(np.max(dat['Wavelength']))


            temp_flux.append({'wave':dat['Wavelength'],
                            'flux':interp1d(dat['Wavelength'],dat['Flux'])})

        if len(np.unique(temp_phase))!=len(temp_phase):
            print('You have more than one file for at least one phase.')
            sys.exit(1)
        low_w=np.max(temp_min_wave)
        high_w=np.min(temp_max_wave)
        final_wavelength=np.arange(low_w,high_w+waveStep/10,waveStep)
        temp_phase=np.array(temp_phase)
        bound_inds=np.where(np.logical_and(np.where(temp_phase>=minPhase)[0],np.where(temp_phase<=maxPhase)[0]))[0]
        final_phase=temp_phase[bound_inds]
        temp_flux=[temp_flux[j] for j in bound_inds]
        phase_inds=np.argsort(final_phase)


        seds[temp_key]=[np.sort(final_phase),final_wavelength,sncosmo.TimeSeriesSource(final_phase[phase_inds],final_wavelength,
                                                              np.array([temp_flux[p]['flux'](final_wavelength)/\
                                                              (np.sum(wave * trans * temp_flux[p]['flux'](wave).flatten()) * \
                                                              dwave / HC_ERG_AA)*base_sed.bandflux(scale_band,final_phase[p]) \
                                                              for p in phase_inds]))]


    return(seds)

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

def line_prepender(filename, line):
    with open(filename, 'r+') as f:
        content = f.read()

        f.seek(0, 0)
        f.write(line.rstrip('\r\n') + '\n' + content)

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
        #Table(gridded).write(filename,format='ascii',header=header)
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

    theta=np.array(gridded[gridded.columns[-1]]).reshape(dim)*scale_factor#-1.


    return([x.upper() for x in gridded.columns][:-1],lambda interp_array:interpn(arrs,theta,xi=interp_array,method='linear',bounds_error=False,fill_value=0))