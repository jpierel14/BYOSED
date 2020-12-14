import numpy as np
import pandas,glob,sncosmo,os,sys,scipy
from scipy.interpolate import interpn,interp1d
from scipy.signal import savgol_filter
from sncosmo.constants import HC_ERG_AA, MODEL_BANDFLUX_SPACING
from sncosmo.utils import integration_grid
from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt

__all__ = ['generate_ND_grids','read_ND_grids','kaepora_to_sed']

def kaepora_to_sed(data_folder,effect_keys,base_sed='hsiao',minWave=0,maxWave=np.inf,minPhase=-np.inf,
                   maxPhase=np.inf,waveStep=10,scale_band='bessellb'):
    waves={}
    trans={}
    dwaves={}
    for b in ['bessellb','bessellv','bessellr']:
        band=sncosmo.get_bandpass(b)    
        wave, dwave = integration_grid(band.minwave(), band.maxwave(),
                               MODEL_BANDFLUX_SPACING)
        waves[b]=wave
        dwaves[b]=dwave
        trans[b]=band(wave)
    if isinstance(base_sed,str):
        base_sed=sncosmo.Model(base_sed)

    base_sed_wave=base_sed._source._wave

    filelists=[]
    for k in effect_keys:
        filelists.append(glob.glob(os.path.join(data_folder,'*%s*'%k)))

    seds={}
    phase_pairs={}
    phase_lists={}
    scale_factors={}
    for i in range(len(filelists)):
        filelist=filelists[i]
        temp_key=effect_keys[i]
        temp_phase=[]
        
        j=0
        for f in filelist:
            dat=ascii.read(f)

            phase=f[f.find('phase=')+6:f.find('phase=')+6+(f[f.find('phase=')+6:]).find('_')]

            phase=-1*float(phase[1:]) if phase[0]=='m' else float(phase[1:])

            temp_phase.append(phase)

            dat=dat[dat['Wavelength']>=minWave]
            dat=dat[dat['Wavelength']<=maxWave]

            if i ==0:
                key_phase=np.round(phase,2)
                phase_pairs[key_phase]=[{'wave':dat['Wavelength'],
                            'flux':dat['Flux'],'interp':interp1d(dat['Wavelength'],dat['Flux'],kind='cubic')}]
                
            else:
                pair_list=np.array(list(phase_pairs.keys()))
                match=pair_list[np.abs(pair_list-phase).argmin()]
                key_phase=np.round(match,2)
                phase_pairs[key_phase].append({'wave':dat['Wavelength'],
                            'flux':dat['Flux'],'interp':interp1d(dat['Wavelength'],dat['Flux'],kind='cubic')})
                
            if np.min(np.abs(dat['Flux']-10))<.001:
                scale_factors[key_phase]=_get_kaepora_scale(phase_pairs[key_phase][i]['interp'],
                    base_sed,waves,trans,dwaves,key_phase)

        if len(np.unique(temp_phase))!=len(temp_phase):
            print('You have more than one file for at least one phase.')
            sys.exit(1)

        phase_lists[temp_key]=temp_phase
        
    final_wavelength=np.arange(base_sed_wave[0],base_sed_wave[-1]+waveStep/10,waveStep)
    to_save={}
    for j in range(len(effect_keys)):
        
        temp_phase=np.array(phase_lists[effect_keys[j]])
        bound_inds=np.where(np.logical_and(temp_phase>=minPhase,temp_phase<=maxPhase))[0]
        final_phase=np.sort(temp_phase[bound_inds])
        pair_list=np.array(list(phase_pairs.keys()))
        match=np.round(np.array([pair_list[np.abs(pair_list-phase).argmin()] for phase in final_phase]),2)
        
        final_flux=[]
        for i in range(len(match)):
            scaled_flux=interp1d(phase_pairs[match[i]][j]['wave'],
                phase_pairs[match[i]][j]['flux']*scale_factors[match[i]],kind='cubic')
            to_save[final_phase[i]]=scale_factors[match[i]]
                                                  
            temp=np.zeros(len(final_wavelength))
            for w in range(len(final_wavelength)):
                if final_wavelength[w]<np.min(phase_pairs[match[i]][j]['wave']) or final_wavelength[w]>np.max(phase_pairs[match[i]][j]['wave']):
                    temp[w]=base_sed._source._flux(final_phase[i],final_wavelength[w])
                else:
                    temp[w]=scaled_flux(final_wavelength[w])
            final_flux.append(temp)
        seds[effect_keys[j]]=[np.sort(final_phase),final_wavelength,sncosmo.TimeSeriesSource(final_phase,final_wavelength,
                                                              np.array(final_flux))]
     
    return(seds)

def _band_min(scale,args):
    """
    Minimization function for _get_kaepora_scale
    """

    ktotal,btotal=args
    return(np.abs(scale*ktotal-btotal))

def _get_kaepora_scale(kinterp,bsed,waves,trans,dwaves,phase):
    """
    Minimizes difference between BVR filters for baseline SED and
    kaepora composite spectra.
    """

    kflux=0
    bflux=0

    bsed['amplitude']=1/np.max(bsed.flux(phase,waves['bessellb']))

    for b in ['bessellb','bessellv','bessellr']:
        kflux+=(np.sum(waves[b] * trans[b] * kinterp(waves[b]).flatten())*dwaves[b] / HC_ERG_AA)
        bflux+=bsed.bandflux(b,phase)
    x=scipy.optimize.minimize(_band_min,
        [.1],
        [kflux,bflux],bounds=[[0,10]])
    final_scale=x.x/bsed['amplitude']
    bsed['amplitude']=1
    return final_scale

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