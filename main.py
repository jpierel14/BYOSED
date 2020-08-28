
#func=lambda x,y,z:np.array(x+y+z)
#_generate_ND_grids(func,'test_host.dat',['phase','mass','velocity','theta'],np.array([0,20]),np.array([2,3]),np.array([4,5]))
#sys.exit()
#p,test=_read_ND_grids('salt2_m0.dat')
#print(test(np.array([[10,5000],[10,6000]])))

import byosed,sncosmo,sys
import numpy as np
import matplotlib.pyplot as plt
# sp,sw,sf=sncosmo.read_griddata_ascii('../salt3_testing/gridded_salt2_stretch.dat')
# eff=sncosmo.TimeSeriesSource(sp,sw,sf)
# hp,hw,hf=sncosmo.read_griddata_ascii('../salt3_testing/test_mangled.dat')
# base=sncosmo.TimeSeriesSource(hp,hw,hf)
# m1=sncosmo.Model('salt2-extended')._source._model['M1']
# #print(base._flux(0,sw).flatten()*(1+eff._flux(0,sw).flatten()/base._flux(0,sw).flatten()))
# plt.plot(sw,eff._flux(0,sw).flatten())#base._flux(0,sw).flatten()*(1+eff._flux(0,sw).flatten()/base._flux(0,sw).flatten()))
# plt.show()
# sys.exit()
def ccm_extinction(wave, ebv, r_v=3.1):
    """
    The extinction (A_lambda) for given wavelength (or vector of wavelengths)
    from the CCM 1989 (and O'Donnell 1994) parameterization. Returns an
    array of extinction values for each wavelength in 'wave'
    INPUTS:
    wave - array of wavelengths (in Angstroms)
    ebv - colour excess E(B-V) float. If a negative ebv is supplied
          fluxes will be reddened rather than dereddened
    OPTIONAL INPUT:
    r_v - float specifying the ratio of total selective
          extinction R(V) = A(V)/E(B-V). If not specified,
          then r_v = 3.1
    OUTPUTS:
    funred - unreddened calibrated flux array, same number of
             elements as wave
    NOTES:
    1. This function was converted from the IDL Astrolib procedure
       last updated in April 1998. All notes from that function
       (provided below) are relevant to this function
    2. (From IDL:) The CCM curve shows good agreement with the Savage & Mathis (1979)
       ultraviolet curve shortward of 1400 A, but is probably
       preferable between 1200 and 1400 A.
    3. (From IDL:) Many sightlines with peculiar ultraviolet interstellar extinction
       can be represented with a CCM curve, if the proper value of
       R(V) is supplied.
    4. (From IDL:) Curve is extrapolated between 912 and 1000 A as suggested by
       Longo et al. (1989, ApJ, 339,474)
    5. (From IDL:) Use the 4 parameter calling sequence if you wish to save the
       original flux vector.
    6. (From IDL:) Valencic et al. (2004, ApJ, 616, 912) revise the ultraviolet CCM
       curve (3.3 -- 8.0 um-1).    But since their revised curve does
       not connect smoothly with longer and shorter wavelengths, it is
       not included here.
    7. For the optical/NIR transformation, the coefficients from
       O'Donnell (1994) are used
    >>> ccm_unred([1000, 2000, 3000], [1, 1, 1], 2 )
    array([9.7976e+012, 1.12064e+07, 32287.1])
    """
    import numpy as np
    scalar = not np.iterable(wave)
    if scalar:
        wave = np.array([wave], float)
    else:
        wave = np.array(wave, float)

    x = 10000.0/wave
    npts = wave.size
    a = np.zeros(npts, float)
    b = np.zeros(npts, float)

    #Infrared
    good = np.where( (x > 0.3) & (x < 1.1) )
    a[good] = 0.574 * x[good]**(1.61)
    b[good] = -0.527 * x[good]**(1.61)

    # Optical & Near IR
    good = np.where( (x  >= 1.1) & (x < 3.3) )
    y = x[good] - 1.82

    c1 = np.array([ 1.0 , 0.104,   -0.609,    0.701,  1.137,
                    -1.718,   -0.827,    1.647, -0.505 ])
    c2 = np.array([ 0.0,  1.952,    2.908,   -3.989, -7.985,
                    11.102,    5.491,  -10.805,  3.347 ] )

    a[good] = np.polyval(c1[::-1], y)
    b[good] = np.polyval(c2[::-1], y)

    # Mid-UV
    good = np.where( (x >= 3.3) & (x < 8) )
    y = x[good]
    F_a = np.zeros(np.size(good),float)
    F_b = np.zeros(np.size(good),float)
    good1 = np.where( y > 5.9 )

    if np.size(good1) > 0:
        y1 = y[good1] - 5.9
        F_a[ good1] = -0.04473 * y1**2 - 0.009779 * y1**3
        F_b[ good1] =   0.2130 * y1**2  +  0.1207 * y1**3

    a[good] =  1.752 - 0.316*y - (0.104 / ( (y-4.67)**2 + 0.341 )) + F_a
    b[good] = -3.090 + 1.825*y + (1.206 / ( (y-4.62)**2 + 0.263 )) + F_b

    # Far-UV
    good = np.where( (x >= 8) & (x <= 11) )
    y = x[good] - 8.0
    c1 = [ -1.073, -0.628,  0.137, -0.070 ]
    c2 = [ 13.670,  4.257, -0.420,  0.374 ]
    a[good] = np.polyval(c1[::-1], y)
    b[good] = np.polyval(c2[::-1], y)

    # Defining the Extinction at each wavelength
    a_v = r_v * ebv
    a_lambda = a_v * (a + b/r_v)
    if scalar:
        a_lambda = a_lambda[0]
    return a_lambda
'''
cp,cw,cl=sncosmo.read_griddata_ascii('byosed/initfiles/salt2_colorlaw.dat')

hp,hw,hf=sncosmo.read_griddata_ascii('../salt3_testing/test_mangled.dat')
base=sncosmo.TimeSeriesSource(hp,hw,hf)

cl=10**(-.4*sncosmo.Model('salt2-extended')._source._colorlaw(cw))#*base._flux(cp,cw)
sncosmo.write_griddata_ascii(cp,cw,np.array([-2.5*np.log10(cl) for p in cp]),'my_effect.dat')
sys.exit()
eff=sncosmo.TimeSeriesSource(cp,cw,np.array([cl for p in cp]))
plt.plot(cw,eff._flux(0,cw).flatten())
c=.9
plt.plot(cw,10**(-.4*c*sncosmo.Model('salt2-extended')._source._colorlaw(cw)))
plt.plot(cw,10**(-.4*(c))*eff._flux(0,cw).flatten())
#plt.plot(cw,cl)
#plt.plot(cw,base._flux(0,cw).flatten())
plt.show()
sys.exit()
'''
hsiao=sncosmo.Model('hsiao').source
hp=hsiao._phase
hw=hsiao._wave
hf=hsiao._flux(hp,hw)
plt.plot(hw,hf[0,:])
plt.show()
sys.exit()
stretched=sncosmo.StretchSource(hp,hw,hf)
stretched.set(s=1.1)

byosed.sed_to_effect(stretched,hsiao,rescale=False,effect_var='stretch',outname='gridded_stretch_hsiao.dat')

#comb=sncosmo.Model(sncosmo.TimeSeriesSource(hp,hw,np.array([hf[i]-stretched._flux(hp[i],hw).flatten() for i in range(len(hp))])))
#sncosmo.write_griddata_ascii(hp,hw,comb.flux(hp,hw),'gridded_hsiao_stretch.dat')
#plt.plot(hp,comb.bandflux('bessellb',hp))
#plt.plot(hp,sncosmo.Model(hsiao).bandflux('bessellb',hp))
#plt.show()
import os
hsiao=sncosmo.Model('hsiao').source
base_repo='/usr/local/WFIRST/ROOT/BYOSED_dev/BYOSEDINPUT'
data_repos=['SALT3_velocity_templates/SALT3_velocity_templates_werrors','SALT3_hostmass','ssfr_composites']
effect_keys=[['highv','lowv'],['lowmass','highmass'],['low','high']]
varname=['velocity','hostmass','ssfr']
eff_vals=[[[-30,-11.0001],[-11,0]],[[0,10.7],[10.7001,30]],[[-17,-10.7],[-10.6999,100]]]
for i in range(len(data_repos)):
    repo=os.path.join(base_repo,data_repos[i])
    seds=byosed.kaepora_to_sed(repo,
                      effect_keys[i],base_sed='hsiao',minWave=3500,maxWave=9000,minPhase=-np.inf,
               maxPhase=np.inf,waveStep=10,scale_band='bessellr')

    eff_dict={}
    for j in range(len(effect_keys[i])):
        for k in range(len(eff_vals[i][j])):
            eff_dict[eff_vals[i][j][k]]=seds[effect_keys[i][j]][2]
    byosed.sed_to_effect(eff_dict,hsiao,rescale=False,effect_var=varname[i],outname='byosed/initfiles/gridded_%s_hsiao.dat'%varname[i])

sys.exit()

salt2=sncosmo.Model(sncosmo.SALT2Source(modeldir='../fix_salt2_extended/'))
#hp,hw,hf=sncosmo.read_griddata_ascii('../salt3_testing/Hsi.dat')
#base=sncosmo.TimeSeriesSource(hp,hw,hf)
base=sncosmo.Model('hsiao').source
salt2.set_source_peakmag(sncosmo.Model(base).source_peakmag('bessellb','ab'),'bessellb','ab')
hp=base._phase
hw=base._wave

#
# salt2.set(x1=.89)
# peff=salt2._source._phase
weff=salt2._source._wave
# plt.plot(weff,ccm_extinction(weff,1)-ccm_extinction(sncosmo.get_bandpass('bessellb').wave_eff,1))
# plt.plot(weff,salt2._source._colorlaw(weff))
# plt.plot([sncosmo.get_bandpass('bessellb').wave_eff,sncosmo.get_bandpass('bessellb').wave_eff],(-.45,.25))
# plt.xlim((3500,20000))
# plt.ylim((-.45,.25))
# plt.show()
# sys.exit()
# feff=salt2._source._model['M1'](peff,weff)*salt2.get('x0')*salt2.get('x1')+base._flux(peff,weff)

sncosmo.write_griddata_ascii(hp,hw,np.array([ccm_extinction(hw,1)-ccm_extinction(sncosmo.get_bandpass('bessellb').wave_eff,1) for p in hp]),'../salt3_testing/gridded_color_hsiao.dat')
sys.exit()
# eff.set(amplitude=10**(-.4*(-19.365+.27)))
# base.set(amplitude=10**(-.4*(-19.365+.27)))
# salt2.set(x0=salt2.get('x0')*10**(-.4*(-19.365+.27)))
# print(np.sum(base._flux(10,weff).flatten()),np.sum(eff._flux(10,weff).flatten()))
#plt.plot(weff,eff._flux(10,weff).flatten())
#plt.plot(peff,sncosmo.Model(eff).bandflux('bessellb',peff))
#stretched=sncosmo.StretchSource(hp,hw,hf)
#stretched.set(s=1.1)

# comb=sncosmo.Model(sncosmo.TimeSeriesSource(hp,hw,np.array([hf[i]-stretched._flux(hp[i],hw).flatten() for i in range(len(hp))])))
# stretched.set(s=1.1)
# salt2.set(x1=1)
# plt.plot(peff,sncosmo.Model(base).bandflux('bessellb',peff))
salt2.set(x1=0)
salt2.set(c=1)
plt.plot(hp,sncosmo.Model(eff).bandflux('bessellr',hp))
plt.plot(hp,salt2.bandflux('bessellr',hp))
plt.show()
sys.exit()
# plt.plot(peff,stretched.bandflux('bessellb',peff))
#
#
# plt.show()
# sys.exit()

byosed.sed_to_effect(eff,base,rescale=True)
sys.exit()
mySED=byosed.GeneralSED(OPTMASK=1,ARGLIST=[],HOST_PARAM_NAMES='HOST_MASS,SFR,AGE,REDSHIFT,METALLICITY')
#byosed.plot_sed(mySED,effect='STRETCH')
flux = mySED.warp_SED(mySED.phase)
mod = mySED.to_sn_model()

print(mod.bandflux('bessellv',0))

#print(mySED.fetchParNames_BYOSED())
#mySED.fetchSED_BYOSED(0,5000,3,2,[2.5,1,1,.5])



#plt.plot(mySED.wave,mySED.sedInterp(0,mySED.wave)/mySED.x0)
#f=mySED.sedInterp(0,mySED.wave).flatten()/mySED.x0
#s=mySED.sn_effects['STRETCH'].flux(0*np.ones(len(mySED.wave)),mySED.wave,[],[])

#plt.plot(mySED.wave,f+s)
#plt.xlim(3400,10000)
#plt.show()
#sys.exit()
# effect='HOST_MASS'
# bounds=[5,20]
# leg_sym='logM'
effect='STRETCH'
leg_sym='$x_1$'
