
#func=lambda x,y,z:np.array(x+y+z)
#_generate_ND_grids(func,'test_host.dat',['phase','mass','velocity','theta'],np.array([0,20]),np.array([2,3]),np.array([4,5]))
#sys.exit()
#p,test=_read_ND_grids('salt2_m0.dat')
#print(test(np.array([[10,5000],[10,6000]])))
import matplotlib.pyplot as plt
import numpy as np
import byosed

#sys.exit()
mySED=byosed.GeneralSED('$WFIRST_ROOT/BYOSED_dev/BYOSEDINPUT/',1,[],'HOST_MASS,SFR,AGE,REDSHIFT,METALLICITY')
#byosed.plot_sed(mySED,effect='STRETCH')
mySED.warp_SED(mySED.phase)
mod = mySED.create_sn_model(mySED)

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
