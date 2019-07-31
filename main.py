
#func=lambda x,y,z:np.array(x+y+z)
#_generate_ND_grids(func,'test_host.dat',['phase','mass','velocity','theta'],np.array([0,20]),np.array([2,3]),np.array([4,5]))
#sys.exit()
#p,test=_read_ND_grids('salt2_m0.dat')
#print(test(np.array([[10,5000],[10,6000]])))
import matplotlib.pyplot as plt
import numpy as np
import byosed
#sys.exit()
mySED=byosed.GeneralSED('$WFIRST_ROOT/BYOSED_dev/BYOSEDINPUT/',2,[],'HOST_MASS,SFR,AGE,REDSHIFT,METALLICITY')

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
# effect='SFR'
# bounds=[.01,5]
# leg_sym='z'
# effect='METALLICITY'
# bounds=[.01,10]
# leg_sym='Z'
fig,ax=plt.subplots(nrows=3,ncols=3,figsize=(15,15),sharex=True)
phases=np.arange(-10,31,5)
k=0
base_params=[9,1,1,.001,1]
for i in range(3):
    for j in range(3):
        if effect in ['VELOCITY','STRETCH']:
            mySED.sn_effects[effect].scale_parameter=0
        else:
            mySED.host_effects[effect].warp_parameter=0
        ax[i][j].plot(mySED.wave,mySED.fetchSED_BYOSED(phases[k],5000,3,3,base_params),label='Hsiao',color='k',linewidth=2)

        for p in range(3):
            if effect not in ['VELOCITY','STRETCH']:
                mySED.host_effects[effect].updateWarp_Param()
                v=np.random.uniform(bounds[0],bounds[1])#mySED.sn_effects[effect].warp_parameter
                print(v)
                #mySED.sn_effects['STRETCH'].updateWarp_Param()
                #s=mySED.sn_effects['STRETCH'].warp_parameter
                if effect!='SFR':
                    ax[i][j].plot(mySED.wave,mySED.fetchSED_BYOSED(phases[k],5000,3,3,
                                                                   [base_params[i] if mySED.host_param_names[i]!=effect else v for i in range(len(base_params))]),label='%s=%.2f'%(leg_sym,v))
                else:
                    ax[i][j].plot(mySED.wave,mySED.fetchSED_BYOSED(phases[k],5000,3,3,
                                                                   [base_params[i] if mySED.host_param_names[i]!='REDSHIFT' else v for i in range(len(base_params))]),label='%s=%.2f'%(leg_sym,v))
            else:
                if effect=='VELOCITY':
                    mySED.sn_effects[effect].scale_parameter=1.
                    mySED.sn_effects[effect].updateWarp_Param()
                    v=mySED.sn_effects[effect].warp_parameter
                else:
                    mySED.sn_effects[effect].updateScale_Param()
                    v=mySED.sn_effects[effect].scale_parameter
                ax[i][j].plot(mySED.wave,mySED.fetchSED_BYOSED(phases[k],5000,3,3,base_params),label='%s=%.2f'%(leg_sym,v))

        ax[i][j].legend(fontsize=14)
        ax[i][j].annotate('Phase='+str(phases[k]),(.5,.05),fontsize=14,xycoords='axes fraction')
        ax[i][j].set_xlim((3000,9500))
        k+=1
        if j==0:
            ax[i][j].set_ylabel('Flux',fontsize=16)
        if i==2 and j==1:
            ax[i][j].set_xlabel('Wavelength ($\AA$)',fontsize=16)
        ax[i][j].tick_params(axis='x', labelsize=14)
        ax[i][j].tick_params(axis='y', labelsize=14)

plt.show()#savefig('/Users/jpierel/rodney/salt3_testing/'+effect+'_byosed.pdf',format='pdf')





