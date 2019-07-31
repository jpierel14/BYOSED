import matplotlib.pyplot as plt
import numpy as np


def plot_sed(gen_sed,effect=None,host_param_bounds=None):
    
    #effect='STRETCH'
    leg_sym='$x_1$'
    # effect='SFR'
    # param_bounds=[.01,5]
    # leg_sym='z'
    # effect='METALLICITY'
    # param_bounds=[.01,10]
    # leg_sym='Z'
    fig,ax=plt.subplots(nrows=3,ncols=3,figsize=(15,15),sharex=True)
    phases=np.arange(-10,31,5)
    k=0

    base_params=[9,1,1,.001,1]
    for i in range(3):
        for j in range(3):
            if effect is not None:
                
                if gen_sed.sn_effects[effect].warp_parameter is None:
                    gen_sed.sn_effects[effect].scale_parameter=0
                else:
                    gen_sed.host_effects[effect].warp_parameter=0
                ax[i][j].plot(gen_sed.wave,gen_sed.fetchSED_BYOSED(phases[k],5000,3,3,base_params),label='Baseline',color='k',linewidth=2)
        
                for p in range(3):
                    if effect in gen_sed.host_effects.keys():
                        if host_param_bounds is None:
                            raise RuntimeError("Please provide bounds for Host parameter, used in effect %s"%effect)
                        gen_sed.host_effects[effect].updateWarp_Param()
                        v=np.random.uniform(host_param_bounds[0],host_param_bounds[1])#gen_sed.sn_effects[effect].warp_parameter
                        #gen_sed.sn_effects['STRETCH'].updateWarp_Param()
                        #s=gen_sed.sn_effects['STRETCH'].warp_parameter
                        ax[i][j].plot(gen_sed.wave,gen_sed.fetchSED_BYOSED(phases[k],5000,3,3,
                                                                           [base_params[i] if gen_sed.host_param_names[i]!=effect else v for i in range(len(base_params))]),label='%s=%.2f'%(leg_sym,v))

                    else:
                        if gen_sed.sn_effects[effect].warp_parameter is not None:
                            gen_sed.sn_effects[effect].scale_parameter=1.
                            gen_sed.sn_effects[effect].updateWarp_Param()
                            v=gen_sed.sn_effects[effect].warp_parameter
                        else:
                            gen_sed.sn_effects[effect].updateScale_Param()
                            v=gen_sed.sn_effects[effect].scale_parameter
                        ax[i][j].plot(gen_sed.wave,gen_sed.fetchSED_BYOSED(phases[k],5000,3,3,base_params),label='%s=%.2f'%(leg_sym,v))
    
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