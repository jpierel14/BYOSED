import numpy as np
import sncosmo

from .effectio import generate_ND_grids,line_prepender

__all__ = ['WarpModel','sed_to_effect','residual_to_effect']

def sed_to_effect(effect_in,base_sed,effect_var='Theta',scale_band='bessellb',rescale=True,outname='my_effect.dat',
                  phase_wave_source='effect',diff_limit=100):
    if isinstance(effect_in,dict):
        effect_dict=effect_in
        out_effect_dict={e:{} for e in effect_in.keys()}
        single=False
    else:
        effect_dict={'single':effect_in}
        single=True

    base_phase=base_sed._phase
    base_wave=base_sed._wave

    overall_min_phase=-np.inf
    overall_max_phase=np.inf
    overall_max_wave=np.inf
    overall_min_wave=-np.inf
    for e in effect_dict.keys():
        effect_sed=effect_dict[e]
        effect_phase=effect_sed._phase
        effect_wave=effect_sed._wave
        if np.min(effect_wave)>overall_min_wave:
            overall_min_wave=np.min(effect_wave)
        if np.max(effect_wave)<overall_max_wave:
            overall_max_wave=np.max(effect_wave)
        if np.min(effect_phase)>overall_min_phase:
            overall_min_phase=np.min(effect_phase)
        if np.max(effect_phase)<overall_max_phase:
            overall_max_phase=np.max(effect_phase)

        #if phase_wave_source=='effect':
        #    base_phase=effect_sed._phase
        #    base_wave=effect_sed._wave

        fluxes=[]
        for phase in base_phase:
            base_flux=base_sed._flux(phase,base_wave).flatten()
            if phase>=np.min(effect_phase) and phase<=np.max(effect_phase):
                scale_factor=effect_sed.bandflux(scale_band,phase)/base_sed.bandflux(scale_band,phase)

                try:
                    effect_flux=effect_sed._flux(phase,base_wave).flatten()
                except:
                    effect_flux=np.zeros(len(base_wave))
                    inds=base_wave[np.where(np.logical_and(base_wave>=np.min(effect_wave),
                                                             base_wave<=np.max(effect_wave)))[0]]
                    effect_flux[inds]=effect_sed._flux(phase,base_wave[inds])
                if not rescale:
                    scale_factor=1.
                effect_flux/=scale_factor
                effect_flux=(effect_flux-base_flux)/base_flux
                effect_flux[np.abs(effect_flux)==np.inf]=0.
                effect_flux[np.isnan(effect_flux)]=0.
                effect_flux[np.abs(effect_flux)>diff_limit]=0.
            else:
                effect_wave=base_wave
                effect_flux=np.zeros(len(base_wave))

            if not single:
                out_effect_dict[e][phase]={base_wave[i]:effect_flux[i] for i in range(len(base_wave))}
            else:
                fluxes.append(effect_flux)

    accept_phase_inds=np.where(np.logical_and(base_phase>=overall_min_phase,base_phase<=overall_max_phase))[0]
    accept_wave_inds=np.where(np.logical_and(base_wave>=overall_min_wave,base_wave<=overall_max_wave))[0]

    if single:

        fluxes=np.array(fluxes).reshape((len(base_phase),len(base_wave)))

        base_phase=base_phase[accept_phase_inds].astype(int)
        base_wave=base_wave[accept_wave_inds].astype(int)
        fluxes=fluxes[np.ix_(accept_phase_inds,accept_wave_inds)]
        sncosmo.write_griddata_ascii(base_phase,base_wave,fluxes,outname)
        line_prepender(outname,'# phase wavelength flux')

    else:
        base_phase=base_phase[accept_phase_inds].astype(int)
        base_wave=base_wave[accept_wave_inds].astype(int)
        def func(velocity,phase,wavelength):
            flux=np.array([out_effect_dict[v][p][w] for v,p,w in zip(velocity,phase,wavelength)])
            return flux
        generate_ND_grids(func,outname,[effect_var,'phase','wavelength','flux'],
                          np.sort(list(out_effect_dict.keys())),base_phase,base_wave)
        #line_prepender(outname,'# %s phase wavelength flux'%effect_var)

def residual_to_effect():
    pass

class WarpModel(object):
    """Base class for anything with parameters.
    Derived classes must have properties ``_param_names`` (list of str)
    and ``_parameters`` (1-d numpy.ndarray).
    """

    def __init__(self, warp_function,parameters,param_names,warp_parameter,warp_distribution,
                 scale_parameter,scale_distribution,scale_type,name):
        self.name = name
        self._parameters = parameters
        self._param_names = [x.upper() for x in param_names]
        self.warp_function=warp_function
        self.warp_parameter=warp_parameter
        self.scale_parameter=scale_parameter
        self.warp_distribution=warp_distribution
        self.scale_distribution=scale_distribution
        self.scale_type=scale_type

    def setWarp_Param(self,param_name,param_value):
        self.set(**{param_name:param_value})

    def updateWarp_Param(self):
        if self.warp_distribution is not None:
            self.warp_parameter=self.warp_distribution()[0]
            if self.name in self._param_names:
                self.set(**{self.name:self.warp_parameter})

    #else:
    #	print("Cannot update warping param, no distribution.")

    def updateScale_Param(self):
        self.scale_parameter=self.scale_distribution()[0]


    def flux(self,phase,wave,host_params,host_param_names):
        phase_wave_dict={'PHASE':phase,'WAVELENGTH':wave}
        self.set(**{p:host_params[host_param_names.index(p)] for p in self._param_names if p in host_param_names})
        parameter_arrays=[np.ones(len(wave))*self._parameters[i] if self._param_names[i] not in ['PHASE','WAVELENGTH']
                          else phase_wave_dict[self._param_names[i]] for i in range(len(self._param_names))]
        return(self.warp_function(np.vstack(parameter_arrays).T).flatten())

    @property
    def param_names(self):
        """List of parameter names."""
        return self._param_names

    @property
    def parameters(self):
        """Parameter value array"""
        return self._parameters

    @parameters.setter
    def parameters(self, value):
        value = np.asarray(value)
        if value.shape != self._parameters.shape:
            raise ValueError("Incorrect number of parameters.")
        self._parameters[:] = value

    def set(self, **param_dict):
        """Set parameters of the model by name."""
        self.update(param_dict)

    def update(self, param_dict):
        """Set parameters of the model from a dictionary."""
        for key, value in param_dict.items():
            self[key] = value

    def __setitem__(self, key, value):
        """Set a single parameter of the model by name."""
        try:
            i = self._param_names.index(key)
        except ValueError:
            raise KeyError("Unknown parameter: " + repr(key))
        self._parameters[i] = value

    def get(self, name):
        """Get parameter of the model by name."""
        return self[name]

    def __getitem__(self, name):
        """Get parameter of the model by name"""
        try:
            i = self._param_names.index(name)
        except ValueError:
            raise KeyError("Model has no parameter " + repr(name))
        return self._parameters[i]


    def __str__(self):
        parameter_lines = [self._headsummary(), 'parameters:']
        if len(self._param_names) > 0:
            m = max(map(len, self._param_names))
            extralines = ['	 ' + k.ljust(m) + ' = ' + repr(v)
                          for k, v in zip(self._param_names, self._parameters)]
            parameter_lines.extend(extralines)
        return '\n'.join(parameter_lines)

    def __copy__(self):
        """Like a normal shallow copy, but makes an actual copy of the
        parameter array."""
        new_model = self.__new__(self.__class__)
        for key, val in self.__dict__.items():
            new_model.__dict__[key] = val
        new_model._parameters = self._parameters.copy()
        return new_model