import numpy as np


__all__ = ['WarpModel']



class WarpModel(object):
    """Base class for anything with parameters.

    Derived classes must have properties ``_param_names`` (list of str)
    and ``_parameters`` (1-d numpy.ndarray).
    """

    def __init__(self,warp_function,
                    parameters,
                    param_names,
                    warp_parameter,
                    warp_distribution,
                    scale_parameter,
                    scale_distribution,
                    name):


        self.name = name
        self._parameters = parameters
        self._param_names = [x.upper() for x in param_names]
        self.warp_function=warp_function
        self.warp_parameter=warp_parameter
        self.scale_parameter=scale_parameter
        self.warp_distribution=warp_distribution
        self.scale_distribution=scale_distribution

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