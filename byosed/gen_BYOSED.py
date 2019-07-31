import numpy as np	
import os, optparse, configparser, sys
from scipy.interpolate import interp2d
from copy import copy

from .effectio import *
from .distributions import *

if not hasattr(sys, 'argv'):
		sys.argv  = ['']

required_keys = ['magsmear']


__mask_bit_locations__={'verbose':1,'dump':2}
__all__ = ['GeneralSED']


class GeneralSED:

	def __init__(self,PATH_VERSION=os.path.join(os.path.abspath(__file__),'initfiles'),
				 OPTMASK=2,ARGLIST='',HOST_PARAM_NAMES=''):
		#print('LIST: ',OPTMASK)
		#print('HOST_PARAM_NAMES: ',HOST_PARAM_NAMES)
		# TODO: write a print statement that warns if
		# HOST_PARAM_NAMES is a variable that the code
		# isn't going to do anything with

		self.verbose = OPTMASK & (1 << __mask_bit_locations__['verbose']) > 0

		self.PATH_VERSION = os.path.expandvars(os.path.dirname(PATH_VERSION))

		self.host_param_names = [x.upper() for x in HOST_PARAM_NAMES.split(',')]
		self.PATH_VERSION = os.path.dirname(PATH_VERSION)

		self.dump = OPTMASK & (1 << __mask_bit_locations__['dump'])>0
		self.sn_id=None

		self.PATH_VERSION = os.path.expandvars(os.path.dirname(PATH_VERSION))
		#self.host_param_names=HOST_PARAM_NAMES

		self.paramfile = os.path.join(self.PATH_VERSION,'BYOSED.params')
		if os.path.exists(self.paramfile):
			config = configparser.ConfigParser()
			config.read(self.paramfile)
		else: raise RuntimeError('param file %s not found!'%self.paramfile)

		parser = self.add_options(usage='',config=config)

		options,  args = parser.parse_args()

		for k in required_keys:
			if k not in options.__dict__.keys():
				raise RuntimeError('key %s not in parameter file'%k)
		self.options = options

		self.warp_effects=self.fetchParNames_CONFIG(config)

		self.sn_effects,self.host_effects=self.fetchWarp_BYOSED(config)

		phase,wave,flux = np.loadtxt(os.path.join(self.PATH_VERSION,self.options.sed_file),unpack=True)


		fluxarr = flux.reshape([len(np.unique(phase)),len(np.unique(wave))])
		self.x0=10**(.4*19.365)
		self.flux = fluxarr*self.x0

		self.phase = np.unique(phase)
		self.wave = np.unique(wave)
		self.wavelen = len(self.wave)

		self.sedInterp=interp2d(self.phase,self.wave,self.flux.T,kind='linear',bounds_error=True)
		print(self.warp_effects)

		return

	def add_options(self, parser=None, usage=None, config=None):

		if parser == None:
			parser = optparse.OptionParser(usage=usage, conflict_handler="resolve")
		# The basics
		parser.add_option('-v', '--verbose', action="count", dest="verbose",default=1)
		parser.add_option('--clobber', default=False, action="store_true",help='clobber output file')
		parser.add_option('--magsmear', default=config.get('MAIN','MAGSMEAR'),
						  type="float",help='amount of Gaussian-random mag smearing (default=%default)')
		parser.add_option('--magoff', default=config.get('MAIN','MAGOFF'),
						  type="float",help='mag offset (default=%default)')
		parser.add_option('--sed_file',default=config.get('MAIN','SED_FILE'),
						  type='str',help='Name of sed file')


		return parser



	def fetchWarp_BYOSED(self,config):
		#read in warp effect data
		#if speed becomes a real issue, we could
		#combine the various additive functions
		#into a single function, which would have
		#to be repeated for every SN but would
		#save looping through all for various
		#phases...not sure if that would actually
		#come out on top though.

		sn_dict=dict([])
		host_dict=dict([])


		for warp in self.warp_effects:
			warp_data={}
			for k in config[warp]:
				try:
					warp_data[k.upper()]=np.array(config.get(warp,k).split()).astype(float)
				except:
					warp_data[k.upper()]=config.get(warp,k)




			if 'SN_FUNCTION' in warp_data:
				if 'DIST' in ' '.join([x for x in warp_data.keys() if 'SN' in x or 'SCALE' in x]):
					distribution=get_distribution(warp,{k:warp_data[k] for k in warp_data.keys() if 'SN' in k or 'SCALE' in k},self.PATH_VERSION,'SN')
				else:
					raise RuntimeError("Did not supply scale distribution information for SN effect %s."%warp)


				if 'SN_FUNCTION_SCALE' not in warp_data:
					scale_factor=1.
				else:
					scale_factor=warp_data['SN_FUNCTION_SCALE']
				try:
					sn_param_names,sn_function=read_ND_grids(os.path.expandvars(os.path.join(self.PATH_VERSION,str(warp_data['SN_FUNCTION']))),scale_factor)
				except RuntimeError:
					raise RuntimeError("Do not recognize format of function for %s SN Function"%warp)
				if warp.upper() in sn_param_names and 'PARAM' not in distribution.keys():
					raise RuntimeError("Must supply parameter distribution for SN effect %s"%warp)
				sn_scale_parameter=distribution['SCALE']()[0]
				warp_parameter=distribution['PARAM']()[0] if 'PARAM' in distribution.keys() else None
				if warp.upper() in sn_param_names and warp_parameter is None:
					raise RuntimeError("Woops, you are not providing a PARAM distribution for your %s effect."%warp.upper())

				sn_dict[warp]=warpModel(warp_function=sn_function,
										param_names=sn_param_names,
										parameters=np.array([0. if sn_param_names[i]!=warp.upper() else warp_parameter for i in range(len(sn_param_names))]),
										warp_parameter=warp_parameter,
										warp_distribution=distribution['PARAM'] if 'PARAM' in distribution.keys() else None,
										scale_parameter=sn_scale_parameter,
										scale_distribution=distribution['SCALE'],
										name=warp)

			if 'HOST_FUNCTION' in warp_data:
				if 'DIST' in ' '.join([x for x in warp_data.keys() if 'HOST' in x or 'SCALE' in x]):
					distribution=get_distribution(warp,{k:warp_data[k] for k in warp_data.keys() if 'HOST' in k or 'SCALE' in k},self.PATH_VERSION,'HOST')
				else:
					raise RuntimeError("Did not supply scale distribution information for HOST effect %s."%warp)
				try:
					host_param_names,host_function=read_ND_grids(os.path.expandvars(os.path.join(self.PATH_VERSION,str(warp_data['HOST_FUNCTION']))))
				except RuntimeError:
					raise RuntimeError("Do not recognize format of function for %s HOST Function"%warp)

				if warp.upper() in host_param_names and 'PARAM' not in distribution.keys() and warp.upper() not in self.host_param_names:
					raise RuntimeError("Must supply parameter distribution for HOST effect %s"%warp)

				host_scale_parameter=distribution['SCALE']()[0]
				warp_parameter=distribution['PARAM']()[0] if 'PARAM' in distribution.keys() else None
				if warp.upper() in host_param_names and warp_parameter is None and warp.upper() not in self.host_param_names:
					raise RuntimeError("Woops, you are not providing a PARAM distribution for your %s effect."%warp.upper())
				host_dict[warp]=warpModel(warp_function=host_function,
										  param_names=host_param_names,
										  parameters=np.array([0. if host_param_names[i]!=warp.upper() else warp_parameter for i in range(len(host_param_names))]),

										  warp_parameter=warp_parameter,
										  warp_distribution=distribution['PARAM'] if 'PARAM' in distribution.keys() else None,
										  scale_parameter=host_scale_parameter,
										  scale_distribution=distribution['SCALE'],
										  name=warp)

		return(sn_dict,host_dict)


	#def updateWarping_Params(self):
	#		for warp in self.warp_effects:
	#			if warp in sn_effects.keys():
	#				self.sn_effects[warp].update(warp_parameter
	#		return self


	def fetchSED_NLAM(self):
		return self.wavelen

	def fetchSED_LAM(self):
		return list(self.wave)


	def fetchSED_BYOSED(self,trest,maxlam,external_id,new_event,hostpars):
		if len(self.wave)>maxlam:
			raise RuntimeError("Your wavelength array cannot be larger than %i but is %i"%(maxlam,len(self.wave)))
		#iPhase = np.where(np.abs(trest-self.phase) == np.min(np.abs(trest-self.phase)))[0][0]
		#print('HOST_PARAMS: ',hostpars)
		if self.sn_id is None:
			self.sn_id=external_id
		fluxsmear=self.sedInterp(trest,self.wave).flatten()
		orig_fluxsmear=copy(fluxsmear)
		if self.options.magsmear!=0.0:
			fluxsmear *= 10**(0.4*(np.random.normal(0,self.options.magsmear)))
		trest_arr=trest*np.ones(len(self.wave))

		for warp in [x for x in self.warp_effects]:# if x!='COLOR']:
			try: #if True:

				if external_id!=self.sn_id:
					if warp in self.sn_effects.keys():
						self.sn_effects[warp].updateWarp_Param()
						self.sn_effects[warp].updateScale_Param()
						if warp in self.host_effects.keys():
							self.host_effects[warp].updateWarp_Param()
							self.host_effects[warp].scale_parameter=1.
					else:
						self.host_effects[warp].updateWarp_Param()
						self.host_effects[warp].updateScale_Param()
					self.sn_id=external_id


				# not sure about the multiplication by x0 here, depends on if SNANA is messing with the
				# absolute magnitude somewhere else
				product=0.
				temp_scale_param = 1.
				if warp in self.sn_effects.keys():
					if self.verbose:
						if self.sn_effects[warp].warp_parameter is not None:
							print('Phase=%.1f, %s: %.2f'%(trest,warp,self.sn_effects[warp].warp_parameter))
						else:
							print('Phase=%.1f, %s: %.2f'%(trest,warp,self.sn_effects[warp].scale_parameter))
					product+=self.sn_effects[warp].flux(trest_arr,self.wave,hostpars,self.host_param_names)

					temp_scale_param*=self.sn_effects[warp].scale_parameter
				#if warp in self.sn_effects[warp]._param_names:
				#	temp_warp_param=1.
				#else:
				#	temp_warp_param=self.sn_effects[warp].warp_parameter
				if warp in self.host_effects.keys():
					if self.verbose:
						if self.host_effects[warp].warp_parameter is not None:
							print('Phase=%.1f, %s: %.2f'%(trest,warp,self.host_effects[warp].warp_parameter))
						else:
							print('Phase=%.1f, %s: %.2f'%(trest,warp,self.host_effects[warp].scale_parameter))

					product+=self.host_effects[warp].flux(trest_arr,self.wave,hostpars,self.host_param_names)
					temp_scale_param*=self.host_effects[warp].scale_parameter

				fluxsmear*=(1.+temp_scale_param*product)
			except RuntimeError:
				import pdb; pdb.set_trace()


		# if 'COLOR' in self.warp_effects:
		# 		if external_id!=self.sn_id:
		# 				self.sn_effects['COLOR'].updateWarping_Params()
		# 				self.sn_id=external_id
		# 		if self.verbose:
		# 				print('Phase=%.1f, %s: %.2f'%(trest,'COLOR',self.sn_effects['COLOR'].warp_parameter))
		# 		fluxsmear*=10**(-0.4*self.sn_effects['COLOR'].warp_parameter*\
		# 				self.sn_effects['COLOR'].flux(trest_arr,self.wave,hostpars,self.host_param_names).flatten())

		return list(fluxsmear)

	def fetchParNames_BYOSED(self):
		return list(self.warp_effects)

	def fetchNParNames_BYOSED(self):
		return len(self.warp_effects)

	def fetchParVals_BYOSED_4SNANA(self,varname):
		if varname in self.sn_effects.keys():
			if self.sn_effects[varname].warp_parameter is not None:
				return self.sn_effects[varname].warp_parameter
			else:
				return self.sn_effects[varname].scale_parameter
		else:
			if self.host_effects[varname].warp_parameter is not None:
				return self.host_effects[varname].warp_parameter
			else:
				return self.host_effects[varname].scale_parameter


	def fetchParVals_BYOSED(self,varname):
		if varname in self.sn_effects.keys():
			return self.sn_effects[varname].warp_parameter
		else:
			return self.host_effects[varname].warp_parameter
	def fetchParNames_CONFIG(self,config):
		if 'FLAGS' in config.sections():
			return([k.upper() for k in list(config['FLAGS'].keys()) if config['FLAGS'][k]=='True'])
		else:
			return([x for x in config.sections() if x not in ['MAIN','FLAGS']])


class warpModel(object):
	"""Base class for anything with parameters.

	Derived classes must have properties ``_param_names`` (list of str)
	and ``_parameters`` (1-d numpy.ndarray).
	"""

	def __init__(self, warp_function,parameters,param_names,warp_parameter,warp_distribution,scale_parameter,scale_distribution,name):
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







	
	
	
	

