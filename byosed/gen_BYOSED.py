import numpy as np
import matplotlib.pyplot as plt

import os, optparse, configparser, sys, sncosmo
from scipy.interpolate import interp2d
from copy import copy

from .effectio import *
from .distributions import *
from .effects import *

if not hasattr(sys, 'argv'):
		sys.argv  = ['']

required_keys = ['magsmear']


__mask_bit_locations__={'verbose':1,'dump':2}
__all__ = ['GeneralSED']


class GeneralSED:

	def __init__(self,PATH_VERSION=os.path.join(os.path.dirname(os.path.abspath(__file__)),'initfiles','BYOSED.params'),
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
		self.norm =float(config['MAIN']['NORM']) if 'NORM' in config['MAIN'].keys() else -19.365
		self.alpha =float(config['MAIN']['ALPHA']) if 'ALPHA' in config['MAIN'].keys() else 0.14
		self.beta =float(config['MAIN']['BETA']) if 'BETA' in config['MAIN'].keys() else 3.1

		#self.x0=10**(.4*19.365)
		self.x0=10**(-.4*self.norm)

		if self.options.magsmear!=0.0:
			self.magsmear = np.random.normal(0,self.options.magsmear)
		else:
			self.magsmear = 0.0
		self.magoff=self.options.magoff

		self.flux = fluxarr*self.x0*10**(-0.4*self.magoff)

		self.phase = np.unique(phase)
		self.wave = np.unique(wave)
		self.wavelen = len(self.wave)

		self.sedInterp=interp2d(self.phase,self.wave,self.flux.T,kind='linear',bounds_error=True)



		#fluxsmear=self.fetchSED_BYOSED([0],5000,1,1,[0 for x in range(len(self.host_param_names))])
		#import pickle
		#with open('/project2/rkessler/SURVEYS/WFIRST/USERS/jpierel/byosed/fluxsmear.dat','wb') as f:
		#	pickle.dump([self.wave,fluxsmear],f)


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
			if 'SCALE_TYPE' not in warp_data.keys():
				warp_data['SCALE_TYPE']='inner'
			elif warp_data['SCALE_TYPE'] not in ['inner','outer']:
				raise RuntimeError("Do not recognize variable SCALE_TYPE, should be 'inner' or 'outer'")



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

				sn_dict[warp]=WarpModel(warp_function=sn_function,
										param_names=sn_param_names,
										parameters=np.array([0. if sn_param_names[i]!=warp.upper() else warp_parameter for i in range(len(sn_param_names))]),
										warp_parameter=warp_parameter,
										warp_distribution=distribution['PARAM'] if 'PARAM' in distribution.keys() else None,
										scale_parameter=sn_scale_parameter,
										scale_distribution=distribution['SCALE'],
										scale_type=warp_data['SCALE_TYPE'],
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
				host_dict[warp]=WarpModel(warp_function=host_function,
										  param_names=host_param_names,
										  parameters=np.array([0. if host_param_names[i]!=warp.upper() else warp_parameter for i in range(len(host_param_names))]),

										  warp_parameter=warp_parameter,
										  warp_distribution=distribution['PARAM'] if 'PARAM' in distribution.keys() else None,
										  scale_parameter=host_scale_parameter,
										  scale_distribution=distribution['SCALE'],
										  scale_type=warp_data['SCALE_TYPE'],
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

	def fetchSED_BYOSED(self,trest,maxlam=5000,external_id=1,new_event=1,hostpars=''):

		try:
			if len(self.wave)>maxlam:
				raise RuntimeError("Your wavelength array cannot be larger than %i but is %i"%(maxlam,len(self.wave)))
			#iPhase = np.where(np.abs(trest-self.phase) == np.min(np.abs(trest-self.phase)))[0][0]
			#print('HOST_PARAMS: ',hostpars)
			if self.sn_id is None:
				self.sn_id=external_id
				newSN=False

			elif external_id!=self.sn_id:
				newSN=True

			else:
				newSN=False

			self.sn_id=external_id
			fluxsmear=self.sedInterp(trest,self.wave).flatten()
			orig_fluxsmear=copy(fluxsmear)

			if self.options.magsmear!=0.0 and (self.sn_id!=external_id or self.magsmear is None):
				self.magsmear=np.random.normal(0,self.options.magsmear)
			else:
				self.magsmear=0.0
			#if self.sn_id!=external_id:
			#	fluxsmear *= 10**(-0.4*self.magoff)
			fluxsmear *= 10**(-0.4*(self.magsmear))

			trest_arr=trest*np.ones(len(self.wave))

			inner_product=np.zeros(len(self.wave))
		except Exception as e:
			print('Python Error :',e)


		outer_product=np.zeros(len(self.wave))

		for warp in [x for x in self.warp_effects]:# if x!='COLOR']:

			try: #if True:

				if newSN:

					if warp in self.sn_effects.keys():
						self.sn_effects[warp].updateWarp_Param()
						self.sn_effects[warp].updateScale_Param()
						if warp in self.host_effects.keys():
							self.host_effects[warp].updateWarp_Param()
							self.host_effects[warp].scale_parameter=1.

					else:
						self.host_effects[warp].updateWarp_Param()
						self.host_effects[warp].updateScale_Param()

				#try:

				#	existing=np.loadtxt('color.dat').reshape(-1,2)
				#	existing=np.append(existing,[hostpars[(self.host_param_names).index('ZCMB')],
				#   self.sn_effects[warp].scale_parameter]).reshape(-1,2)

				#	np.savetxt('color.dat',existing)
				#except:
				#	np.savetxt('color.dat',np.array([hostpars[(self.host_param_names).index('ZCMB')],
				#	 self.sn_effects[warp].scale_parameter]).reshape(-1,2))

				# not sure about the multiplication by x0 here, depends on if SNANA is messing with the
				# absolute magnitude somewhere else
				product=np.ones(len(self.wave))
				temp_scale_param = 0
				temp_outer_product=np.ones(len(self.wave))
				outer_scale_param=0

				if warp in self.sn_effects.keys():
					if self.verbose:
						if self.sn_effects[warp].warp_parameter is not None:
							print('Phase=%.1f, %s: %.2f'%(trest,warp,self.sn_effects[warp].warp_parameter))
						else:
							print('Phase=%.1f, %s: %.2f'%(trest,warp,self.sn_effects[warp].scale_parameter))

					if self.sn_effects[warp].scale_type=='inner':

						product*=self.sn_effects[warp].flux(trest_arr,self.wave,hostpars,self.host_param_names)
						if temp_scale_param ==0:
							temp_scale_param = self.sn_effects[warp].scale_parameter
						else:
							temp_scale_param*=self.sn_effects[warp].scale_parameter

					else:

						temp_outer_product*=self.sn_effects[warp].flux(trest_arr,self.wave,hostpars,self.host_param_names)
						if outer_scale_param ==0:
							outer_scale_param = self.sn_effects[warp].scale_parameter
						else:
							outer_scale_param*=self.sn_effects[warp].scale_parameter

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
					if self.host_effects[warp].scale_type=='inner':
						if temp_scale_param==0:
							temp_scale_param=self.host_effects[warp].scale_parameter
						else:
							temp_scale_param*=self.host_effects[warp].scale_parameter
						product*=self.host_effects[warp].flux(trest_arr,self.wave,hostpars,self.host_param_names)
						
					else:
						temp_outer_product*=self.host_effects[warp].flux(trest_arr,self.wave,hostpars,self.host_param_names)
						if outer_scale_param==0:
							outer_scale_param=self.host_effects[warp].scale_parameter
						else:
							outer_scale_param*=self.host_effects[warp].scale_parameter
						outer_scale_param*=self.host_effects[warp].scale_parameter
	
				inner_product+=product*temp_scale_param

				outer_product+=temp_outer_product*outer_scale_param

			except Exception as e:
				print('Python Error :',e)


		
		fluxsmear*=((1+inner_product)*10**(-0.4*outer_product))

		fluxsmear*=self.brightness_correct_Ia()

		return list(fluxsmear)
		
	def brightness_correct_Ia(self):
		if 'COLOR' in self.sn_effects.keys():
			c=self.sn_effects['COLOR'].scale_parameter
		else:
			c=0
		if 'STRETCH' in self.sn_effects.keys():
			s=self.sn_effects['STRETCH'].scale_parameter
		else:
			s=0

		return(10**(-.4*(self.beta*c-self.alpha*s)))

	def warp_SED(self,trest=None,hostpars=None):
		if trest is None:
			trest=self.phase
		if hostpars is not None:
			self.host_param_names = ','.join(list(hostpars.keys()))
			hostpar_vals = ','.join(list(hostpars.values()))
		else:
			hostpar_vals = ''
		warped_flux = []

		for p in trest:
			warped_flux=np.append(warped_flux,self.fetchSED_BYOSED(p,hostpars=hostpar_vals))

		self.warped_phase = trest
		self.warped_wave = self.wave
		self.warped_flux = warped_flux.reshape([len(self.warped_phase),len(self.warped_wave)])

		self.warpedSED_Interp=interp2d(self.warped_phase,self.warped_wave,self.warped_flux.T,kind='linear',bounds_error=True)

		return (self.warped_flux)


	def to_sn_model(self):
		try:
			source = sncosmo.TimeSeriesSource(self.warped_phase,self.warped_wave,self.warped_flux)
		except:
			print("Could not create sncosmo model, have you created warped data yet?")
		self.sn_model = sncosmo.Model(source)
		return(self.sn_model)

	def plot_sed(self,effect=None,host_param_bounds=None):
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

					if self.sn_effects[effect].warp_parameter is None:
						self.sn_effects[effect].scale_parameter=0
					else:
						self.host_effects[effect].warp_parameter=0
					ax[i][j].plot(self.wave,self.fetchSED_BYOSED(phases[k],5000,3,3,base_params),label='Baseline',color='k',linewidth=2)

					for p in range(3):
						if effect in self.host_effects.keys():
							if host_param_bounds is None:
								raise RuntimeError("Please provide bounds for Host parameter, used in effect %s"%effect)
							self.host_effects[effect].updateWarp_Param()
							v=np.random.uniform(host_param_bounds[0],host_param_bounds[1])#gen_sed.sn_effects[effect].warp_parameter
							#gen_sed.sn_effects['STRETCH'].updateWarp_Param()
							#s=gen_sed.sn_effects['STRETCH'].warp_parameter
							ax[i][j].plot(self.wave,self.fetchSED_BYOSED(phases[k],5000,3,3,
																			   [base_params[i] if self.host_param_names[i]!=effect else v for i in range(len(base_params))]),label='%s=%.2f'%(leg_sym,v))

						else:
							if self.sn_effects[effect].warp_parameter is not None:
								self.sn_effects[effect].scale_parameter=1.
								self.sn_effects[effect].updateWarp_Param()
								v=self.sn_effects[effect].warp_parameter
							else:
								self.sn_effects[effect].updateScale_Param()
								v=self.sn_effects[effect].scale_parameter
							ax[i][j].plot(self.wave,self.fetchSED_BYOSED(phases[k],5000,3,3,base_params),label='%s=%.2f'%(leg_sym,v))

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










	
	
	
	

