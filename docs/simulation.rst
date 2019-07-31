***********************
Simulating with BYOSED
***********************


The BYOSED framework allows any spectrophotometric model to be used
as the underlying template to generate simulated Type Ia light curve data
with SNANA. By default, this model is the Hsiao+07 model (initfiles/Hsiao07.dat).
This can be replaced by any model.


Param File Basics
=================

The only file to set up is the BYOSED.params file. This contains the general aspects
of the simulated SN you want to create using BYOSED, and any warping effects you
want to add in. This file is separated into the following required and optional sections:

[MAIN]
------
**(Required)**

This section contains **SED_FILE** (name of SED file), as well as **MAGSMEAR** (magnitude 
smearing) and **MAGOFF** (magnitude offsets) definitions to be applied to the base SED defined by
sed_file. You may also define **CLOBBER** and **VERBOSE** flags here as well. This section may look
like the following:

::
	
	[MAIN]

	SED_FILE: Hsiao07.dat
	MAGSMEAR: 0.0
	MAGOFF: 0.0


[FLAGS]
-------
**(Optional)**

This section allows you to simply turn warping effects defined in the next section(s) on and off. If
this section exists, then it supersedes later sections and defines the warping effects to be used. 
If it does not exist, all defined warping effects are used. Adding this onto the **[MAIN]** section,
the params file might now look like the following:

::

	[MAIN]

	SED_FILE: Hsiao07.dat
	MAGSMEAR: 0.1
	MAGOFF: 0.0

	[FLAGS]

	COLOR: False
	STRETCH: True
	HOST_MASS: True


In this case, a magnitude smearing of 0.1 would be applied to the Hsiao model at all wavelengths,
and some host mass and stretch effects are applied as well based on functions you will 
define in the next sections. 

Warping Effects
===============

The following sections contain all of the various wavelength/phase dependent effects that you want
to apply to your SED. In this case, based on the **[FLAGS]** section, you must have a "HOST_MASS" section
and a "STRETCH" section. You can name effects whatever you want, as long as the name of your section and the corresponding
name in the **[FLAGS]** section are identical. Effects are either **"HOST"** effects, or they are **"SN"** effects.
They both require the same basic variables, but have different names. In all cases, you will have the following variables
defined:

1. SCALE_DIST_PEAK
	* The PEAK of an (a)symmetric Gaussian that will define the distribution for the scale parameter
2. SCALE_DIST_SIGMA
  	* The "low" and "high" standard deviations of the same distribution
3. SCALE_DIST_LIMITS
  	* The lower and upper cutoff you would like for the same distribution 

HOST Effects
============

Creating a HOST warping effect section requires the following
variable:

HOST_FUNCTION

  * A file name to be read that contains a gridded list parameters and values like the following:

::

	#phase wavelength host_mass value
	-20 	1000 	9.0		25.75805 
	-20 	1010 	9.0		25.64852
	-20 	1020 	9.0		25.53899
	-20 	1030 	9.0		25.42946
	-20 	1040 	9.0		25.31993
	-20 	1050 	9.0		25.2104
	     		...

You may pass relevant host galaxy parameters (e.g. host mass, redshift, etc.) to the warpSED function,
or you may create a distribution here to randomize the relevant parameter. You can do this by 
defining another (a)symmetric Gaussian:

1. HOST_PARAM_DIST_PEAK
	* The PEAK of an (a)symmetric Gaussian that will define the distribution for the host galaxy parameter
2. HOST_PARAM_DIST_SIGMA
	* The "low" and "high" standard deviations of the same distribution
3. HOST_PARAM_DIST_LIMITS
	* The lower and upper cutoff you would like for the same distribution 

OR by passing a filename that contains a distribution of the relevant host parameter:

HOST_PARAM_DIST_FILE: host_mass_distribution.txt

SN Effects
==========

These are exactly the same as host effects, with different labels. To create a SN effect, follow
the same directions listed for the HOST effect, but replace "HOST" with "SN" in each variable name. 


Creating Effects in the Params File
===================================

You must now define a section for each warping effect, with these variables. For our current example,
where I have defined host_mass and stretch effects in my **[FLAGS]** section, I must define these two
sections. If I do not define a **[FLAGS]** section, then whatever sections that exist apart from
the **[MAIN]** section are assumed to be warping effects. One such section might look like the
following:


::

	[STRETCH]

	SN_FUNCTION: color_func.dat
	SCALE_DIST_PEAK: 0.5
	SCALE_DIST_SIGMA: 1.0 0.7
	SCALE_DIST_LIMITS: -2.5 2.5

All together, after adding in the HOST_MASS section as well, a **BYOSED.params** file might look something like this:

::

	[MAIN]

	SED_FILE: Hsiao07.dat
	MAGSMEAR: 0.0
	MAGOFF: 0.0

	[FLAGS]

	COLOR: False
	STRETCH: True
	VELOCITY: False
	SFR: False
	METALLICITY: False
	HOST_MASS: True

	[HOST_MASS]

	HOST_FUNCTION: host_mass_func.dat

	SCALE_DIST_PEAK: 1
	SCALE_DIST_SIGMA: .000001 .000001
	SCALE_DIST_LIMITS: .99999 1.00001

	HOST_PARAM_DIST_PEAK: 10
	HOST_PARAM_DIST_SIGMA: 2 2
	HOST_PARAM_DIST_LIMITS: 5 20

	[STRETCH]

	SN_FUNCTION: stretch_func.dat

	SCALE_DIST_PEAK: 0.5
	SCALE_DIST_SIGMA: 1.0 0.7
	SCALE_DIST_LIMITS: -2.5 2.5

Or, if you do not define a flags section, host_mass and stretch will automatically be used as 
warping effects with the following **BYOSED.params** file:

::

	[MAIN]

	SED_FILE: Hsiao07.dat
	MAGSMEAR: 0.0
	MAGOFF: 0.0

	[HOST_MASS]

	HOST_FUNCTION: host_mass_func.dat

	SCALE_DIST_PEAK: 1
	SCALE_DIST_SIGMA: .000001 .000001
	SCALE_DIST_LIMITS: .99999 1.00001

	HOST_PARAM_DIST_PEAK: 10
	HOST_PARAM_DIST_SIGMA: 2 2
	HOST_PARAM_DIST_LIMITS: 5 20


	[STRETCH]

	SN_FUNCTION: stretch_func.dat

	SCALE_DIST_PEAK: 0.5
	SCALE_DIST_SIGMA: 1.0 0.7
	SCALE_DIST_LIMITS: -2.5 2.5

Final Notes
===========

Now you can replace the Hsiao template with your own template SED, and start adding in warping
effects. This warping process is designed so that as many effects as you would like can be
included. Each effect is applied multiplicatively to the baseline model. For the example file 
above, the final flux would look like this 

.. math::

   F(\lambda,\phi)=A\times H(\lambda,\phi)\Big[1+S(\lambda,\phi)s+M(\lambda,\phi,M)m\Big]

Where here F is the final flux, H is the Hsiao template, S is the defined stretch function,
M is the defined host mass function, s is the scale parameter pulled from the distribution defined
for the stretch function, m is the scale parameter pulled from the distribution defined 
for the host mass function, and M is the host mass itself, pulled from the parameter 
distribution defined for the host mass effect. 
In principle this could look like the following if you had N such effects:

.. math::

   F(\lambda,\phi)=A\times H(\lambda,\phi)\Big[1+X_1(\lambda,\phi)x_1+X_2(\lambda,\phi)x_2+...+X_N(\lambda,\phi)x_N\Big]


Combining HOST and SN Effects
=============================

You can also define an effect that involves both HOST and SN parameters. Perhaps you want an effect that combines host mass
and velocity. You might then have a params file that looks like this (the entire effect still only gets one scale parameter):

::

	[MAIN]

	SED_FILE: Hsiao07.dat
	MAGSMEAR: 0.0
	MAGOFF: 0.0

	[HOST_MASS_VELOCITY]

	HOST_FUNCTION: host_mass_func.dat

	SCALE_DIST_PEAK: 1
	SCALE_DIST_SIGMA: .000001 .000001
	SCALE_DIST_LIMITS: .99999 1.00001

	HOST_PARAM_DIST_PEAK: 10
	HOST_PARAM_DIST_SIGMA: 2 2
	HOST_PARAM_DIST_LIMITS: 5 20

	SN_FUNCTION: gridded_velocity.dat

	SN_PARAM_DIST_FILE: velocity_hist_data.txt


In this case, the final flux would be calculated in the following way:

.. math::

   F(\lambda,\phi)=A\times H(\lambda,\phi)\Big[1+V(\lambda,\phi,v)sM(\lambda,\phi,M)\Big]

Where here F is the final flux, H is the Hsiao template, V is the velocity component of the HOST_MASS_VELOCITY effect,
s is the scale factor, and M is the host mass component of the HOST_MASS_VELOCITY effect. This generalizes to N such
effects in the following way:

.. math::
	
	F(\lambda,\phi)=A\times H(\lambda,\phi)\Big[1+SN_1(\theta_{SN})s_1G_1(\theta_{SN},\theta_{HOST})+SN_2(\theta_{SN})s_2G_2(\theta_{SN},\theta_{HOST})+...+SN_N(\theta_{SN})s_NG_N(\theta_{SN},\theta_{HOST})\Big]

Example Files
=============

These are example files that can be used for your :download:`sed_file <./example_files/Hsiao07.dat>` and :download:`BYOSED.params <./example_files/BYOSED.params>`.
The host mass and stretch functions are defined by accompanying :download:`host mass <./example_files/gridded_mass.dat>` and :download:`stretch <./example_files/stretch_func.dat>` files.






















