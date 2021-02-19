import os
from setuptools.command.test import test as TestCommand
from distutils.core import setup
import numpy.distutils.misc_util
import sys


if sys.version_info < (3,0):
	sys.exit('Sorry, Python 2 is not supported')

class BYOSEDTest(TestCommand):

	def run_tests(self):
		import byosed
		errno = byosed.test()
		sys.exit(errno)

AUTHOR = 'Justin Pierel'
AUTHOR_EMAIL = 'jr23@email.sc.edu'
VERSION = '0.0.3'
LICENSE = 'BSD'
URL = 'byosed.readthedocs.org'

def recursive_glob(basedir, pattern):
    matches = []
    for root, dirnames, filenames in os.walk(basedir):
        for filename in fnmatch.filter(filenames, pattern):
            matches.append(os.path.join(root, filename))
    return matches

PACKAGENAME='sntd'
# Add the project-global data
pkgdatadir = os.path.join(PACKAGENAME, 'kaepora')
initfiles = os.path.join(PACKAGENAME, 'initfiles')

data_files = []
data_files.extend(recursive_glob(pkgdatadir, '*'))
data_files.extend(recursive_glob(initfiles, '*'))


data_files = [f[len(PACKAGENAME)+1:] for f in data_files]

setup(
	name='BYOSED',
	version=VERSION,
	cmdclass={'test': BYOSEDTest},
	packages=['byosed'],
        package_data={'sntd':data_files},
	include_package_data=True,
	author=AUTHOR,
	author_email=AUTHOR_EMAIL,
	license=LICENSE,
	long_description=open('README.md').read(),
        include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
	install_requires=['numpy>=1.5.0',
					  'scipy>=0.9.0',
					  'astropy',
					  'pytest-astropy',
					  'pandas',
					  'matplotlib',
					  'sncosmo'],
	)
