#!/usr/bin/python3
# -*- coding: utf-8 -*-

# Note: To use the 'upload' functionality of this file, you must:
#   $ pipenv install twine --dev

import io
import os
import sys
import glob
from shutil import rmtree

from setuptools import find_packages, setup, Command

# Package meta-data.
NAME = 'loopflapper'
DESCRIPTION = 'Alter the angle of a CDRH3 loop model. Predict optimum angle'
URL = ''
EMAIL = 'zcbtlm0@ucl.ac.uk, lilian_denzler@protonmail.com'
AUTHOR = 'Lilian Denzler'
REQUIRES_PYTHON = '>=3.6.0'
VERSION = '0.0.0'

# What packages are required for this module to be executed?
REQUIRED = [
	'biopandas',
	'keras',
	'seaborn',
	'tensorflow',
	'yellowbrick',
	'xgboost',
	'statsmodels',
	'scikit-learn'

]
TEST_REQUIRED=[]

# What packages are optional?
EXTRAS = {
	# 'fancy feature': ['django'],
}

# The rest you shouldn't have to touch too much :)
# ------------------------------------------------
# Except, perhaps the License and Trove Classifiers!
# If you do change the License, remember to change the Trove Classifier for that!

here = os.path.abspath(os.path.dirname(__file__))
os.system("echo 'export QUALILOOP={}' >> ~/.bashrc".format(os.path.join(here)))
# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
try:
	with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
		long_description = '\n' + f.read()
except FileNotFoundError:
	long_description = DESCRIPTION

# Load the package's __version__.py module as a dictionary.
about = {}
if not VERSION:
	project_slug = NAME.lower().replace("-", "_").replace(" ", "_")
	with open(os.path.join(here, project_slug, '__version__.py')) as f:
		exec(f.read(), about)
else:
	about['__version__'] = VERSION


class UploadCommand(Command):
	"""Support setup.py upload."""

	description = 'Build and publish the package.'
	user_options = []

	@staticmethod
	def status(s):
		"""Prints things in bold."""
		print('\033[1m{0}\033[0m'.format(s))

	def initialize_options(self):
		pass

	def finalize_options(self):
		pass

	def run(self):
		try:
			self.status('Removing previous builds…')
			rmtree(os.path.join(here, 'dist'))
		except OSError:
			pass

		self.status('Building Source and Wheel (universal) distribution…')
		os.system('{0} setup.py sdist bdist_wheel --universal'.format(sys.executable))

		self.status('Uploading the package to PyPI via Twine…')
		os.system('twine upload dist/*')

		self.status('Pushing git tags…')
		os.system('git tag v{0}'.format(about['__version__']))
		os.system('git push --tags')

		sys.exit()

class InstallExternal(Command):
	"""Support setup.py external."""

	description = '''Install and build external libraries and tools needed for the qualiloop package. The following packages are installed:
	ProFit, ECalc, Bioptools, Bioplib. All code belongs to Prof. Andrew C.R. Martin and is accessible at:
	https://github.com/orgs/ACRMGroup/repositories '''
	user_options = []

	@staticmethod
	def status(s):
		"""Prints things in bold."""
		print('\033[1m{0}\033[0m'.format(s))

	def initialize_options(self):
		pass

	def finalize_options(self):
		pass

	def run(self):
		try:
			self.status('Removing previous builds…')
			tar_files=glob.glob('*.tar.gz*')
			ProFit_files=glob.glob('ProFit*')
			profit_files=glob.glob('profit*')
			bioptools_files=glob.glob('bioptools*')
			for rm_list in [tar_files,ProFit_files,profit_files, bioptools_files]:
				for file in rm_list:
				  os.remove(file)

		except OSError:
			pass
		self.status('Pulling Bioplib')
		os.system('wget -P {} https://github.com/ACRMGroup/bioplib/archive/refs/heads/master.zip'.format(os.path.join(here,"bin")))
		os.system('unzip {} -d {}'.format(os.path.join(here,"bin","master.zip"),os.path.join(here,"bin")))
		os.system('rm -r {}'.format(os.path.join(here,"bin","master.zip")))
		self.status('Installing Bioplib')
		makefile=open(os.path.join(here,"bin","bioplib-master","src","Makefile"),"r")
		list_of_lines=makefile.readlines()
		list_of_lines[1]=f"DEST	=	{os.path.join(here,'bin','bioplib-master')}	"
		self.status('HERE')
		makefile=open(os.path.join(here,"bin","bioplib-master","src","Makefile"),"w")
		makefile.writelines(list_of_lines)
		makefile.close()

		makefile=open(os.path.join(here,"bin","bioplib-master","pdbtagvars","Makefile"),"r")
		list_of_lines=makefile.readlines()
		list_of_lines[0]=f"CFLAGS = -g -ansi -pedantic -Wall -I{os.path.join(here,'bin','bioplib-master')}/include	"
		list_of_lines[2]=f"LFLAGS = -L{os.path.join(here,'bin','bioplib-master')}/lib -lbiop -lgen -lxml2	"	
		self.status('HERE')
		makefile=open(os.path.join(here,"bin","bioplib-master","pdbtagvars","Makefile"),"w")
		makefile.writelines(list_of_lines)
		makefile.close()
		

		os.system(f'(cd {os.path.join(here,"bin","bioplib-master/src")} && make && make doxygen && make install && make installdata )')
		

		self.status('Pulling Bioptools')
		os.system('sudo apt-get install libxml2 libxml2-dev')
		os.system('wget -P {} https://github.com/ACRMGroup/bioptools/archive/V1.9.tar.gz'.format(os.path.join(here,"bin")))
		os.system('mkdir {}'.format(os.path.join(here,"bin", "bioptools")))
		os.system('tar -xf {} -C {}'.format(os.path.join(here,"bin","V1.9.tar.gz"),os.path.join(here,"bin", "bioptools")))
		os.system('rm -r {}'.format(os.path.join(here,"bin","V1.9.tar.gz")))
		self.status('Installing Bioptools')
		os.system('(cd {} && ./makemake.pl -prefix={} -bioplib && make && make install)'.format(os.path.join(here,"bin","bioptools/bioptools-1.9/src"),os.path.join(here,"bin", "bioptools")))
		os.system("echo 'export PATH=$PATH:{}' >> ~/.bashrc".format(os.path.join(here,"bin")))
		os.system("echo 'export DATADIR={}' >> ~/.bashrc".format(os.path.join(here,"bin","bioptools/bioptools-1.9","data")))
	

		self.status('Pulling ProFit')
		os.system('wget -P {} http://www.bioinf.org.uk/software/profit/235216/profit.tar.gz'.format(os.path.join(here,"bin")))
		os.system('tar -xf {} -C {}'.format(os.path.join(here,"bin","profit.tar.gz"),os.path.join(here,"bin")))
		os.system('rm -r {}'.format(os.path.join(here,"bin","profit.tar.gz")))
		self.status('Installing ProFit')
		os.system('(cd {} && make)'.format(os.path.join(here,"bin","ProFit_V3.3","src")))
		#os.system("echo 'export DATADIR={}' >> ~/.bashrc".format(os.path.join(here,"bin","ProFit_V3.3","data")))

		self.status('Pulling ECalc')
		os.system('wget -P {} https://github.com/AndrewCRMartin/ecalc/archive/refs/heads/master.zip'.format(os.path.join(here,"bin")))
		os.system('unzip {} -d {}'.format(os.path.join(here,"bin","master.zip"),os.path.join(here,"bin")))
		os.system('rm -r {}'.format(os.path.join(here,"bin","master.zip")))
		self.status('Installing Ecalc')
		makefile=open(os.path.join(here,"bin","ecalc-master","src","Makefile"),"r")
		list_of_lines=makefile.readlines()
		list_of_lines[1]=f"INCFLAGS = -I{os.path.abspath(os.path.join(here,'bin','bioptools/bioptools-1.9/src/libsrc/bioplib'))}"
		list_of_lines[2]=f"LIBFLAGS = -L{os.path.abspath(os.path.join(here,'bin','bioptools/bioptools-1.9/src/libsrc/bioplib'))}"
		makefile=open(os.path.join(here,"bin","ecalc-master","src","Makefile"),"w")
		makefile.writelines(list_of_lines)
		makefile.close()
		os.system('(cd {} && make LIBFLAGS={} INCFLAGS={} && make clean)'.format(os.path.join(here,"bin","ecalc-master","src"),os.path.abspath(os.path.join(here,'bin','bioptools/bioptools-1.9/src/libsrc/bioplib')), os.path.abspath(os.path.join(here,'bin','bioptools/bioptools-1.9/src/libsrc/bioplib'))))
		os.system("echo 'export ECALCDATA={}' >> ~/.bashrc".format(os.path.join(here,"bin","ecalc","data")))

		sys.exit()

# Where the magic happens:
setup(
	name=NAME,
	version=about['__version__'],
	description=DESCRIPTION,
	long_description=long_description,
	long_description_content_type='text/markdown',
	author=AUTHOR,
	author_email=EMAIL,
	python_requires=REQUIRES_PYTHON,
	url=URL,
	packages=find_packages(include=['qualiloop', 'qualiloop.*'],exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
	# If your package is a single module, use this instead of 'packages':
	# py_modules=['mypackage'],
	entry_points={
		 'console_scripts': ['makemodel=qualiloop.model_experimenter:main',
							'makeprediction=qualiloop.final_predictor:param_parser',
							'loopflapper=qualiloop.flapper_full:param_parser'],

		 'gui_scripts': []
	 },
	scripts=['bin/RMS_calclocal_CA','bin/RMS_calclocal_AA','bin/RMS_calcglobal_CA','bin/RMS_calcglobal_AA'],
	install_requires=REQUIRED,
	tests_require=TEST_REQUIRED,
	extras_require=EXTRAS,
	include_package_data=True,
	license='MIT',
	classifiers=[
		# Trove classifiers
		# Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
		'License :: OSI Approved :: MIT License',
		'Programming Language :: Python',
		'Programming Language :: Python :: 3',
		'Programming Language :: Python :: 3.6',
		'Programming Language :: Python :: Implementation :: CPython',
		'Programming Language :: Python :: Implementation :: PyPy'
	],
	# $ setup.py publish support.
	cmdclass={
		'upload': UploadCommand, 'external': InstallExternal
	},
)
