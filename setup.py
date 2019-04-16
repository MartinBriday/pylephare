#! /usr/bin/env python
#

DESCRIPTION = "sedkcorr: Apply k-correction from a given SED."
LONG_DESCRIPTION = """  sedkcorr: Apply k-correction from a given SED. """

DISTNAME = 'sedkcorr'
AUTHOR = 'Martin Briday'
MAINTAINER = 'Martin Briday' 
MAINTAINER_EMAIL = 'briday@ipnl.in2p3.fr'
URL = 'https://github.com/MartinBriday/sedkcorr'
LICENSE = 'Apache 2.0'
DOWNLOAD_URL = 'https://github.com/MartinBriday/sedkcorr'
VERSION = '0.1.0'

try:
    from setuptools import setup, find_packages
    _has_setuptools = True
except ImportError:
    from distutils.core import setup

def check_dependencies():
    install_requires = []

    try:
        import propobject
    except ImportError:
        install_requires.append('propobject')

    try:
        import pyifu
    except ImportError:
        install_requires.append('pyifu')
    
    try:
        import sncosmo
    except ImportError:
        install_requires.append('sncosmo')
    
    try:
        import astropy
    except ImportError:
        install_requires.append('astropy')
       
     
    return install_requires

if __name__ == "__main__":

    install_requires = check_dependencies()

    if _has_setuptools:
        packages = find_packages()
        print(packages)
    else:
        # This should be updated if new submodules are added
        packages = ['sedkcorr']

    setup(name=DISTNAME,
          author=AUTHOR,
          author_email=MAINTAINER_EMAIL,
          maintainer=MAINTAINER,
          maintainer_email=MAINTAINER_EMAIL,
          description=DESCRIPTION,
          long_description=LONG_DESCRIPTION,
          license=LICENSE,
          url=URL,
          version=VERSION,
          download_url=DOWNLOAD_URL,
          install_requires=install_requires,
          packages=packages,
          package_data={"sedkcorr":["k_correction/filter_bandpass/*.dat"]},
          classifiers=[
              'Intended Audience :: Science/Research',
              'Programming Language :: Python :: 2.7',
              'Programming Language :: Python :: 3.5',              
              'License :: OSI Approved :: BSD License',
              'Topic :: Scientific/Engineering :: Astronomy',
              'Operating System :: POSIX',
              'Operating System :: Unix',
              'Operating System :: MacOS'],
      )










