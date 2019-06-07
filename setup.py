from distutils.core import setup

from setuptools import find_packages

from staramr import __version__

classifiers = """
Development Status :: 4 - Beta
Environment :: Console
License :: OSI Approved :: Apache Software License
Intended Audience :: Science/Research
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bio-Informatics
Programming Language :: Python :: 3.5
Programming Language :: Python :: 3.6
Operating System :: POSIX :: Linux
""".strip().split('\n')

setup(name='staramr',
      version=__version__,
      description='Scans genome contigs against ResFinder, PlasmidFinder, and PointFinder databases',
      author='Aaron Petkau',
      author_email='aaron.petkau@gmail.com',
      url='https://github.com/phac-nml/staramr',
      license='Apache v2.0',
      classifiers=classifiers,
      install_requires=[
          'biopython>=1.70',
          'pandas>=0.23.0',
          'GitPython>=2.1.3',
          'xlsxwriter>=1.0.2',
          'numpy>=1.12.1',
          'green>=2.13.0',
          'coloredlogs>=10.0'
      ],
      packages=find_packages(),
      include_package_data=True,
      scripts=['bin/staramr'],
      setup_requires = ['green']
      )
