from distutils.core import setup

from setuptools import find_packages

from staramr import __version__

setup(name='staramr',
      version=__version__,
      description='Scans genome contigs against ResFinder and PointFinder databases',
      author='Aaron Petkau',
      author_email='aaron.petkau@gmail.com',
      url='https://github.com/phac-nml/staramr',
      install_requires=[
          'biopython>=1.70',
          'pandas>=0.22.0',
          'six>=1.11.0',
          'tox>=2.9.1',
          'GitPython>=2.1.3',
          'xlsxwriter>=1.0.2',
      ],
      packages=find_packages(exclude=['tests']),
      include_package_data=True,
      scripts=['bin/staramr']
      )
