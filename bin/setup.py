"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

from setuptools import setup, find_packages

NAME="redundans"

here = path.abspath(path.dirname(__file__))
# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

py_modules = [NAME, "fasta2homozygous", "fastq2insert_size", "fastq2mates", "fastq2sspace", "filterReads", ]

setup(name=NAME,
      version='0.13a4',
      description='',
      long_description=long_description,
      author='Leszek Pryszcz',
      author_email='l.p.pryszcz+distutils@gmail.com',
      url='https://github.com/lpryszcz/%s'%NAME,
      license='GPLv3',

      install_requires=["FastaIndex", "pyScaf", ],

      packages=find_packages(exclude=['test', ]), 
      py_modules=py_modules,

      entry_points={"console_scripts":
                    ["%s = %s:main"%(m, m) for m in py_modules]}, 
      
      classifiers=[
          # How mature is this project? Common values are  3 - Alpha 4 - Beta 5 - Production/Stable
          'Development Status :: 4 - Beta',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
      ],

      keywords='heterozygous genome assembly scaffolding reduction paired-end long-reads synteny',
      
     )
