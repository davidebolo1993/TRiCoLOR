from setuptools import setup,find_packages



setup(name='TRiCoLOR',      
  version='1.1.0',
  description='Tandem Repeats Caller for LOng Reads',
  url='https://github.com/davidebolo1993/TRiCoLOR',
  requires=['python (>= 3.6)'],
  author='Davide Bolognini, Tobias Rausch',
  author_email='davidebolognini7@gmail.com, rausch@embl.de',
  license='LICENSE.txt',
  install_requires=['pyfaidx>=0.5.5.2', 'pysam>=0.15.0', 'editdistance>=0.5.2', 'pandas>=0.22.0', 'numpy>=1.15.3', 'plotly==3.10.0', 'cyvcf2>=0.10.1', 'pybedtools>=0.8.0', 'mappy>=2.17'],
  zip_safe=False,
  packages=find_packages(),
  include_package_data=True,
  entry_points={'console_scripts': ['TRiCoLOR=TRiCoLOR.TRiCoLOR:main']}          
)
