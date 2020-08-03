# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 17:35:08 2020

@author: sven
"""


from setuptools import setup
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='zonarPy',
      version='0.0.1',
      description='Read and process Zooglider Zonar data',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Sven Gastauer',
      url='https://github.com/SvenGastauer/zonarPy',
      download_url = 'https://github.com/user/zonarPy/archive/0.0.1.tar.gz',
      author_email='sgastauer@ucsd.edu',
      license='MIT',
      packages=['zonarPy'],
      keywords = ['Zooglider', 'Zonar', 'Python', 'acoustics','oceanography'],
      install_requires=[
          'xarray',
          'numpy',
          'pandas',
          'datetime',
          'matplotlib',
          'time',
          'os',
          'astral',
          'simplekml',
          'requests',
          'datetime',
          'glob',
          'mayavi',
          'collections',
          'plotly'
        ],
      classifiers=[
        'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        'Intended Audience :: Developers',      # Define that your audience are developers
        'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: MIT License',   # Again, pick a license
        'Programming Language :: Python :: 3',      #Specify which pyhton versions that you want to support
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
      ],
      include_package_data=True,
      package_data={'': ['data/*.MT']},
)