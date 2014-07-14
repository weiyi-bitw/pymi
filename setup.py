from distutils.core import setup
from distutils.extension import Extension
import numpy

cmodule = Extension('_c_bsplinemi', include_dirs = ['pymi', numpy.get_include()], sources = ['pymi/utils.c'])

setup(
    name='pymi',
    version='0.22',
    author='Wei-Yi Cheng',
    author_email='wei-yi.cheng@mssm.edu',
    ext_modules = [cmodule],
    packages=['pymi'],
    include_package_data = True,
    scripts=['bin/getAllMIWith'],
    url='http://hidysabc.com',
    license='LICENSE.txt',
    description='Python library for MI calculation.',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy==1.8.1"
    ],
)

