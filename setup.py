from distutils.core import setup

setup(
    name='pymibspline',
    version='0.22',
    author='Wei-Yi Cheng',
    author_email='wei-yi.cheng@mssm.edu',
    packages=['pymibspline'],
    include_package_data = True,
    scripts=['bin/getAllMIWith'],
    url='http://hidysabc.com',
    license='LICENSE.txt',
    description='Python library for MI calculation using B-spline method.',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy==1.8.1"
    ],
)

