from setuptools import setup, find_packages
 
classifiers = [
  'Development Status :: 2 - Pre-Alpha',
  'Intended Audience :: Science/Research',
  'License :: OSI Approved :: MIT License',
  'Programming Language :: Python :: 3'
]
 
setup(
  name='NllocPy',
  version='0.0.1',
  description='Nlloc python wrapper',
  long_description=open('README.md').read(),
  url='',  
  author='Luigi Sante Zampa',
  author_email='zampaluigis@gmail.com',
  license='MIT', 
  classifiers=classifiers,
  keywords='calculator', 
  packages=find_packages(),
  install_requires=["nllgrid" ],
)