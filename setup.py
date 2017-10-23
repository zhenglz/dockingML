from setuptools import setup

"""
Description of how to make a python package

https://python-packaging.readthedocs.io/en/latest/command-line-scripts.html

"""

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='dockingML',
      version='1.0',
      long_description=readme(),
      description='Docking and Machine Learning for virtual screening',
      url='https://github.com/zhenglz/dockingML',
      author='zhenglz',
      author_email='zhenglz@outlook.com',
      license='NA',
      packages=['dockingML'],
      install_requires=[
          'numpy',
          'pandas',
          'mpi4py',
      ],
      include_package_data=True,
      zip_safe=False,

      entry_points = {'console_scripts': ['dockingML=dockingML.cmd:main'],},
      )