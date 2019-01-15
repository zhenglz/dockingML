from setuptools import setup

"""
Description of how to make a python package

https://python-packaging.readthedocs.io/en/latest/command-line-scripts.html

"""


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='dockingML',
      version='1.2.0',
      long_description=readme(),
      description='Molecular dynamics analysis, Docking and Machine Learning for Computer Aid Drug discovery',
      url='https://github.com/zhenglz/dockingML',
      author='zhenglz',
      author_email='zhenglz@outlook.com',
      license='GPL-3.0',
      packages=['dockml', 'mdanaly', 'automd'],
      install_requires=[
          'numpy',
          'pandas',
          'scikit-learn',
          'matplotlib',
          'networkx',
          'mdtraj',
      ],
      include_package_data=True,
      zip_safe=False,
      python_requires='>=3.5',
      )
