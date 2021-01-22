from setuptools import setup, find_packages


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='fepanalysis',
      version='1.0',
      description='Alchemical transformations MD analisys',
      long_description=long_description,
      long_description_content_type="text/markdown",
      classifiers=[
        'Development Status :: Release',
        'License :: MIT License',
        'Programming Language :: Python :: 3.7',
        "Operating System :: OS Linux",
        'Topic :: Molecular Mechanics :: Force Field parameterization',
      ],
      keywords='FEP ',
      url='https://github.com/Isra3l/ligpargen',
      author='Israel Cabeza de Vaca Lopez',
      author_email='israel.cabezadevaca@icm.uu.se',
      license='MIT License',
      packages=find_packages(),
      install_requires=[
          'markdown', 'numpy'
      ],
      entry_points={
          'console_scripts': ['feptool=fepanalysis.feptool:main'],
      },
      include_package_data=True,
      python_requires='==3.7.*',
      zip_safe=False)