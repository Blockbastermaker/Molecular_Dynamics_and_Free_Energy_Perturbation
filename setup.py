from setuptools import setup, find_packages


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='fepanalysis',
      version='1.0',
      description='Alchemical transformations MD/FEP analysis',
      long_description=long_description,
      long_description_content_type="text/markdown",
      classifiers=[
        'Development Status :: Release',
        'License :: MIT License',
        'Programming Language :: Python :: 3.7.4',
        "Operating System :: OS Linux",
        'Topic :: Alchemical transformations analysis',
      ],
      keywords='FEP ',
      url='https://github.com/nouraldinkahlous/Molecular_Dynamics_and_Free_Energy_Perturbation',
      author='Nour Aldin Kahlous & Israel Cabeza de Vaca Lopez',
      author_email='Nouraldinkahlous@gmail.com',
      license='MIT License',
      packages=find_packages(),
      install_requires=[
          'markdown','numpy',"seaborn"
      ],
      entry_points={
          'console_scripts': ['feptool=fepanalysis.feptool:main'],
      },
      include_package_data=True,
      python_requires='==3.7.4',
      zip_safe=False)