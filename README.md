# **FEP Analysis v1.0**

**Author:** &nbsp;&nbsp;Nour Aldin Kahlous & Israel Cabeza de Vaca Lopez</br>
**Email:**  &nbsp;&nbsp;&nbsp;Nouraldinkahlous@gmail.com </br>
**Place:** &nbsp;&nbsp;&nbsp; Jens Carlsson Lab at Uppsala University</br>
**Date:** &nbsp;&nbsp;  2020-2021

**Description:** Comprehensive MD/FEP Analysis 

**Note:** The

## **INSTALATION**

MD/FEP Analysis requires the following:

1 - Python 3.7.4

    conda create --name py37 python=3.7.4

2 - activate the conda environment:
    conda activate py37

3 - Download and install FEP Analysis toolkit:
    git clone https://github.com/nouraldinkahlous/Molecular_Dynamics_and_Free_Energy_Perturbation.git

    pip -e install fepanalysis



Optional:
 Check your installation by running the tests included in the fepanalysis folder.

cd fepanalysis; python -m unittest test_fepanalysis/test_fepanalysis.py


TIP: Do not forget to activate your py37 environment before using fepanalysis.

conda activate py37