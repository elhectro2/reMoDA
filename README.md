# reMoDA
Automated calculation of unfolding detection metrics for relaxation Molecular Dynamics

## Introduction
reMoDA is a Python workflow designed for the detection of unfolding events in Molecular Dynamics simulations, especially tuned for relaxation Molecular Dynamics. It calculates four different kinds of metrics, built in four different modules: Clustering, energetics, multianalysis and PCA. It has been developed by Helena García-Cebollada at the [Javier Sancho lab](https://javiersancholab.bifi.es/)

## Installation
reMoDA is presented as a Python 3 library, which performs several analysis using other libraries and software. To successfully run reMoDA, those dependencies must be installed and their execution commands must be correctly indicated (see [Setting up the environment](#setting-up-the-environment)).
### Dependencies
Please make sure all the listed dependencies are installed and correctly running. The versions used for development and testing are indicated in parenthesis.
#### Python libraries
Required Python (3.10.12) libraries, as listed in [requirements.txt](requirements.txt), are:
* matplotlib(3.6.2)
* scikit-learn(1.2.0)
* numpy(1.23.5)
* pandas(1.5.2)
* statsmodels(0.13.5)
* scipy(1.9.3)
* tsmoothie(1.0.4)
* mdtraj(1.9.7)
#### External software
* TM-score (downloadable from the [Zhang group webpage](https://zhanggroup.org/TM-score/))
* GROMACS (2021.4-2, the latest version can be downloaded from the [GROMACS webpage](https://manual.gromacs.org/current/download.html))
* DSSP (executable provided with reMoDA at Multianalysis/Analysis/dssp)
### Setting up the environment
The environment can be set by changing the four variables in [command_dirs.py](command_dirs.py):
* analysis_dir: The path in which reMoDA is located (i.e., where the [main_analysis.py](main_analysis.py) script is located)
* dssp: Directory in which the dssp executable is located (by default, the location of the provided executable at [Multianalysis/Analysis/dssp](Multianalysis/Analysis/dssp))
* gromacs: Alias for calling GROMACS in the bash terminal
* TMscore: Alias for calling TMScore in the bash terminal

## Usage
### Preparing the input folders

### Setting the parameters

### Launching reMoDA
