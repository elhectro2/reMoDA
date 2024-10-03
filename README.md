# reMoDA
Automated calculation of unfolding detection metrics for relaxation Molecular Dynamics

## Introduction
reMoDA is a Python workflow designed for the detection of unfolding events in Molecular Dynamics simulations, especially tuned for relaxation Molecular Dynamics. It calculates four different kinds of metrics, built in four different modules: Clustering, energetics, multianalysis and PCA. It has been developed by Helena García-Cebollada at the [Javier Sancho lab](https://javiersancholab.bifi.es/) and is distributed under a [MIT license](LICENSE) through [this GitHub repository](https://github.com/elhectro2/reMoDA).

## Installation
reMoDA is presented as a Python 3 library, which performs several analysis using other libraries and software. To successfully run reMoDA, those dependencies must be installed and their execution commands must be correctly indicated (see [Setting up the environment](#setting-up-the-environment)).
### Dependencies
Please make sure all the listed dependencies are installed and correctly running. The versions used for development and testing are indicated in parenthesis.
#### Python libraries
Required Python (3.10.12) libraries, as listed in [requirements.txt](requirements.txt), are:
* matplotlib(3.6.2)
* scikit-learn(1.2.0)
* numpy(1.23.5)
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
reMoDA can be used by simply editing the [main_analysis.py](main_analysis.py) script, to tune the parameters and then launching it. However, the input folders must comply to a given structure.
### Preparing the input folders
There must be a single folder for each sample (one or more) and, inside each folder, one for each replica of the sample. Even if only one replica is available, a folder for that replica must be created. For each replica, the folder must contain the following files:
* A .xtc file, containing the trajectory
* A .edr file, containing the saved energetic properties along the trajectory
* A .tpr file, containing the topology.
* A .gro file (optional). If it is not present in the folder, it will be created from the first frame of the .xtc trajectory.
### Setting the parameters
The main_analysis.py script is mainly composed by tunable variables to set the correct parameters for the calculation. The parameters can be grouped based on their definition.
#### Directories and samples parameters
* **input_dir:** Directory in which the folders of all samples are saved
* **output_dir:** Directory in which all the resulting files will be created
* **samples_replicas_dict:** Dictionary indicating all replicas for all samples. It can be introduced manually or, in case all samples have the same replicas (same number and name of replicas), using the **samples** and **r_list** parameters.
* **samples:** (Optional) Name for the samples that have the same replicas
* **r_list:** (Optional) Name for the replicas of the previously defined **samples**
* **reference_sample:** The sample to which all mutants are to be compared (usually the WT one or the one where no unfolding is expected)
* **labels_dict:** (optional) Label for each of the samples. If the dictionary is empty, it takes the name of the sample as label.
* **position_dict:** (optional) Position of interest for each sample. If none is given, the program will automatically check if the name of the sample is in the O#M format, where O is the 1-letter code for the original amino acid, # is the position of the mutation and M is the 1-letter code for the mutated amino acid. If this is the format of the name of the sample, # will be used as position. If not, it will return None, and no local calculations will be performed.
#### Simulation parameters
If these parameters are not known, they can be obtained by running `gmx rms` on the xtc file and inspecting the resulting file using `vim` or any other text processor.
* **time_step:** Time between saved frames of the simulation (in ns)
* **duration:** Timelength of the full trajectory (in ns)
* **final:** Time of the last frame of the trajectory (in ns)
* **starting:** (automatically calculated) Starting time of the trajectory (in ns). It may or may not be 0.
* **n_groups:** Number of groups defined in the .tpr file (used when automating gmx make_ndx for local analyses)
#### Multianalysis parameters
* **protein_type:** Protein folding according to the secondary structure: alpha (a), beta (b), or alpha+beta(a+b). If the value of protein_type is not a or b, it will be assigned a+b. Used for secondary structure calculations.
#### Clustering parameters
* **clustering_global_threshold:** Minimum distance between two different clusters in global clustering (in nm). 0.3 is recommended.
* **clustering_local_threshold:** Minimum distance between two different clusters in local clustering (in nm), if performed. 0.075 is recommended.
### Launching reMoDA
After setting the correct parameters, reMoDA can be launched by executing the [main_analysis.py](main_analysis.py) script. This can be done using an IDE or directly on the terminal. We strongly recommend using the terminal, as the calculations are computationally expensive and the IDE may need more resources, hampering reMoDA calculations.
### Examples of usage
To illustrate the usage of reMoDA, the most common scenarios are shown as examples of usage. The usage of reMoDA in a single sample is not recommended for the Clustering and PCA modules, as their units are arbitrary and the results cannot be interpreted with certainty, so those modules only work by comparing samples. Leaving that scenario out, the most common ones are samples with the same number of replicas, equally named between samples; samples with one replica each; and the most complex case, with more than two replicas and different number of replicas for each one.
#### Two samples with the same number of replicas
We'll use as an example the hypothetical case of a protein with a mutation of interest, so the WT and the A237V mutant are simulated, with three replicas each. The input directory should be like this:
```
input_dir  
├── WT
|   ├── r1
|   |   ├── traj.xtc
|   |   ├── energies.edr
|   |   └── top.tpr
|   ├── r2
|   |   ├── traj.xtc
|   |   ├── energies.edr
|   |   └── top.tpr
|   └── r3
|       ├── traj.xtc
|       ├── energies.edr
|       └── top.tpr
└── A273V
    ├── r1
    |   ├── traj.xtc
    |   ├── energies.edr
    |   └── top.tpr
    ├── r2
    |   ├── traj.xtc
    |   ├── energies.edr
    |   └── top.tpr
    └── r3
        ├── traj.xtc
        ├── energies.edr
        └── top.tpr
```
As the three replicas for each sample are named equally, the simplified version for defining the replicas dictionary can be used:
```python
samples = ["WT", "A273V"]
r_list = ["r1", "r2", "r3"]
samples_replicas_dict = {x: r_list for x in samples}
```
The reference sample would be the WT in this case, even though in two samples analysis it makes no clear difference, it is a required variable. As the labels will be the same as the sample names and the mutation follows the defined name, there's no need to define any of theses dictionaries. However, for illustrative purposes, we will define them in this example:
```python
reference_sample = "WT"
labels_dict = {"WT": "WT", "A273V": "A273V"}
position_dict = {"WT": None, "A273V": 273}
```
For this example, let's assume the simulation parameters are the original ones set in main_analysis.py, but we are using a mainly alpha protein:
```python
time_step = 1  # Time between saved frames of the simulation (in ns)
duration = 1000  # Timelength of the full trajectory (in ns)
final = 1000  # Time of the last frame of the trajectory (in ns)
starting = final - duration  # Starting time of the trajectory (in ns). It may or may not be 0.
n_groups = 18  # Number of groups defined in the .tpr file (used when automating gmx make_ndx for local analyses)

# Parameters for multianalysis.
protein_type = "a"  # Protein folding according to the secondary structure: alpha (a), beta (b), or alpha+beta(a+b).
# If the value of protein_type is not a or b, it will be assigned a+b. Used for secondary structure calculations.
protein_type = autocomplete_protein_type(protein_type)

# Parameters for clustering
clustering_global_threshold = 0.3  # In nm
clustering_local_threshold = 0.075  # In nm
```
Then, the only remaining step is to launch the main script (we recommend launching via terminal):
```bash
python3 main_analysis.py
```

#### Two samples with one replica each
We'll use as an example the hypothetical case of a protein at two temperatures, so the 298K and the 347K temperatures are simulated, with one replica each. The input directory should be like this:
```
input_dir  
├── 298K
|   └── r1
|       ├── traj.xtc
|       ├── energies.edr
|       └── top.tpr
└── 347K
    └── r1
        ├── traj.xtc
        ├── energies.edr
        └── top.tpr
```
Please, make sure that a folder for the unique replica is created inside the sample folder. If the input files are located directly in the sample folder, they will not be found by reMoDA.

As the unique replica folder for each sample are named equally, the simplified version for defining the replicas dictionary can be used:
```python
samples = ["298K", "347K"]
r_list = ["r1"]
samples_replicas_dict = {x: r_list for x in samples}
```
In this case, sample names and labels are identical, so there is no need defininf them. However, this example uses different temperatures, so there is not a residue of interest. In this case, we can define the residue of interest as ```None``` and no local calculations will be performed.
```python
reference_sample = "298K"
labels_dict = {}
autocomplete_labels(samples_replicas_dict, labels_dict)
position_dict = {"298K": None, "347K": None}
```
Considering the simulation parameters are correct, the only step left is to launch the main_analysis.py script as in the previous example:
```bash
python3 main_analysis.py
```

#### Three samples with different number of replicas
We'll use as an example the hypothetical case of a protein with two mutations of interest, located at different residues, so the WT, the D43E and the L123I mutants are simulated, with different number and names of replicas for each sample. The input directory should be like this:
```
input_dir  
├── WT
|   ├── r1
|   |   ├── traj.xtc
|   |   ├── energies.edr
|   |   └── top.tpr
|   ├── r2
|   |   ├── traj.xtc
|   |   ├── energies.edr
|   |   └── top.tpr
|   └── r3
|       ├── traj.xtc
|       ├── energies.edr
|       └── top.tpr
├── D43E
|   ├── r4
|   |   ├── traj.xtc
|   |   ├── energies.edr
|   |   └── top.tpr
|   └── r5
|       ├── traj.xtc
|       ├── energies.edr
|       └── top.tpr
└── L123I
    └── r6
        ├── traj.xtc
        ├── energies.edr
        └── top.tpr
```
In this case, the simplified version for defining the dictionary cannot be used, so the full dictionary is indicated manually:
```python
samples_replicas_dict = {"WT": ["r1", "r2", "r3"],
                         "D43E": ["r4", "r5",],
                         "L123I": ["r6"]}
```
As there are more than two samples, in this case selecting a different reference sample will make a difference. If WT is selected, the comparisons WT_vs_D43E and WT_vs_L123I will be generated, while if D43E is selected, the generated comparisons will be D43E_vs_WT and D43E_vs_L123I. In this case, we select the WT as reference sample and leave the label and position dictionaries empty so that they are autocompleted:
```python
reference_sample = WT
labels_dict = {}
autocomplete_labels(samples_replicas_dict, labels_dict)
position_dict = {}
autocomplete_positions(samples_replicas_dict, position_dict)
```
Assuming the parameters of the simulation are correct, the only step left is launching the script:
```bash
python3 main_analysis.py
```
