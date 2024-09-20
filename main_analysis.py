from helper_functions import *
from Clustering.main_clustering import main_clustering
from Energetics.main_energetics import main_energetics
from Multianalysis.main_multianalysis import main_multianalysis
from pca.main_pca import main_pca

# Directories on which the input folder can be found and on which the results will be saved. The input and output
# directories are checked in order to homogenize their format, to make sure they end in /. The folder structure for
# intermediate calculations and for the results will be generated after the samples and replicas are defined.
input_dir = "/media/elhectro2/TOSHIBA_EXT/Trabajo_Unizar/reMoDA/Test_sets/CHEK2/Simulations/Trajectories-to-analyze_SMenao"
output_dir = "/media/elhectro2/TOSHIBA_EXT/Trabajo_Unizar/reMoDA/Test_sets/CHEK2/Simulations/Trajectories-to-analyze_SMenao/output"
input_dir = check_dir(input_dir)
output_dir = check_dir(output_dir)

# Dictionary with all the samples and the replicas present in each one. Each sample may or may not have different
# number of replicas with different names. The replicas do not necessarily need to follow the "r1" nomenclature.
# samples_replicas_dict = {"WT": ["r1", "r2", "r3", "r4", "r5"]}
samples_replicas_dict = {# "A392V-358K": ["r1_A", "r1_B", "r2_A", "r2_B"],
                         # "A392V-378K": ["r1_A", "r1_B", "r3_A", "r3_B"],
                         # "A392V-398K": ["r1_A", "r1_B", "r2_A", "r2_B", "r3_A", "r3_B"],
                         # "H345Y-378K": ["r1_A", "r1_B", "r2_A", "r2_B"],
                         # "H345Y-398K": ["r2_A", "r2_B"],
                         # "P471L-358K": ["r1_A", "r1_B", "r2_A", "r2_B", "r3_A", "r3_B"],
                         # "P471L-378K": ["r1_A", "r1_B", "r2_A", "r2_B", "r3_A", "r3_B"],
                         # "P471L-398K": ["r1_A", "r1_B", "r2_A", "r2_B", "r3_A", "r3_B"],
                         "R474C-378K": ["r1_A", "r1_B", "r2_A", "r2_B", "r3_A", "r3_B"],
                         "R474C-398K": ["r1_A", "r1_B", "r2_A", "r2_B", "r3_A", "r3_B"],
                         "WT-358K": ["r1_A", "r1_B"]}


# If the replicas are named the same on all samples, this piece of code can be easier to modify
"""samples = ["298K", "310K"]
r_list = ["r1", "r2", "r3", "r4", "r5"]
samples_replicas_dict = {x: r_list for x in samples}  # Generates the combinations of all samples with all replicas
"""
# The sample to which all mutants are to be compared (usually the WT one) is defined
reference_sample = "WT-358K"

# Names for each of the samples. If you want to use the same name as in samples_replicas_dict, don't introduce it 
# in the dictionary, it will autocomplete automatically
labels_dict = {# "A392V-358K": "A392V-358K",
               # "A392V-378K": "A392V-378K",
               # "A392V-398K": "A392V-398K",
               # "H345Y-378K": "H345Y-378K",
               # "H345Y-398K": "H345Y-398K",
               # "P471L-358K": "P471L-358K",
               # "P471L-378K": "P471L-378K",
               # "P471L-398K": "P471L-398K",
               "R474C-378K": "R474C-378K",
               "R474C-398K": "R474C-398K",
               "WT-358K": "WT-358K"}

autocomplete_labels(samples_replicas_dict, labels_dict)

# Position of the mutated residue for each of the samples. If you don't want to perform local analyses or the sample is 
# the WT simulation, introduce the value None. If no value is introduced for a sample, the program will try to obtain 
# the position by trimming out the first and last character of the sample string. If the result is an integer, it will 
# be used as the position. If not, the value None will be assigned and no local analysis would be performed on that 
# sample.
position_dict = {"WT-358K": None, "A392V-358K": 392, "A392V-378K": 392, "A392V-398K": 392, "H345Y-378K": 345,
                 "H345Y-398K": 345, "P471L-358K": 471, "P471L-378K": 471, "P471L-398K": 471, "R474C-378K": 474,
                 "R474C-398K": 474}
# position_dict = {"WT": None, "R158C": 158}
# position_dict = {"298K": None, "310K": None}
autocomplete_positions(samples_replicas_dict, position_dict)

# The folder structure for intermediate calculations and for saving the results is generated
create_folder_structure(input_dir, output_dir, samples_replicas_dict, reference_sample, position_dict)

# Parameters of the simulation
time_step = 1  # Time between saved frames of the simulation (in ns)
duration = 1000  # Timelength of the full trajectory (in ns)
final = 1000  # Time of the last frame of the trajectory (in ns)
starting = final - duration  # Starting time of the trajectory (in ns). It may or may not be 0.
n_groups = 18  # Number of groups defined in the .tpr file (used when automatizing gmx make_ndx for local analyses)

# Parameters for multianalysis.
protein_type = "a+b"  # Protein folding according to the secondary structure: alpha (a), beta (b), or alpha+beta(a+b).
# If the value of protein_type is not a or b, it will be assigned a+b. Used for secondary structure calculations.
protein_type = autocomplete_protein_type(protein_type)

# Parameters for clustering
clustering_global_threshold = 0.3  # In nm
clustering_local_threshold = 0.075  # In nm. Only used if local analyses are performed.

# Main functions for analysis: Multiple metrics calculation, 2D-RMSD clustering, PCA of metrics, energetic disaggregate
main_multianalysis(input_dir, output_dir, samples_replicas_dict, reference_sample, labels_dict, position_dict,
                   time_step, final, starting, protein_type)
main_pca(input_dir, output_dir, samples_replicas_dict, reference_sample, labels_dict, position_dict, time_step, final,
         starting, protein_type, output_dir)
main_clustering(input_dir, output_dir, samples_replicas_dict, reference_sample, labels_dict, position_dict, time_step,
                n_groups, clustering_global_threshold, clustering_local_threshold)
main_energetics(input_dir, output_dir, samples_replicas_dict, reference_sample, labels_dict)

# Cleans the original directories for avoiding duplication of the resulting files.
# final_cleaning(input_dir, samples_replicas_dict)
