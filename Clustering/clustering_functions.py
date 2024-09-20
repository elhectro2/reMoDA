import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess
from command_dirs import command_dir
from gromacs.fileformats.xvg import XVG
from sklearn.cluster import AgglomerativeClustering
gmx = command_dir["gromacs"]


def calculate_rmsd(input_dir, output_dir, sample, replica):
    """
    Calculates the Root Mean Square Deviation (RMSD) of a trajectory (a replica of a sample)

    :param input_dir: Directory in which the input folders are located
    :param output_dir: Directory in which the results will be saved
    :param sample: Sample on which the calculation is to be performed
    :param replica: Replica of the sample on which the calculation is to be performed
    :return: Generates the xvg file with the RMSD calculated for each frame
    """
    specific_input_dir = input_dir + sample + "/" + replica + "/"
    specific_output_dir = output_dir + sample + "/Multianalysis/Raw_results/"
    if not os.path.isfile(specific_output_dir + "rmsd_%s.xvg" % replica):
        order = "echo 'Backbone Backbone' | %s rms -f %s/*.xtc -s %s/*.gro -o %s/rmsd_%s.xvg -tu ns" \
                % (gmx, specific_input_dir, specific_input_dir, specific_input_dir, replica)
        subprocess.call(order, shell=True)


def read_rmsd(rmsd_dir):
    """
    Using an already generated RMSD file in xvg format, reads its content

    :param rmsd_dir: Path to the RMSD xvg file.
    :return: A list of the time for each frame in the simulation and another list for the calculated RMSD for that frame
    """
    data = []
    time_series = []
    values = []
    file = open(rmsd_dir)
    for line in file:  # Reads all the file and does things
        if line[0] not in "#@":  # Non-comment lines
            line = line.split()  # Splits the line in a list of "words" (elements separated by spaces in a string)
            data.append(line)  # Data for each time is stored in a list of lists
    for item in data:
        time_series.append(float(item[0]))
        values.append(float(item[1]))
    return time_series, values


def read_xpm(file_dir):
    """
    Converts the xpm 2D-RMSD encoded matrix into a numerical matrix, which is saved as a variable

    :param file_dir: Path in which the vpm 2D-RMSD matrix file is located
    :return: The 2D-RMSD matrix translated into numbers
    """
    in_file = open(file_dir)
    line = in_file.readline()
    while line[0] != "\"":
        line = in_file.readline()
    line = in_file.readline()
    code_dict = {}
    while line[0] == "\"":
        line = line[1:].split()
        code_dict[line[0]] = float(line[5][1:-1])
        line = in_file.readline()
    result = []
    for line in in_file.readlines():
        if line[0] == "\"":
            line_results = []
            number = -2
            if line[-2] == ",":
                number -= 1
            for c in line[1:number]:
                line_results.append(code_dict[c])
            result.append(line_results)
    return result


def sele_backbone(residue, clustering_ref_in_dir, reference_sample, ref_replica):
    """
    Generates an index file with all the residues within the sphere of 1 nm around the residue of interest

    :param residue: Residue of interest (usually the mutated one)
    :param clustering_ref_in_dir: Directory in which the clustering files for the reference sample are located
    :param reference_sample: Sample to which all the others are compared (usually WT)
    :param ref_replica: Replica of the reference sample that is being compared
    :return: Generates an index ndx file with the residues for local clustering
    """
    subprocess.call("%s select -f %s%s_%s_backbone.pdb -s %s%s_%s_backbone.pdb -on %sindex_back_resid%s.ndx "
                    "-resnr number -select  \"(group \"Backbone\" and not resid %i) and same residue as within 1 "
                    "of com of resid %i\""
                    % (gmx, clustering_ref_in_dir, reference_sample, ref_replica, clustering_ref_in_dir,
                       reference_sample, ref_replica, clustering_ref_in_dir, residue, residue, residue), shell=True)


def get_clustering(distance_matrix, threshold):
    """
    Calculates the clustering from a precalculated distance matrix (usually from read_xpm)

    :param distance_matrix: Precalculated distance matrix
    :param threshold: Distance threshold to separate the clusters (in nm)
    :return: List with the cluster each frame is at.
    """
    clustering = AgglomerativeClustering(n_clusters=None, metric="precomputed", linkage="single",
                                         distance_threshold=threshold).fit_predict(distance_matrix)
    return clustering


def clustering_global(input_dir, output_dir, reference_sample, sample, ref_replica, replica,
                      clustering_global_threshold, time_step):
    """
    Calculates global clustering and returns the data

    :param input_dir: Directory in which the input folders are located
    :param output_dir: Directory in which the results will be saved
    :param reference_sample: Sample to which all the others are compared (usually WT)
    :param sample: Sample on which the calculation is to be performed
    :param ref_replica: Replica of the reference sample on which the calculation is to be performed
    :param replica: Replica of the sample on which the calculation is to be performed
    :param clustering_global_threshold: Threshold for global clustering (in nm, usually 0.6)
    :param time_step: Time step between frames of the input .xtc trajectory (in ns)
    :return: Returns the resulting data (time and cluster) for both samples
    """

    # Defining directories
    specific_ref_in_dir = input_dir + reference_sample + "/" + ref_replica + "/"
    clustering_ref_in_dir = specific_ref_in_dir + "Clustering/"
    specific_in_dir = input_dir + sample + "/" + replica + "/"
    clustering_in_dir = specific_in_dir + "Clustering/"
    specific_out_dir_clustering = "%s%s_vs_%s/Clustering/Global/%s-%s_vs_%s-%s/" % (output_dir, reference_sample, sample,
                                                                                   reference_sample, ref_replica,
                                                                                   sample, replica)

    # Calculating backbone trajectories and references
    if not os.path.exists("%s%s_%s_backbone.pdb" % (clustering_ref_in_dir, reference_sample, ref_replica)):
        subprocess.call("echo 4 | %s trjconv -f %s*.xtc -s %s*.gro -o %s%s_%s_backbone.pdb -dump 0" %
                        (gmx, specific_ref_in_dir, specific_ref_in_dir, clustering_ref_in_dir, reference_sample,
                         ref_replica), shell=True)
    if not os.path.exists("%s%s_%s_backbone.xtc" % (clustering_ref_in_dir, reference_sample, ref_replica)):
        subprocess.call("echo 4 | %s trjconv -f %s*.xtc -s %s*.gro -o %s%s_%s_backbone.xtc" %
                        (gmx, specific_ref_in_dir, specific_ref_in_dir, clustering_ref_in_dir, reference_sample,
                         ref_replica), shell=True)
    if not os.path.exists("%s%s_%s_backbone.pdb" % (clustering_in_dir, sample, replica)):
        subprocess.call("echo 4 | %s trjconv -f %s*.xtc -s %s*.gro -o %s%s_%s_backbone.pdb -dump 0" %
                        (gmx, specific_in_dir, specific_in_dir, clustering_in_dir, sample, replica), shell=True)
    if not os.path.exists("%s%s_%s_backbone.xtc" % (clustering_in_dir, sample, replica)):
        subprocess.call("echo 4 | %s trjconv -f %s*.xtc -s %s*.gro -o %s%s_%s_backbone.xtc" %
                        (gmx, specific_in_dir, specific_in_dir, clustering_in_dir, sample, replica), shell=True)

    # Calculating 2D-RMSD matrices
    if not os.path.exists("%s%s_%s_global_matrix.xpm" % (clustering_ref_in_dir, reference_sample, ref_replica)):
        subprocess.call("echo \'4 4\' | %s rms -f %s%s_%s_backbone.xtc -s %s%s_%s_backbone.pdb -o "
                        "%s%s_%s_global_matrix.xvg -m %s%s_%s_global_matrix.xpm" % (gmx, clustering_ref_in_dir,
                                                                                    reference_sample, ref_replica,
                                                                                    clustering_ref_in_dir,
                                                                                    reference_sample, ref_replica,
                                                                                    clustering_ref_in_dir,
                                                                                    reference_sample, ref_replica,
                                                                                    clustering_ref_in_dir,
                                                                                    reference_sample, ref_replica),
                        shell=True)
    if not os.path.exists("%s%s_%s_global_matrix.xpm" % (clustering_in_dir, sample, replica)):
        subprocess.call("echo \'4 4\' | %s rms -f %s%s_%s_backbone.xtc -s %s%s_%s_backbone.pdb -o "
                        "%s%s_%s_global_matrix.xvg -m %s%s_%s_global_matrix.xpm" % (gmx,
                                                                                    clustering_in_dir, sample, replica,
                                                                                    clustering_in_dir, sample, replica,
                                                                                    clustering_in_dir, sample, replica,
                                                                                    clustering_in_dir, sample, replica),
                        shell=True)
    if not os.path.exists("%s%s-%s_vs_%s-%s.xpm" % (specific_out_dir_clustering, reference_sample, ref_replica, sample,
                                                    replica)):
        subprocess.call("echo \'4 4\' | %s rms -f %s%s_%s_backbone.xtc -f2 %s%s_%s_backbone.xtc "
                        "-s %s%s_%s_backbone.pdb -o %s%s-%s_vs_%s-%s.xvg -m %s%s-%s_vs_%s-%s.xpm" %
                        (gmx, clustering_ref_in_dir, reference_sample, ref_replica, clustering_in_dir, sample, replica,
                         clustering_ref_in_dir, reference_sample, ref_replica, specific_out_dir_clustering,
                         reference_sample, ref_replica, sample, replica, specific_out_dir_clustering, reference_sample,
                         ref_replica, sample, replica),
                        shell=True)

    # Translate matrices into numbers
    matrix_1 = read_xpm("%s%s_%s_global_matrix.xpm" % (clustering_ref_in_dir, reference_sample, ref_replica))
    matrix_2 = read_xpm("%s%s_%s_global_matrix.xpm" % (clustering_in_dir, sample, replica))
    matrix_1_2 = read_xpm("%s%s-%s_vs_%s-%s.xpm" % (specific_out_dir_clustering, reference_sample, ref_replica, sample,
                                                    replica))
    matrix_2_1 = np.array(matrix_1_2[::-1]).T.tolist()[::-1]

    # Joins matrices and gets the appropriate format
    matrix_a = []
    for i in range(len(matrix_2)):
        result = matrix_1_2[i] + matrix_2[i]
        matrix_a.append(result)
    n_2 = len(matrix_a)
    matrix_b = []
    for i in range(len(matrix_1)):
        result = matrix_1[i] + matrix_2_1[i]
        matrix_b.append(result)
    n_1 = len(matrix_b)
    matrix = matrix_a + matrix_b
    matrix = matrix[::-1]

    # Performs clustering and assigns cluster numbers in order of appearance in the metatrajectory.
    result = get_clustering(matrix, clustering_global_threshold)
    numbering = 1
    numbering_dict = {}

    # Gets cluster number (y) at each time of simulation (x) and saves it as a .xvg file for both trajectories
    y_1 = []
    x_1 = []
    time_1 = []
    for i in range(n_1):
        time_1.append(i * time_step * 1000)
        x_1.append(i * time_step)
        if result[i] not in numbering_dict:
            numbering_dict[result[i]] = numbering
            numbering += 1
        y_1.append(numbering_dict[result[i]])
    a = XVG(array=[time_1, y_1])
    a.write("%sclusters_%s_%s.xvg" % (specific_out_dir_clustering, reference_sample, ref_replica))
    y_2 = []
    x_2 = []
    time_2 = []
    for i in range(n_2):
        x_2.append(i * time_step)
        time_2.append(i * time_step * 1000)
        if result[i + n_1] not in numbering_dict:
            numbering_dict[result[i + n_1]] = numbering
            numbering += 1
        y_2.append(numbering_dict[result[i + n_1]])
    a = XVG(array=[time_2, y_2])
    a.write("%sclusters_%s_%s.xvg" % (specific_out_dir_clustering, sample, replica))
    return x_1, y_1, x_2, y_2


def clustering_local(input_dir, output_dir, reference_sample, sample, ref_replica, replica, clustering_local_threshold,
                     residue, time_step, n_groups):
    """
    Calculates local clustering around the residue of interest and returns the data

    :param input_dir: Directory in which the input folders are located
    :param output_dir: Directory in which the results will be saved
    :param reference_sample: Sample to which all the others are compared (usually WT)
    :param sample: Sample on which the calculation is to be performed
    :param ref_replica: Replica of the reference sample on which the calculation is to be performed
    :param replica: Replica of the sample on which the calculation is to be performed
    :param clustering_local_threshold: Threshold for local clustering (in nm, usually 0.6)
    :param residue: Residue of interest (usually the mutated one)
    :param time_step: Time step between frames of the input .xtc trajectory (in ns)
    :param n_groups: Higher number of the groups present in the .gro file
    :return: Returns the resulting data (time and cluster) for both samples
    """
    # Defining directories
    specific_ref_in_dir = input_dir + reference_sample + "/" + ref_replica + "/"
    clustering_ref_in_dir = specific_ref_in_dir + "Clustering/"
    specific_in_dir = input_dir + sample + "/" + replica + "/"
    clustering_in_dir = specific_in_dir + "Clustering/"
    specific_out_dir_clustering = "%s%s_vs_%s/Clustering/Local/%s-%s_vs_%s-%s/" % (output_dir, reference_sample, sample,
                                                                                   reference_sample, ref_replica,
                                                                                   sample, replica)

    # Calculating backbone trajectories and references
    if not os.path.exists("%s%s_%s_backbone.pdb" % (clustering_ref_in_dir, reference_sample, ref_replica)):
        subprocess.call("echo 4 | %s trjconv -f %s*.xtc -s %s*.gro -o %s%s_%s_backbone.pdb -dump 0" %
                        (gmx, specific_ref_in_dir, specific_ref_in_dir, clustering_ref_in_dir, reference_sample,
                         ref_replica), shell=True)
    if not os.path.exists("%s%s_%s_backbone.xtc" % (clustering_ref_in_dir, reference_sample, ref_replica)):
        subprocess.call("echo 4 | %s trjconv -f %s*.xtc -s %s*.gro -o %s%s_%s_backbone.xtc" %
                        (gmx, specific_ref_in_dir, specific_ref_in_dir, clustering_ref_in_dir, reference_sample,
                         ref_replica), shell=True)
    if not os.path.exists("%s%s_%s_backbone.pdb" % (clustering_in_dir, sample, replica)):
        subprocess.call("echo 4 | %s trjconv -f %s*.xtc -s %s*.gro -o %s%s_%s_backbone.pdb -dump 0" %
                        (gmx, specific_in_dir, specific_in_dir, clustering_in_dir, sample, replica), shell=True)
    if not os.path.exists("%s%s_%s_backbone.xtc" % (clustering_in_dir, sample, replica)):
        subprocess.call("echo 4 | %s trjconv -f %s*.xtc -s %s*.gro -o %s%s_%s_backbone.xtc" %
                        (gmx, specific_in_dir, specific_in_dir, clustering_in_dir, sample, replica), shell=True)

    # Generating and selecting index for local calculations, 1 nm around the residue of interest in the first frame
    if not os.path.exists("%sindex_back_resid%s.ndx" % (clustering_ref_in_dir, residue)):
        sele_backbone(residue, clustering_ref_in_dir, reference_sample, ref_replica)

    # Calculate local matrices
    if not os.path.exists("%s%s_%s_resid%s_matrix.xpm" % (clustering_ref_in_dir, reference_sample, ref_replica,
                                                          str(residue))):
        subprocess.call("echo \'4 %s\' | %s rms -f %s%s_%s_backbone.xtc -s %s%s_%s_backbone.pdb -o "
                        "%s%s_%s_resid%s_matrix.xvg -m %s%s_%s_resid%s_matrix.xpm -n %sindex_back_resid%s.ndx"
                        % (str(n_groups + 1), gmx, clustering_ref_in_dir, reference_sample, ref_replica,
                           clustering_ref_in_dir, reference_sample, ref_replica, clustering_ref_in_dir,
                           reference_sample, ref_replica, str(residue), clustering_ref_in_dir, reference_sample,
                           ref_replica, str(residue), clustering_ref_in_dir, str(residue)), shell=True)
    if not os.path.exists("%s%s_%s_resid%s_matrix.xpm" % (clustering_in_dir, sample, replica, str(residue))):
        subprocess.call("echo \'4 %s\' | %s rms -f %s%s_%s_backbone.xtc -s %s%s_%s_backbone.pdb -o "
                        "%s%s_%s_resid%s_matrix.xvg -m %s%s_%s_resid%s_matrix.xpm -n %sindex_back_resid%s.ndx"
                        % (str(n_groups + 1), gmx, clustering_in_dir, sample, replica, clustering_in_dir, sample, replica,
                           clustering_in_dir, sample, replica, str(residue), clustering_in_dir, sample, replica,
                           str(residue), clustering_ref_in_dir, str(residue)), shell=True)
    if not os.path.exists("%s%s-%s_vs_%s-%s_resid%s.xpm" % (specific_out_dir_clustering, reference_sample, ref_replica,
                                                            sample, replica, str(residue))):
        subprocess.call("echo \'4 %s\' | %s rms -f %s%s_%s_backbone.xtc -f2 %s%s_%s_backbone.xtc "
                        "-s %s%s_%s_backbone.pdb -o %s%s-%s_vs_%s-%s_resid%s.xvg -m %s%s-%s_vs_%s-%s_resid%s.xpm "
                        "-n %sindex_back_resid%s.ndx" %
                        (str(n_groups + 1), gmx, clustering_ref_in_dir, reference_sample, ref_replica, clustering_in_dir,
                         sample, replica, clustering_ref_in_dir, reference_sample, ref_replica,
                         specific_out_dir_clustering, reference_sample, ref_replica, sample, replica, str(residue),
                         specific_out_dir_clustering, reference_sample, ref_replica, sample, replica, str(residue),
                         clustering_ref_in_dir, str(residue)), shell=True)

    # Translate matrices into numbers
    matrix_1 = read_xpm("%s%s_%s_resid%s_matrix.xpm" % (clustering_ref_in_dir, reference_sample, ref_replica,
                                                        str(residue)))
    matrix_2 = read_xpm("%s%s_%s_resid%s_matrix.xpm" % (clustering_in_dir, sample, replica, str(residue)))
    matrix_1_2 = read_xpm(
        "%s%s-%s_vs_%s-%s_resid%s.xpm" % (specific_out_dir_clustering, reference_sample, ref_replica, sample,
                                          replica, str(residue)))
    matrix_2_1 = np.array(matrix_1_2[::-1]).T.tolist()[::-1]

    # Joins matrices and gets the appropriate format
    matrix_a = []
    for i in range(len(matrix_2)):
        result = matrix_1_2[i] + matrix_2[i]
        matrix_a.append(result)
    n_2 = len(matrix_a)
    matrix_b = []
    for i in range(len(matrix_1)):
        result = matrix_1[i] + matrix_2_1[i]
        matrix_b.append(result)
    n_1 = len(matrix_b)
    matrix = matrix_a + matrix_b
    matrix = matrix[::-1]

    # Performs clustering and assigns cluster numbers in order of appearance in the metatrajectory.
    result = get_clustering(matrix, clustering_local_threshold)
    numbering = 1
    numbering_dict = {}

    # Gets cluster number (y) at each time of simulation (x) and saves it as a .xvg file for both trajectories
    y_1 = []
    x_1 = []
    time_1 = []
    for i in range(n_1):
        time_1.append(i * time_step * 1000)
        x_1.append(i * time_step)
        if result[i] not in numbering_dict:
            numbering_dict[result[i]] = numbering
            numbering += 1
        y_1.append(numbering_dict[result[i]])
    a = XVG(array=[time_1, y_1])
    a.write("%sclusters_%s_%s_resid%s.xvg" % (specific_out_dir_clustering, reference_sample, ref_replica,
                                              str(residue)))
    y_2 = []
    x_2 = []
    time_2 = []
    for i in range(n_2):
        x_2.append(i * time_step)
        time_2.append(i * time_step * 1000)
        if result[i + n_1] not in numbering_dict:
            numbering_dict[result[i + n_1]] = numbering
            numbering += 1
        y_2.append(numbering_dict[result[i + n_1]])
    a = XVG(array=[time_2, y_2])
    a.write("%sclusters_%s_%s_resid%s.xvg" % (specific_out_dir_clustering, sample, replica, str(residue)))
    return x_1, y_1, x_2, y_2


def plot_clustering_global(input_dir, output_dir, reference_sample, sample, ref_replica, replica, ref_label, label,
                           clustering_global_threshold, time_step):
    """
    Plots the global clustering compared results as Time vs cluster for both samples

    :param input_dir: Directory in which the input folders are located
    :param output_dir: Directory in which the results will be saved
    :param reference_sample: Sample to which all the others are compared (usually WT)
    :param sample: Sample on which the calculation is to be performed
    :param ref_replica: Replica of the reference sample on which the calculation is to be performed
    :param replica: Replica of the sample on which the calculation is to be performed
    :param ref_label: Label for the reference replica of the reference sample
    :param label: Label for the replica of the sample
    :param clustering_global_threshold: Threshold for global clustering (in nm, usually 0.6)
    :param time_step: Time step between frames of the input .xtc trajectory (in ns)
    :return: Saves the plot of the resulting global data (time vs cluster) for both samples
    """
    specific_out_dir_clustering = "%s%s_vs_%s/Clustering/Global/%s-%s_vs_%s-%s" % (output_dir, reference_sample, sample,
                                                                                   reference_sample, ref_replica,
                                                                                   sample, replica)
    # Calculate and obtain RMSD data for comparison plot
    # Debería estar ya calculate_rmsd(input_dir, output_dir, reference_sample, ref_replica)
    # Debería estar ya calculate_rmsd(input_dir, output_dir, sample, replica)
    ref_rmsd_dir = output_dir + reference_sample + "/Multianalysis/Raw_results/rmsd_%s.xvg" % ref_replica
    rmsd_dir = output_dir + sample + "/Multianalysis/Raw_results/rmsd_%s.xvg" % replica
    rmsd_1 = read_rmsd(ref_rmsd_dir)
    rmsd_2 = read_rmsd(rmsd_dir)

    # Plot the RMSD data and set correct format for the plots
    matplotlib.rc('xtick', labelsize=12)
    matplotlib.rc('ytick', labelsize=12)
    fig = plt.figure(figsize=(5, 7.5))
    axes = []
    ax = fig.add_subplot(2, 1, 1)
    ax.plot(rmsd_1[0], rmsd_1[1], label=str(ref_label), c="blue", alpha=1)
    ax.plot(rmsd_2[0], rmsd_2[1], label=str(label), c="orange", alpha=0.7)
    ax.set_xlabel("Time (ns)", fontsize=16)
    ax.set_ylabel("RMSD (nm)", fontsize=16)
    ax.legend(loc="best", fontsize=14)
    ax.margins(0.1)
    axes.append(ax)
    ax = fig.add_subplot(2, 1, 2)

    # Calculate global clustering and plot it.
    x_1, y_1, x_2, y_2 = clustering_global(input_dir, output_dir, reference_sample, sample, ref_replica, replica,
                                           clustering_global_threshold, time_step)
    ymax = (max([max(y_1), max(y_2)]))
    step_size = int((ymax - 1) / 4)
    start = int(((ymax - 1) % 4) / 2) + 1
    if step_size < 1:
        step_size = 1
    ylabels = range(start, ymax + 1, step_size)
    new_yticks = ylabels
    plt.yticks(new_yticks, ylabels)
    ax.plot(x_1, y_1, label=str(ref_label), c="blue", alpha=1)
    ax.plot(x_2, y_2, label=str(label), c="orange", alpha=0.7)
    ax.set_xlabel("Time (ns)", fontsize=16)
    ax.set_ylabel("Cluster number", fontsize=16)
    ax.legend(loc="best", fontsize=14)
    ax.margins(0.1)
    axes.append(ax)
    fig.tight_layout()
    fig.savefig(specific_out_dir_clustering +
                "%s-%s_vs_%s-%s global_clustering_%snm_min.png" % (reference_sample, ref_replica, sample, replica,
                                                               str(clustering_global_threshold)))
    plt.close(fig)


def plot_clustering_local(input_dir, output_dir, reference_sample, sample, ref_replica, replica, ref_label, label,
                          clustering_local_threshold, residue, time_step, n_groups):
    """
    Plots the local clustering compared results as Time vs cluster for both samples

    :param input_dir: Directory in which the input folders are located
    :param output_dir: Directory in which the results will be saved
    :param reference_sample: Sample to which all the others are compared (usually WT)
    :param sample: Sample on which the calculation is to be performed
    :param ref_replica: Replica of the reference sample on which the calculation is to be performed
    :param replica: Replica of the sample on which the calculation is to be performed:param ref_label: Label for the reference replica of the reference sample
    :param label: Label for the replica of the sample
    :param clustering_local_threshold: Threshold for local clustering (in nm, usually 0.6)
    :param residue: Residue of interest (usually the mutated one)
    :param time_step: Time step between frames of the input .xtc trajectory (in ns)
    :param n_groups: Higher number of the groups present in the .gro file
    :return: Saves the plot of the resulting local data (time vs cluster) for both samples
    """
    specific_out_dir_clustering = "%s%s_vs_%s/Clustering/Local/%s-%s_vs_%s-%s" % (output_dir, reference_sample, sample,
                                                                                  reference_sample, ref_replica,
                                                                                  sample, replica)

    # Obtains and plots local clustering data
    x_1, y_1, x_2, y_2 = clustering_local(input_dir, output_dir, reference_sample, sample, ref_replica, replica,
                                          clustering_local_threshold, residue, time_step, n_groups)
    plt.plot(x_1, y_1, label=ref_label)
    plt.plot(x_2, y_2, label=label)
    ymax = (max([max(y_1), max(y_2)]))
    ylabels = range(1, ymax + 1)
    new_yticks = ylabels
    plt.yticks(new_yticks, ylabels)
    plt.xlabel("Time(ns)")
    plt.ylabel("Cluster number")
    plt.legend(loc="best")
    plt.savefig(specific_out_dir_clustering + "%s-%s_vs_%s-%s_resid%s_clustering_%snm_min.png"
                % (reference_sample, ref_replica, sample, replica, str(residue), str(clustering_local_threshold)))
    plt.close()
