from Clustering.clustering_functions import calculate_rmsd, plot_clustering_global, plot_clustering_local


def main_clustering(input_dir, output_dir, samples_replicas_dict, reference_sample, labels_dict, position_dict,
                    time_step, n_groups, clustering_global_threshold, clustering_local_threshold):
    """
    Plots global and, if possible (residue of interest available), local clusterings of a group of simulations

    :param input_dir: Directory in which the input folders are located
    :param output_dir: Directory in which the results will be saved
    :param samples_replicas_dict: Dictionary with all samples and their respective replicas
    :param reference_sample: Sample to which all the others are compared (usually WT)
    :param labels_dict: Dictionary with the label for each sample
    :param position_dict: Dictionary with the mutation position for each sample
    :param time_step: Time step between frames of the input .xtc trajectory (in ns)
    :param n_groups: Higher number of the groups present in the .gro file
    :param clustering_global_threshold: Threshold for global clustering (in nm, usually 0.6)
    :param clustering_local_threshold: Threshold for local clustering (in nm, usually 0.3)
    :return: Nothing, just generates and saves clustering plots in the output directory.
    """
    for replica in samples_replicas_dict[reference_sample]:
        calculate_rmsd(input_dir, output_dir, reference_sample, replica)
    ref_label = labels_dict[reference_sample]
    for sample in samples_replicas_dict:
        if sample != reference_sample:
            label = labels_dict[sample]
            for replica in samples_replicas_dict[sample]:
                calculate_rmsd(input_dir, output_dir, sample, replica)
            for replica in samples_replicas_dict[sample]:
                for ref_replica in samples_replicas_dict[reference_sample]:
                    print(ref_replica, replica)
                    plot_clustering_global(input_dir, output_dir, reference_sample, sample, ref_replica, replica,
                                           ref_label, label, clustering_global_threshold, time_step)
                    if position_dict[sample]:
                        residue = position_dict[sample]
                        print(input_dir, output_dir, reference_sample, sample, ref_replica, replica,
                                              ref_label, label, clustering_local_threshold, residue, time_step,
                                              n_groups)
                        plot_clustering_local(input_dir, output_dir, reference_sample, sample, ref_replica, replica,
                                              ref_label, label, clustering_local_threshold, residue, time_step,
                                              n_groups)
