from Multianalysis.multianalysis_functions import global_multianalysis, local_multianalysis
from Multianalysis.Plotting.plot_global import plot_global_compared


def main_multianalysis(input_dir, output_dir, samples_replicas_dict, reference_sample, labels_dict, position_dict,
                       time_step, final, starting, protein_type):
    """
    Calculates and plots several glabal and (if possible) local metrics related to rMD.

    :param input_dir: Directory in which the input folders are located
    :param output_dir: Directory in which the results will be saved
    :param samples_replicas_dict: Dictionary with all samples and their respective replicas
    :param reference_sample: Sample to which all the others are compared (usually WT)
    :param labels_dict: Dictionary with the label for each sample
    :param position_dict: Dictionary with the mutation position for each sample
    :param time_step: Time step between frames of the input .xtc trajectory (in ns)
    :param final: Ending time of the simulation (in ns)
    :param starting: Starting time of the simulation (in ns, usually 0)
    :param protein_type: Folding according to secondary structure. Expected inputs are a/alpha, b/beta or a+b/alpha+beta
    :return: Saves all single and compared plots with the reference sample, as well as the raw results
    """
    all_data_dict = {}
    ref_label = labels_dict[reference_sample]
    all_data_dict[reference_sample] = global_multianalysis(input_dir, output_dir, reference_sample,
                                                           samples_replicas_dict[reference_sample],
                                                           time_step, final, starting, protein_type, ref_label)
    for sample in samples_replicas_dict:
        if sample != reference_sample:
            sample_label = labels_dict[sample]
            all_data_dict[sample] = global_multianalysis(input_dir, output_dir, sample, samples_replicas_dict[sample],
                                                         time_step, final, starting, protein_type, sample_label)

            plot_global_compared(all_data_dict, labels_dict, [reference_sample, sample], output_dir, protein_type)
            if position_dict[sample]:
                reference_replicas = samples_replicas_dict[reference_sample]
                replicas = samples_replicas_dict[sample]
                residue = position_dict[sample]
                local_multianalysis(input_dir, output_dir, reference_sample, sample, reference_replicas, replicas,
                                    labels_dict, residue)
