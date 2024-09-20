from Energetics.energetics_functions import energetics
from Energetics.Plotting.plot_global import plot_energetics_compared


def main_energetics(input_dir, output_dir, samples_replicas_dict, reference_sample, labels_dict):
    """
    Calculates and plots several glabal and (if possible) local metrics related to rMD.

    :param input_dir: Directory in which the input folders are located
    :param output_dir: Directory in which the results will be saved
    :param samples_replicas_dict: Dictionary with all samples and their respective replicas
    :param reference_sample: Sample to which all the others are compared (usually WT)
    :param labels_dict: Dictionary with the label for each sample
    :return: Saves all single and compared plots with the reference sample, as well as the raw results
    """
    all_data_dict = {}
    ref_label = labels_dict[reference_sample]
    all_data_dict[reference_sample] = energetics(input_dir, output_dir, reference_sample,
                                                 samples_replicas_dict[reference_sample], ref_label)
    for sample in samples_replicas_dict:
        if sample != reference_sample:
            sample_label = labels_dict[sample]
            all_data_dict[sample] = energetics(input_dir, output_dir, sample, samples_replicas_dict[sample],
                                               sample_label)
            plot_energetics_compared(all_data_dict, labels_dict, [reference_sample, sample], output_dir, "Energetics")
