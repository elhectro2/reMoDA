import os
import shutil


def check_dir(directory):
    """
    Makes sure the directory ends with a slash(/), to simplify further steps.

    :param directory: The string of the directory to check
    :return: The string of the directory ending in one single slash(/)
    """
    if directory[-1] == "/":
        return directory
    else:
        return directory + "/"


def check_makedirs(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def create_folder_structure(input_dir, output_dir, samples_replicas_dict, reference_sample, position_dict):
    """
    Creates the correct folder structure for both the input and output directories.
    
    :param input_dir: Directory in which the input folders are located
    :param output_dir: Directory in which the results will be saved
    :param samples_replicas_dict: Dictionary with all samples and their respective replicas
    :param reference_sample: Sample to which all the others are compared (usually WT)
    :param position_dict: Dictionary with the mutation position for each sample
    :return: Nothing, only generates folders
    """
    check_makedirs(output_dir + "All_PCA_global")
    for sample in samples_replicas_dict:
        specific_out_dir = output_dir + sample + "/"
        check_makedirs(specific_out_dir + "Multianalysis/Plots")
        check_makedirs(specific_out_dir + "Multianalysis/Raw_results")
        check_makedirs(specific_out_dir + "Energetics/Plots")
        check_makedirs(specific_out_dir + "Energetics/Raw_results")
        if sample != reference_sample:
            specific_out_dir_compared = output_dir + reference_sample + "_vs_" + sample + "/"
            check_makedirs(specific_out_dir_compared + "Multianalysis/Plots")
            check_makedirs(specific_out_dir_compared + "Energetics/Plots")
        for replica in samples_replicas_dict[sample]:
            specific_in_dir = input_dir + sample + "/" + replica + "/"
            check_makedirs(specific_in_dir + "Multianalysis")
            check_makedirs(specific_in_dir + "Clustering")
            check_makedirs(output_dir + "Original_input/" + sample + "/" + replica)
            # check_makedirs(output_dir + sample + "/PCA_global/" + replica)
            check_makedirs(output_dir + "All_PCA_global/" + sample + "_" + replica)
            if sample != reference_sample:
                specific_out_dir_clustering = output_dir + reference_sample + "_vs_" + sample + "/Clustering/Global/"
                specific_out_dir_pca = output_dir + reference_sample + "_vs_" + sample + "/PCA/Global/"
                check_makedirs("%s%s_%s" % (specific_out_dir_pca, sample, replica))
                check_makedirs("%sAll" % specific_out_dir_pca)
                for reference_replica in samples_replicas_dict[reference_sample]:
                    check_makedirs("%s%s-%s_vs_%s-%s" % (specific_out_dir_clustering, reference_sample,
                                                         reference_replica, sample, replica))
                    check_makedirs("%s%s_%s" % (specific_out_dir_pca, reference_sample, reference_replica))
                if position_dict[sample]:
                    specific_out_dir_clustering = output_dir + reference_sample + "_vs_" + sample + \
                                                  "/Clustering/Local/"
                    specific_out_dir_pca = output_dir + reference_sample + "_vs_" + sample + "/PCA/Local/"
                    check_makedirs("%s%s_%s" % (specific_out_dir_pca, sample, replica))
                    check_makedirs("%sAll" % specific_out_dir_pca)
                    for reference_replica in samples_replicas_dict[reference_sample]:
                        check_makedirs("%s%s-%s_vs_%s-%s" % (specific_out_dir_clustering, reference_sample,
                                                             reference_replica, sample, replica))
                        check_makedirs("%s%s_%s" % (specific_out_dir_pca, reference_sample, reference_replica))


def autocomplete_labels(samples_replicas_dict, labels_dict):
    """
    Checks if the dictionary with the labels contains all samples and autocompletes it with the names of the samples.

    :param samples_replicas_dict: Dictionary with all samples and their respective replicas
    :param labels_dict: Dictionary with the label for each sample
    :return: Autocompleted dictionary with labels for all the samples
    """
    for sample in samples_replicas_dict:
        if sample not in labels_dict:
            labels_dict[sample] = sample
    return labels_dict


def autocomplete_positions(samples_replicas_dict, position_dict):
    """
    Checks if the dictionary with the mutation positions contains all samples and autocompletes it. If the sample is
    named using one-letter standard mutation nomenclature, it will take the position from there. If not, it will
    return None for that sample.

    :param samples_replicas_dict: Dictionary with all samples and their respective replicas
    :param position_dict: Dictionary with the mutation position for each sample
    :return: Autocompleted dictionary with mutation positions (or None) for all the samples
    """
    for sample in samples_replicas_dict:
        if sample not in position_dict:
            try:
                position_dict[sample] = int(sample[1:-1])
            except ValueError:
                position_dict[sample] = None
    return position_dict


def autocomplete_protein_type(protein_type):
    """
    Converts the input of the user for the protein folding type into the format the program requires.

    :param protein_type: Folding according to secondary structure. Expected inputs are a/alpha, b/beta or a+b/alpha+beta
    :return: returns the correct protein_type in a format that the main program recognizes. If the input is unexpected,
    it will return alpha+beta.
    """
    if protein_type.lower() == "a" or protein_type.lower() == "alpha":
        return "alpha"
    elif protein_type.lower() == "b" or protein_type.lower() == "beta":
        return "beta"
    else:
        return "alpha+beta"


def final_cleaning(input_dir, samples_replicas_dict):
    """
    Removes duplicated files of the intermediate calculations of the input directory, leaving just the initial files

    :param input_dir: Directory in which the input folders are located
    :param samples_replicas_dict: Dictionary with all samples and their respective replicas
    :return: Nothing, only removes folders
    """
    for sample in samples_replicas_dict:
        for replica in samples_replicas_dict[sample]:
            specific_in_dir = input_dir + sample + "/" + replica + "/"
            shutil.rmtree(specific_in_dir + "Multianalysis")
            shutil.rmtree(specific_in_dir + "Clustering")
