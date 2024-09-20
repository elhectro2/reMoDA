from Energetics.energetics_functions import energetics
from Energetics.Plotting.filtering import *
from Energetics.Plotting.plot_global import plot_energetics_compared
import datetime


def calculate_filtered_dict(all_data_dict, filter_list):
    filtered_dict = {}
    for filter in filter_list:
        print(filter, datetime.datetime.now())
        filtered_dict[filter] = {}
        for sample in all_data_dict:
            filtered_dict[filter][sample] = {}
            for replica in all_data_dict[sample]:
                filtered_dict[filter][sample][replica] = {}
                for metric in all_data_dict[sample][replica]:
                    time_series = all_data_dict[sample][replica][metric][0]
                    values = all_data_dict[sample][replica][metric][1]
                    if filter == "lfilter":
                        values = get_y_lfilter(values, n=300, a=1)
                    elif filter == "savgol":
                        values = get_y_sav_gol(values, win_len=2500, polyorder=2, mode="nearest")
                    elif filter == "conv":
                        values = get_y_convolution(values, win_len=2500, win_type="ones")
                    elif filter == "lowess":
                        values, time_series = get_y_lowess(values, time_series, fraction=0.1)
                    elif filter == "rolling_mean":
                        values = get_y_rolling_mean(values, time_series, win_len=2500)
                    filtered_dict[filter][sample][replica][metric] = [time_series, values]
    return filtered_dict


def main_energetics(input_dir, output_dir, samples_replicas_dict, reference_sample, labels_dict, filter_list):
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
            filtered_dict = calculate_filtered_dict(all_data_dict, filter_list)
            for filter in filtered_dict:
                print(filter, datetime.datetime.now())
                plot_energetics_compared(filtered_dict[filter], labels_dict, [reference_sample, sample], output_dir,
                                         "Energetics_%s" % filter)


input_dir = "/media/elhectro2/TOSHIBA_EXT/Trabajo_Unizar/MD_Analysis/optimized_analysis_scripts/Viviana/input_skip_10/"
output_dir = "/media/elhectro2/TOSHIBA_EXT/Trabajo_Unizar/MD_Analysis/optimized_analysis_scripts/Viviana/output_skip_10/"
samples_replicas_dict = {"WT": ["r1", "r2", "r3", "r4", "r5"], "R158C": ["r1", "r2", "r3", "r4", "r5"]}
reference_sample = "WT"
labels_dict = {"WT": "WT", "R158C": "R158C"}
filter_list = ["lfilter", "savgol", "conv", "lowess", "rolling_mean"]
main_energetics(input_dir, output_dir, samples_replicas_dict, reference_sample, labels_dict, filter_list)
