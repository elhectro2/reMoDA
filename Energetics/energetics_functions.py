
from Energetics.Plotting.plot_global import plot_energetics_simple
from Energetics.Analysis.energetics import *
gmx = command_dir["gromacs"]


def energetics(input_dir, output_dir, sample, replicas, label):
    """
    Calculates all global metrics of a sample, plots them and returns the data for all metrics to compute compared plots

    :param input_dir: Directory in which the input folders are located
    :param output_dir: Directory in which the results will be saved
    :param sample: Sample on which the calculation is to be performed
    :param replicas: List of replicas available for the sample
    :param time_step: Time step between frames of the input .xtc trajectory (in ns)
    :param final: Ending time of the simulation (in ns)
    :param starting: Starting time of the simulations (in ns, usually 0)
    :param protein_type: Folding according to secondary structure. Expected inputs are a/alpha, b/beta or a+b/alpha+beta
    :param label: Label for the sample
    :return: All data for all global metrics to compute compared plots
    """
    specific_output_dir = output_dir + sample + "/Energetics/Raw_results/"
    specific_output_plots = output_dir + sample + "/Energetics/Plots/"
    sample_data_dict = {}
    for replica in replicas:
        replica_data_dict = {}
        specific_in_dir = input_dir + sample + "/" + replica + "/"

        # Calculating all global functions
        calculate_lj_pn(specific_in_dir, specific_output_dir, replica)
        calculate_lj_pp(specific_in_dir, specific_output_dir, replica)
        calculate_improp(specific_in_dir, specific_output_dir, replica)
        calculate_prop(specific_in_dir, specific_output_dir, replica)
        calculate_coul_pn(specific_in_dir, specific_output_dir, replica)
        calculate_coul_pp(specific_in_dir, specific_output_dir, replica)

        # Reading all global functions
        replica_data_dict["coul_pn"] = read_energetics(specific_output_dir, "coul_pn", replica)
        replica_data_dict["coul_pp"] = read_energetics(specific_output_dir, "coul_pp", replica)
        replica_data_dict["improper"] = read_energetics(specific_output_dir, "improper", replica)
        replica_data_dict["proper"] = read_energetics(specific_output_dir, "proper", replica)
        replica_data_dict["LJ_pp"] = read_energetics(specific_output_dir, "LJ_pp", replica)
        replica_data_dict["LJ_pn"] = read_energetics(specific_output_dir, "LJ_pn", replica)

        sample_data_dict[replica] = replica_data_dict

    # Plotting function
    plot_energetics_simple(sample_data_dict, label, specific_output_plots, "Energetics")

    # Returns all data for comparison
    return sample_data_dict
