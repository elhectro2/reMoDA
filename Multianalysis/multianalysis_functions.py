import glob
import os
import shutil
import subprocess
from command_dirs import command_dir
from Multianalysis.Analysis.energy import calculate_energy, read_energy
from Multianalysis.Analysis.gyration import calculate_gyration, read_gyration
from Multianalysis.Analysis.h_bonds_global_inter import calculate_h_bonds_global_inter, read_h_bonds_global_inter
from Multianalysis.Analysis.h_bonds_global_intra import calculate_h_bonds_global_intra, read_h_bonds_global_intra
from Multianalysis.Analysis.h_bonds_local_inter import calculate_h_bonds_local_inter, read_h_bonds_local_inter
from Multianalysis.Analysis.h_bonds_local_intra import calculate_h_bonds_local_intra, read_h_bonds_local_intra
from Multianalysis.Analysis.native_contacts_global import calculate_native_contacts_global, read_native_contacts_global
from Multianalysis.Analysis.native_contacts_local import calculate_native_contacts_local, read_native_contacts_local
from Multianalysis.Analysis.pressure import calculate_pressure, read_pressure
from Multianalysis.Analysis.rmsd import calculate_rmsd, read_rmsd
from Multianalysis.Analysis.rmsdist import calculate_rmsdist, read_rmsdist
from Multianalysis.Analysis.rmsf import calculate_rmsf, read_rmsf
from Multianalysis.Analysis.sasa import calculate_sasa, read_sasa, read_sasa_local
from Multianalysis.Analysis.ss import calculate_ss, read_ss, read_ss_coil
from Multianalysis.Analysis.TM_score import calculate_tm_score, read_tm_score
from Multianalysis.Plotting.plot_global import plot_global_simple
from Multianalysis.Plotting.plot_local import plot_local_simple, plot_local_compared
gmx = command_dir["gromacs"]


def copy_input_files(input_dir, save_dir):
    """
    Prepares a folder with the original input files in the output directory for possible replication

    :param input_dir: Directory in which the input folders are located
    :param save_dir: Path of the folder in the output directory where the input files are copied
    :return: Copies the edr, xtc and tpr file from the input to the output folder
    """
    edr = glob.glob(input_dir + "*.edr")
    xtc = glob.glob(input_dir + "*.xtc")
    tpr = glob.glob(input_dir + "*.tpr")
    gro = glob.glob(input_dir + "*.gro")
    for filetype in [edr, xtc, tpr, gro]:
        for item in filetype:
            item_type = item[-4:]
            if not glob.glob(save_dir + "*%s" % item_type):
                shutil.copy2(item, save_dir)


def global_multianalysis(input_dir, output_dir, sample, replicas, time_step, final, starting, protein_type, label):
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
    specific_output_dir = output_dir + sample + "/Multianalysis/Raw_results/"
    specific_output_plots = output_dir + sample + "/Multianalysis/Plots/"
    sample_data_dict = {}
    for replica in replicas:
        replica_data_dict = {}
        specific_in_dir = input_dir + sample + "/" + replica + "/"
        input_save = output_dir + "Original_input/" + sample + "/" + replica + "/"

        # Copy input to the resulting folder
        copy_input_files(specific_in_dir, input_save)

        # Generating .gro for reference
        if not glob.glob(specific_in_dir + "*.gro"):
            subprocess.call("echo 'System' | %s trjconv -f %s*.xtc -o %sreference.gro -s %s*.tpr "
                            "-dump %i -pbc nojump" % (gmx, specific_in_dir, specific_in_dir, specific_in_dir,
                                                      int(starting * 1000)),
                            shell=True)
        ini_gro = glob.glob(specific_in_dir + "*.gro")[0]
        if not glob.glob(specific_in_dir + "Multianalysis/*.gro"):
            shutil.copy2(ini_gro, specific_in_dir + "Multianalysis/")
        # Calculating all global functions
        calculate_energy(specific_in_dir, specific_output_dir, replica)
        calculate_gyration(specific_in_dir, specific_output_dir, replica)
        calculate_h_bonds_global_inter(specific_in_dir, specific_output_dir, replica)
        calculate_h_bonds_global_intra(specific_in_dir, specific_output_dir, replica)
        calculate_native_contacts_global(specific_in_dir, specific_output_dir, replica)
        calculate_pressure(specific_in_dir, specific_output_dir, replica)
        calculate_rmsd(specific_in_dir, specific_output_dir, replica)
        calculate_rmsdist(specific_in_dir, specific_output_dir, replica)
        calculate_rmsf(specific_in_dir, specific_output_dir, replica)  # Includes B-factors
        calculate_sasa(specific_in_dir, specific_output_dir, replica)
        calculate_ss(specific_in_dir, specific_output_dir, replica, protein_type)
        calculate_tm_score(specific_in_dir, specific_output_dir, replica, starting, time_step, final)

        # Reading all global functions
        replica_data_dict["energy"] = read_energy(specific_output_dir, replica)
        replica_data_dict["gyration"] = read_gyration(specific_output_dir, replica)
        replica_data_dict["hbonds_inter"] = read_h_bonds_global_inter(specific_output_dir, replica)
        replica_data_dict["hbonds_intra"] = read_h_bonds_global_intra(specific_output_dir, replica)
        replica_data_dict["native_contacts"] = read_native_contacts_global(specific_output_dir, replica)
        replica_data_dict["pressure"] = read_pressure(specific_output_dir, replica)
        replica_data_dict["rmsd"] = read_rmsd(specific_output_dir, replica)
        replica_data_dict["rmsdist"] = read_rmsdist(specific_output_dir, replica)
        replica_data_dict["rmsf"] = read_rmsf(specific_output_dir, replica)
        replica_data_dict["sasa"] = read_sasa(specific_output_dir, replica)
        replica_data_dict["ss"] = read_ss(specific_output_dir, replica)
        replica_data_dict["coil"] = read_ss_coil(specific_output_dir, replica)
        replica_data_dict["tm_score"] = read_tm_score(specific_output_dir, replica)

        sample_data_dict[replica] = replica_data_dict

    # Plotting function
    plot_global_simple(sample_data_dict, label, specific_output_plots, protein_type)

    # Returns all data for comparison
    return sample_data_dict


def local_multianalysis(input_dir, output_dir, reference_sample, sample, reference_replicas, replicas, labels_dict,
                        residue):
    """
    Calculates all local metrics around the residue of interest for the sample and the reference, and plots the single
    and compared metrics for both samples.

    :param input_dir: Directory in which the input folders are located
    :param output_dir: Directory in which the results will be saved
    :param reference_sample: Sample to which all the others are compared (usually WT)
    :param sample: Sample on which the calculation is to be performed
    :param reference_replicas: Replicas for the reference sample
    :param replicas: Replicas for the sample
    :param labels_dict: Dictionary with the label for each sample
    :param residue: Residue of interest (usually the mutated one)
    :return: Saves single and compared local plots for both samples
    """
    all_data_dict = {}
    # Calculating reference
    specific_ref_output_dir = output_dir + reference_sample + "/Multianalysis/Raw_results/"
    specific_ref_output_plots = output_dir + reference_sample + "/Multianalysis/Plots/"
    sample_data_dict = {}
    for replica in reference_replicas:
        replica_data_dict = {}
        specific_ref_in_dir = input_dir + "/" + reference_sample + "/" + replica + "/"
        calculate_h_bonds_local_inter(specific_ref_in_dir, specific_ref_output_dir, replica, residue)
        calculate_h_bonds_local_intra(specific_ref_in_dir, specific_ref_output_dir, replica, residue)
        calculate_native_contacts_local(specific_ref_in_dir, specific_ref_output_dir, replica, residue)
        calculate_sasa(specific_ref_in_dir, specific_ref_output_dir, replica)
        replica_data_dict["hbonds_inter_l"] = read_h_bonds_local_inter(specific_ref_output_dir, replica, residue)
        replica_data_dict["hbonds_intra_l"] = read_h_bonds_local_intra(specific_ref_output_dir, replica, residue)
        replica_data_dict["native_contacts_l"] = read_native_contacts_local(specific_ref_output_dir, replica, residue)
        replica_data_dict["sasa_l"] = read_sasa_local(specific_ref_output_dir, replica, residue)
        sample_data_dict[replica] = replica_data_dict
    all_data_dict[reference_sample] = sample_data_dict

    # Plotting reference
    plot_local_simple(all_data_dict[reference_sample], labels_dict[reference_sample], specific_ref_output_plots, residue)

    # Calculating sample
    specific_output_dir = output_dir + sample + "/Multianalysis/Raw_results/"
    specific_output_plots = output_dir + sample + "/Multianalysis/Plots/"
    sample_data_dict = {}
    for replica in replicas:
        replica_data_dict = {}
        specific_in_dir = input_dir + "/" + sample + "/" + replica + "/"
        calculate_h_bonds_local_inter(specific_in_dir, specific_output_dir, replica, residue)
        calculate_h_bonds_local_intra(specific_in_dir, specific_output_dir, replica, residue)
        calculate_native_contacts_local(specific_in_dir, specific_output_dir, replica, residue)
        calculate_sasa(specific_in_dir, specific_output_dir, replica)
        replica_data_dict["hbonds_inter_l"] = read_h_bonds_local_inter(specific_output_dir, replica, residue)
        replica_data_dict["hbonds_intra_l"] = read_h_bonds_local_intra(specific_output_dir, replica, residue)
        replica_data_dict["native_contacts_l"] = read_native_contacts_local(specific_output_dir, replica, residue)
        replica_data_dict["sasa_l"] = read_sasa_local(specific_output_dir, replica, residue)
        sample_data_dict[replica] = replica_data_dict
    all_data_dict[sample] = sample_data_dict

    # Plotting sample
    plot_local_simple(all_data_dict[sample], labels_dict[sample], specific_output_plots, residue)

    # Plotting compared
    samples = [reference_sample, sample]
    plot_local_compared(all_data_dict, labels_dict, samples, output_dir, residue)

    return all_data_dict
