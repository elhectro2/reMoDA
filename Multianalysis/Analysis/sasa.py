import os
import glob
import mdtraj as md
import numpy as np
from Multianalysis.Analysis.rmsf import calculate_rmsf, read_rmsf


def calculate_sasa(input_dir, output_dir, replica):
    calculate_rmsf(input_dir, output_dir, replica)
    rmsf = read_rmsf(output_dir, replica)
    if not os.path.exists("%ssasa_%s.txt" % (output_dir, replica)):
        traj_dir = glob.glob(input_dir + "*.xtc")
        top_dir = glob.glob(input_dir + "Multianalysis/*.gro")
        traj = md.load_xtc(traj_dir[0], top=top_dir[0])
        trajectory = traj.remove_solvent()
        time_list = np.divide(trajectory.time, 1000).tolist()
        txt_time_list = []
        for element in time_list:
            txt_time_list.append(str(round(element, 3)))
        include = []
        for atom in trajectory.topology.atoms:
            if atom.element.symbol != "VS":
                include.append(atom.index)
        traj = trajectory.atom_slice(include)
        print("sasa calculation about to start")
        sasa = md.shrake_rupley(traj)
        print("sasa calculation finished")
        total_sasa = sasa.sum(axis=1)
        sasa_file = open("%ssasa_%s.txt" % (output_dir, replica), "w")
        data = []
        for line in sasa:
            data.append(line.tolist())
        sasa_file.write("Time(ps) Total_SASA resid%s\n" % ' resid'.join(str(x) for x in rmsf[0])) # TODO revisar esto cuando esté read_rmsf
        for j in range(len(txt_time_list)):
            sasa_file.write(txt_time_list[j] + " " + str(total_sasa[j]))
            sasa_file.write(" %s\n" % ' '.join(str(x) for x in data[j]))
        sasa_file.close()


def read_sasa(output_dir, replica):
    file_dir = "%ssasa_%s.txt" % (output_dir, replica)
    data = []
    time_series = []
    values = []
    file = open(file_dir)
    for line in file:  # Reads all the file and does things
        if line[0] != "T":  # Non-comment lines
            line = line.split()  # Splits the line in a list of "words" (elements separated by spaces in a string)
            data.append(line)  # Data for each time is stored in a list of lists
    for item in data:
        time_series.append(float(item[0]))
        values.append(float(item[1]))
    return [time_series, values]


def read_sasa_local(output_dir, replica, residue):
    file_dir = "%ssasa_%s.txt" % (output_dir, replica)
    data = []
    time_series = []
    values = []
    file = open(file_dir)
    for line in file:  # Reads all the file and does things
        if line[0] != "T":  # Non-comment lines
            line = line.split()  # Splits the line in a list of "words" (elements separated by spaces in a string)
            data.append(line)  # Data for each time is stored in a list of lists
        else:
            line = line.split()
            index = line.index("resid%i" % residue)
    for item in data:
        time_series.append(float(item[0]))
        values.append(float(item[index]))
    return [time_series, values]
