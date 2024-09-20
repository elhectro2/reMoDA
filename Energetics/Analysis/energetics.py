import os
import subprocess
from command_dirs import command_dir
gmx = command_dir["gromacs"]


def calculate_coul_pp(input_dir, output_dir, replica):
    if not os.path.exists("%scoul_pp_%s.xvg" % (output_dir, replica)):
        order = "echo 'Coul-SR:Protein-Protein' | %s energy -f %s*.edr -o %scoul_pp_%s.xvg" % (gmx, input_dir,
                                                                                               output_dir, replica)
        subprocess.call(order, shell=True)


def calculate_coul_pn(input_dir, output_dir, replica):
    if not os.path.exists("%scoul_pn_%s.xvg" % (output_dir, replica)):
        order = "echo 'Coul-SR:Protein-non-Protein' | %s energy -f %s*.edr -o %scoul_pn_%s.xvg" % (gmx,  input_dir,
                                                                                                   output_dir, replica)
        subprocess.call(order, shell=True)


def calculate_improp(input_dir, output_dir, replica):
    if not os.path.exists("%simproper_%s.xvg" % (output_dir, replica)):
        order = "echo 'Improper-Dih.' | %s energy -f %s*.edr -o %simproper_%s.xvg" % (gmx, input_dir, output_dir,
                                                                                      replica)
        subprocess.call(order, shell=True)


def calculate_prop(input_dir, output_dir, replica):
    if not os.path.exists("%sproper_%s.xvg" % (output_dir, replica)):
        order = "echo 'Proper-Dih.' | %s energy -f %s*.edr -o %sproper_%s.xvg" % (gmx, input_dir, output_dir, replica)
        subprocess.call(order, shell=True)


def calculate_lj_pp(input_dir, output_dir, replica):
    if not os.path.exists("%sLJ_pp_%s.xvg" % (output_dir, replica)):
        order = "echo 'LJ-SR:Protein-Protein' | %s energy -f %s*.edr -o %sLJ_pp_%s.xvg" % (gmx, input_dir, output_dir,
                                                                                           replica)
        subprocess.call(order, shell=True)


def calculate_lj_pn(input_dir, output_dir, replica):
    if not os.path.exists("%sLJ_pn_%s.xvg" % (output_dir, replica)):
        order = "echo 'LJ-SR:Protein-non-Protein' | %s energy -f %s*.edr -o %sLJ_pn_%s.xvg" % (gmx,  input_dir,
                                                                                               output_dir, replica)
        subprocess.call(order, shell=True)


def read_energetics(output_dir, file_name, replica):
    file_dir = "%s%s_%s.xvg" % (output_dir, file_name, replica)
    data = []
    time_series = []
    values = []
    file = open(file_dir)
    for line in file:  # Reads all the file and does things
        if line[0] not in "#@":  # Non-comment lines
            line = line.split()  # Splits the line in a list of "words" (elements separated by spaces in a string)
            data.append(line)  # Data for each time is stored in a list of lists
    for item in data:
        time_series.append(float(item[0])/1000)
        values.append(float(item[1]))
    return [time_series, values]


def read_all_energies(output_dir, replica):
    analysis_list = ["coul_pp", "coul_pn", "improper", "proper", "LJ_pp", "LJ_pn"]
    return [read_energetics(output_dir, x, replica) for x in analysis_list]
