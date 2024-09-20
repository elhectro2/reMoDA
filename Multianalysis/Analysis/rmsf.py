import os
import subprocess
from command_dirs import command_dir
gmx = command_dir["gromacs"]


def calculate_rmsf(input_dir, output_dir, replica):
    if not os.path.exists("%srmsf_%s.xvg" % (output_dir, replica)):
        order = "echo 'Backbone' | %s rmsf -f %s*.xtc -s %s*.tpr -o %srmsf_%s.xvg " \
                "-ox %sBfactors_on_averaged_structure_%s.pdb -res" % (gmx, input_dir, input_dir, output_dir, replica,
                                                                      output_dir, replica)
        subprocess.call(order, shell=True)


def read_rmsf(output_dir, replica):
    file_dir = "%srmsf_%s.xvg" % (output_dir, replica)
    data = []
    res_series = []
    values = []
    file = open(file_dir)
    for line in file:  # Reads all the file and does things
        if line[0] not in "#@":  # Non-comment lines
            line = line.split()  # Splits the line in a list of "words" (elements separated by spaces in a string)
            data.append(line)  # Data for each time is stored in a list of lists
    for item in data:
        res_series.append(int(item[0]))
        values.append(float(item[1]))
    return [res_series, values]
