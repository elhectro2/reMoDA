import os
import subprocess
from command_dirs import command_dir
gmx = command_dir["gromacs"]


def calculate_energy(input_dir, output_dir, replica):
    if not os.path.exists("%senergy_%s.xvg" % (output_dir, replica)):
        order = "echo 'Total-Energy' | %s energy -f %s*.edr -o %senergy_%s.xvg" % (gmx, input_dir, output_dir, replica)
        subprocess.call(order, shell=True)


def read_energy(output_dir, replica):
    file_dir = "%senergy_%s.xvg" % (output_dir, replica)
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
