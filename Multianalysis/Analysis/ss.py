import os
import subprocess
from command_dirs import command_dir
gmx = command_dir["gromacs"]
dssp = command_dir["dssp"]


def calculate_ss(input_dir, output_dir, replica):
    if not os.path.exists("%ssscount_%s.xvg" % (output_dir, replica)):
        order = "export DSSP=%s\necho Protein | %s dssp -f %s*.xtc -s %s*.gro -hmode dssp" \
                "-clear -tu ns -sc %ssscount_%s.xvg -o %ssscount_%s.dat" % (dssp, gmx, input_dir, input_dir, 
                                                                             output_dir, replica, output_dir, replica)
        print(order)
        subprocess.call(order, shell=True)


def read_ss(output_dir, replica, protein_type):
    if protein_type and protein_type.lower() == "alpha":
        indices = [10]
    elif protein_type and protein_type.lower() == "beta":
        indices = [8, 9]
    else:
        indices = [8, 9, 10]
    file_dir = "%ssscount_%s.xvg" % (output_dir, replica)
    data = []
    time_series = []
    values = []
    indices = []
    file = open(file_dir, errors="ignore")
    for line in file:  # Reads all the file and does things
        print(line)
        if line[0] not in "#@":  # Non-comment lines
            line = line.split()  # Splits the line in a list of "words" (elements separated by spaces in a string)
            data.append(line)  # Data for each time is stored in a list of lists
    for item in data:
        time_series.append(float(item[0]))
        if indices:
            values.append(sum([float(item[index]) for index in indices]))
        else:
            values.append(0)
    return [time_series, values]


def read_ss_coil(output_dir, replica):
    file_dir = "%ssscount_%s.xvg" % (output_dir, replica)
    data = []
    time_series = []
    values = []
    file = open(file_dir, errors="ignore")
    for line in file:  # Reads all the file and does things
        if line[0] not in "#@":  # Non-comment lines
            line = line.split()  # Splits the line in a list of "words" (elements separated by spaces in a string)
            data.append(line)  # Data for each time is stored in a list of lists
        elif line[0:3] == "@ s" and "Loops" in line[3:]:
            index = int(line[3:].split()[0]) + 1
    print(index)
    for item in data:
        time_series.append(float(item[0]))
        values.append(float(item[1]))
    return [time_series, values]
