import os
import subprocess
from command_dirs import command_dir
gmx = command_dir["gromacs"]


def calculate_h_bonds_local_intra(input_dir, output_dir, replica, residue):
    if not os.path.isfile("%sH_bonds_resid%s_intra_%s.xvg" % (output_dir, str(residue), replica)):
        order = "%s make_ndx -f %s*.gro -o %sMultianalysis/index_%s.ndx << EOF\n1 & r %s\n1 & ! r %s\nq\nEOF" \
                % (gmx, input_dir, input_dir, str(residue), str(residue), str(residue))
        subprocess.call(order, shell=True)
        order = "echo 'Protein_&_r_%s Protein_&_!r_%s' | %s hbond -f %s*.xtc -s %s*.tpr -tu ns " \
                "-n %sMultianalysis/index_%s.ndx -num %sH_bonds_resid%s_intra_%s.xvg"\
                % (str(residue), str(residue), gmx, input_dir, input_dir, input_dir, str(residue), output_dir,
                   str(residue), replica)
        subprocess.call(order, shell=True)


def read_h_bonds_local_intra(output_dir, replica, residue):
    file_dir = "%sH_bonds_resid%s_intra_%s.xvg" % (output_dir, str(residue), replica)
    data = []
    time_series = []
    values = []
    file = open(file_dir)
    for line in file:  # Reads all the file and does things
        if line[0] not in "#@":  # Non-comment lines
            line = line.split()  # Splits the line in a list of "words" (elements separated by spaces in a string)
            data.append(line)  # Data for each time is stored in a list of lists
    for item in data:
        time_series.append(float(item[0]))
        values.append(float(item[1]))
    return [time_series, values]
