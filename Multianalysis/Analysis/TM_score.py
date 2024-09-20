import os
import subprocess
from command_dirs import command_dir
gmx = command_dir["gromacs"]
tmscore = command_dir["TMscore"]


def calculate_tm_score(input_dir, output_dir, replica, starting, time_step, final):
    if not os.path.exists("%sTMscore_%s.xvg" % (output_dir, replica)):
        if not os.path.exists("%sMultianalysis/pdbs" % input_dir):
            os.makedirs("%sMultianalysis/pdbs" % input_dir)
        order = "echo 1 | %s trjconv -f %s*.xtc -s %s*.gro -o %sMultianalysis/pdbs/pdb_.pdb -sep"\
                % (gmx, input_dir, input_dir, input_dir)
        subprocess.call(order, shell=True)
        file = open("%sTMscore_%s.xvg" % (output_dir, replica), "w")
        file.write("# Time(ns) TM-score\n%s 1.000\n" % str(starting))
        for i in range(int((final - starting) / time_step)):
            order = "TMscore %sMultianalysis/pdbs/*_%i.pdb %sMultianalysis/pdbs/*_0.pdb" % (input_dir, i + 1, input_dir)
            temp_tm = subprocess.check_output(order, shell=True, universal_newlines=True)
            temp_tm = temp_tm.split("\n")
            for line in temp_tm:
                line = line.split()
                if len(line) > 0 and line[0] == "TM-score":
                    file.write(str(round(float(starting + (i + 1) * time_step), 1)) + " " + line[2] + "\n")
        file.close()


def read_tm_score(output_dir, replica):
    file_dir = "%sTMscore_%s.xvg" % (output_dir, replica)
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
