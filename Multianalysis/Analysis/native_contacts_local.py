import numpy as np
import mdtraj as md
from itertools import combinations
import os
import shutil
import glob
from Multianalysis.Analysis.rmsf import calculate_rmsf, read_rmsf


def best_hummer_q(traj, native, resid, start):
    """Compute the fraction of native contacts according the definition from
    Best, Hummer and Eaton [1], only if they involve the residue of interest.

    Parameters
    ----------
    traj : md.Trajectory
        The trajectory to do the computation for
    native : md.Trajectory
        The 'native state'. This can be an entire trajectory, or just a single frame.
        Only the first conformation is used

    Returns
    -------
    q : np.array, shape=(len(traj),)
        The fraction of native contacts in each frame of `traj`

    References
    ----------
    ..[1] Best, Hummer, and Eaton, "Native contacts determine protein folding
          mechanisms in atomistic simulations" PNAS (2013)
    """

    BETA_CONST = 50  # 1/nm
    LAMBDA_CONST = 1.8
    NATIVE_CUTOFF = 0.45  # nanometers

    # get the indices of all of the heavy atoms
    heavy = native.topology.select_atom_indices('heavy')
    # get the pairs of heavy atoms which are farther than 3
    # residues apart
    heavy_pairs = np.array(
        [(i, j) for (i, j) in combinations(heavy, 2)
         if (abs(native.topology.atom(i).residue.index - native.topology.atom(j).residue.index) > 3) and
         (native.topology.atom(i).residue.index == resid - start or native.topology.atom(j).residue.index == resid - start)])
    # compute the distances between these pairs in the native state
    heavy_pairs_distances = md.compute_distances(native[0], heavy_pairs)[0]
    # and get the pairs s.t. the distance is less than NATIVE_CUTOFF
    native_contacts = heavy_pairs[heavy_pairs_distances < NATIVE_CUTOFF]

    # now compute these distances for the whole trajectory
    r = md.compute_distances(traj, native_contacts)
    # and recompute them for just the native state
    r0 = md.compute_distances(native[0], native_contacts)

    q = np.mean(1.0 / (1 + np.exp(BETA_CONST * (r - LAMBDA_CONST * r0))), axis=1)
    return q


def calculate_native_contacts_local(input_dir, output_dir, replica, residue):
    if not os.path.isfile("%snative_contacts_resid%s_%s.txt" % (output_dir, str(residue), replica)):
        calculate_rmsf(input_dir, output_dir, replica)
        start = int(read_rmsf(output_dir, replica)[0][0])
        traj_dir = glob.glob(input_dir + "*.xtc")
        top_dir = glob.glob(input_dir + "Multianalysis/*.gro")
        traj = md.load_xtc(traj_dir[0], top=top_dir[0])
        nat = md.load(top_dir[0])
        q = best_hummer_q(traj, nat, residue, start)
        out_file = open("%snative_contacts_resid%s_%s.txt" % (output_dir, str(residue), replica), "w")
        out_file.write("Time(ns) fraction_of_native_contacts\n")
        for j in range(len(q)):
            out_file.write("%f %f\n" % (round(float(traj.time[j])/1000, 3), q[j]))
        out_file.close()


def read_native_contacts_local(output_dir, replica, residue):
    file_dir = "%snative_contacts_resid%s_%s.txt" % (output_dir, str(residue), replica)
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
