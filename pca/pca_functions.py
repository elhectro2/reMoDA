from Multianalysis.Plotting.colors import colors_replicas
import cv2
from helper_functions import  check_makedirs
import matplotlib.pyplot as plt
import matplotlib
from Multianalysis.multianalysis_functions import global_multianalysis, local_multianalysis
import os
from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
global_pca_analysis_functions = ["gyration", "hbonds_inter", "hbonds_intra", "native_contacts", "rmsd",
                                 "rmsdist", "sasa", "ss", "coil", "tm_score"]
local_pca_analysis_functions = ["hbonds_inter_l", "hbonds_intra_l", "native_contacts_l", "sasa_l"]
os.environ["XDG_SESSION_TYPE"] = "minimal"  # TODO check if it is not necessary


def get_all_data_global(input_dir, output_dir, samples_replicas_dict, reference_sample, labels_dict,
                        time_step, final, starting, protein_type):
    """
    Gets all data for the global analysis as a dictionary.

    :param input_dir: Directory in which the input folders are located
    :param output_dir: Directory in which the results will be saved
    :param samples_replicas_dict: Dictionary with all samples and their respective replicas
    :param reference_sample: Sample to which all the others are compared (usually WT)
    :param labels_dict: Dictionary with the label for each sample
    :param time_step: Time step between frames of the input .xtc trajectory (in ns)
    :param final: Ending time of the simulation (in ns)
    :param starting: Starting time of the simulations (in ns, usually 0)
    :param protein_type: Folding according to secondary structure. Expected inputs are a/alpha, b/beta or a+b/alpha+beta

    :return: A dictionary with all data resulting from global analysis on all replicas of all samples
    """

    global_data_dict = {}
    ref_label = labels_dict[reference_sample]
    global_data_dict[reference_sample] = global_multianalysis(input_dir, output_dir, reference_sample,
                                                           samples_replicas_dict[reference_sample],
                                                           time_step, final, starting, protein_type, ref_label)
    for sample in samples_replicas_dict:
        if sample != reference_sample:
            sample_label = labels_dict[sample]
            global_data_dict[sample] = global_multianalysis(input_dir, output_dir, sample, samples_replicas_dict[sample],
                                                         time_step, final, starting, protein_type, sample_label)
    return global_data_dict


def get_additional_data_local(input_dir, output_dir, samples_replicas_dict, sample, reference_sample, labels_dict,
                              position):
    """
    Gets the data for the local properties of a given residue of the sample of interest and the reference sample.

    :param input_dir: Directory in which the input folders are located
    :param output_dir: Directory in which the results will be saved
    :param samples_replicas_dict: Dictionary with all samples and their respective replicas
    :param sample: Sample on which the calculation is to be performed
    :param reference_sample: Sample of reference for the analysis
    :param labels_dict: Dictionary with the label for each sample
    :param position: Position of the residue of interest

    :return: A dictionary with all data resulting from local analysis on the residue of interest on a sample and
    the reference sample.
    """
    reference_replicas = samples_replicas_dict[reference_sample]
    replicas = samples_replicas_dict[sample]
    data_dict_local = local_multianalysis(input_dir, output_dir, reference_sample, sample, reference_replicas,
                                              replicas, labels_dict, position)
    return data_dict_local


def get_pca_model(data_dict, samples_replicas_dict, out_pca, include_local=False):
    """
    Obtains data for PCA for all data.

    :param data_dict: A dictionary with all data resulting from (global or local) multianalysis on all replicas of all
    samples
    :param samples_replicas_dict: Dictionary with all samples and their respective replicas
    :param out_pca: Directory in which to save the statistics and weights of the PCA model

    :return: pc_per_replica: Dictionary of principal components 1 and 2 for each replica to be analyzed
    :return: total_xlim: Minimum and maximum value for PC1, as a duple
    :return: total_ylim: Minimum and maximum value for PC2, as a duple
    """
    # Organizing data from the original format to a PCA-compatible dictionary, with canonic times.
    canonic_x = []
    y_dict = {}
    analysis_functions = global_pca_analysis_functions
    if include_local:
        analysis_functions += local_pca_analysis_functions
    for sample in samples_replicas_dict:
        if sample not in y_dict:
            y_dict[sample] = {}
        for replica in data_dict[sample]:
            if replica not in y_dict[sample]:
                y_dict[sample][replica] = {}
            for function in analysis_functions:
                if not canonic_x:
                    canonic_x = [float(round(i, 2)) for i in data_dict[sample][replica][function][0]]
                new_y_data = ["NA" for i in canonic_x]
                for i in range(len(data_dict[sample][replica][function][0])):
                    time = float(round(data_dict[sample][replica][function][0][i], 2))
                    data_point = data_dict[sample][replica][function][1][i]
                    if time in canonic_x:
                        list_index = canonic_x.index(time)
                        new_y_data[list_index] = data_point
                y_dict[sample][replica][function] = new_y_data

    # Getting all data in the same list for passing it to the PCA fitting algorithm
    y_per_replica = []
    y = []
    for sample in y_dict:
        for replica in y_dict[sample]:
            name = sample + "_" + replica
            for i in range(len(canonic_x)):
                y_at_time = []
                for function in global_pca_analysis_functions:
                    y_at_time.append(y_dict[sample][replica][function][i])
                y.append(y_at_time)
            y_per_replica.append(name)

    # Fitting the PCA algorithm
    pca_pipe = make_pipeline(StandardScaler(), PCA())
    pca_pipe.fit(y)
    pca_model = pca_pipe.named_steps['pca']

    # Getting PCA weights and explained variance
    pc_weights_header = global_pca_analysis_functions
    # check_makedirs(out_pca)
    file_out = open(f"{out_pca}PC_stats.txt", "w")
    file_out.write("PC  %variance " + " ".join(pc_weights_header) + "\n")
    comps = pca_model.components_
    variances = pca_model.explained_variance_ratio_
    for i in range(len(variances)):
        pc = f"PC{(i + 1)} "
        variance = str(round(variances[i] * 100, 2)) + "% "
        comp = [str(j) for j in comps[i]]
        file_out.write(pc + variance + " ".join(comp) + "\n")
    file_out.close()

    # Getting the PC values in a single file.
    x = canonic_x * len(y_per_replica)
    pc_data = pca_pipe.transform(y)
    file_out = open(f"{out_pca}PC_values.txt", "w")
    file_out.write("IID " + " ".join("PC%i" % (j + 1) for j in range(len(pc_data[0]))) + "\n")
    for i in range(len(pc_data)):
        file_out.write("%s %s\n" % (x[i], " ".join([str(j) for j in pc_data[i]])))
    file_out.close()

    # Getting the PC values in a dictionary for further use
    pc_per_replica = {}
    i = 0
    for sample in y_dict:
        for replica in y_dict[sample]:
            pc_replica = []
            name = sample + "_" + replica
            for item in range(len(canonic_x)):
                pc_at_time = [x[i], pc_data[i][0], pc_data[i][1]]
                pc_replica.append(pc_at_time)
                i += 1
            pc_per_replica[name] = pc_replica
    total_xlim = [min(pc_data[:, 0]), max(pc_data[:, 0])]
    total_ylim = [min(pc_data[:, 1]), max(pc_data[:, 1])]
    return pc_per_replica, total_xlim, total_ylim


def plot_pca_parts(r_data, total_xlim, total_ylim, parts_output_dir, parts=20, time_step=0.1):
    """
    Obtains plots for a full trajectory, divided in a given number of parts of equal time length

    :param r_data: Dictionary of principal components 1 and 2 for each replica to be analyzed
    :param total_xlim: Minimum and maximum value for PC1, as a duple
    :param total_ylim: Minimum and maximum value for PC2, as a duple
    :param parts_output_dir: Directory in which to save the plots of each part
    :param parts: Number of parts in which to divide the trajectory
    :param time_step: Time step between frames of the input .xtc trajectory (in ns)

    :return: None
    """
    color_cycle = colors_replicas(len(r_data))
    color_cycle = [next(color_cycle) for i in range(len(r_data))]
    # check_makedirs(f"{parts_output_dir}")
    # Get individual plots for each frame of the video
    for part in range(parts):
        if part % 10 == 0:
            print(f"Part {part}")
        n_replicas = len(r_data)
        i = 0
        for replica in r_data:
            name = replica
            x = [r_data[replica][j][1] for j in range(len(r_data[replica]))]
            y = [r_data[replica][j][2] for j in range(len(r_data[replica]))]
            time = [r_data[replica][j][0] for j in range(len(r_data[replica]))]
            alpha_factor = 1 - (i / (n_replicas - 1)) * 0.3

            start_prev = int(round((part - 1) * len(y) / parts, 0))
            start = int(part * len(y) / parts)
            end = (part + 1) * len(y) / parts
            if int(end) != round(end):
                end = int(end) + 1
            else:
                end = int(end)
            start_prev_bound = 0

            if start_prev >= 0:
                old_x = x[start_prev:start]
                old_y = y[start_prev:start]
                start_prev_bound = start_prev
                plt.scatter(old_x, old_y, alpha=0.5 * alpha_factor, edgecolors="none",
                            c=color_cycle[i % len(color_cycle)], s=15)

            new_x = x[start:end]
            new_y = y[start:end]
            plt.scatter(new_x, new_y, label=name, alpha=alpha_factor, edgecolors="none",
                        c=color_cycle[i % len(color_cycle)], s=15)
            i += 1
        bar_cmap = (matplotlib.colors.ListedColormap(["white", ("grey", 0.5), "grey", "white"]))
        bounds = [0, start_prev_bound * time_step, start * time_step, end * time_step, max(time) + 1]
        print("bounds " + str(bounds))
        norm = matplotlib.colors.BoundaryNorm(bounds, bar_cmap.N)
        plt.colorbar(matplotlib.cm.ScalarMappable(cmap=bar_cmap, norm=norm), location="top", spacing="proportional",
                     label="Time(ns)", ax=plt.gca(), ticks=[(max(time) / 5) * i for i in range(6)])
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.legend(loc="best")
        plt.xlim(total_xlim)
        plt.ylim(total_ylim)
        plt.savefig(f"{parts_output_dir}/part{part + 1}.png")
        plt.close()


def plot_pca_all(pc_per_replica, total_xlim, total_ylim, all_output_dir, time_step=0.1):
    """
    Obtains the PCA plot for a full trajectory

    :param pc_per_replica: Dictionary of principal components 1 and 2 for each replica to be analyzed
    :param total_xlim: Minimum and maximum value for PC1, as a duple
    :param total_ylim: Minimum and maximum value for PC2, as a duple
    :param all_output_dir: Directory in which to save the plots of the whole trajectory
    :param time_step: Time step between frames of the input .xtc trajectory (in ns)

    :return: None
    """
    plot_pca_parts(pc_per_replica, total_xlim, total_ylim, all_output_dir, parts=1, time_step=time_step)


def plot_pca_video(frames_dir, video_output_dir, video_frames=100):
    """
    Obtains a video from a sequence of plots located in a folder

    :param frames_dir: The directory in which the plots for the video are located
    :param video_output_dir: The directory in which the final video is saved
    :param video_frames: The number of frames that are available for the video

    :return: None
    """
    img_folder = f"{frames_dir}"
    video_name = f'{video_output_dir}video_PCA.avi'

    images = [f"{img_folder}/part{part + 1}.png" for part in range(video_frames)]
    print(images)
    frame = cv2.imread(images[0])
    height, width, layers = frame.shape

    video = cv2.VideoWriter(video_name, 0, 10, (width, height))

    for image in images:
        video.write(cv2.imread(image))

    cv2.destroyAllWindows()
    video.release()


def pca_global(input_dir, output_dir, samples_replicas_dict, reference_sample, labels_dict, time_step, final, starting,
               protein_type, parts=20, video_frames=100):
    """
    The full function for calculating all PCA plots with global properties from multianalysis

    :param input_dir: Directory in which the input folders are located
    :param output_dir: Directory in which the results will be saved
    :param samples_replicas_dict: Dictionary with all samples and their respective replicas
    :param reference_sample: Sample to which all the others are compared (usually WT)
    :param labels_dict: Dictionary with the label for each sample
    :param time_step: Time step between frames of the input .xtc trajectory (in ns)
    :param final: Ending time of the simulation (in ns)
    :param starting: Starting time of the simulations (in ns, usually 0)
    :param protein_type: Folding according to secondary structure. Expected inputs are a/alpha, b/beta or a+b/alpha+beta
    :param parts: Number of parts in which to divide the trajectory for the parts plot
    :param video_frames: The number of frames that are generated for making the video

    :return: A dictionary with all data resulting from global analysis on all replicas of all samples
    """
    print("get all data")
    global_data_dict = get_all_data_global(input_dir, output_dir, samples_replicas_dict, reference_sample, labels_dict,
                                        time_step, final, starting, protein_type)
    print("get pca model")
    output_pca_data = f"{output_dir}All/PCA/Global/"
    pc_per_replica, global_xlim, global_ylim = get_pca_model(global_data_dict, samples_replicas_dict, output_pca_data)
    print("get full PCA plot")
    all_global_output_dir = f"{output_dir}All/PCA/Global/all/"
    plot_pca_all(pc_per_replica, global_xlim, global_ylim, all_global_output_dir, time_step=0.1)
    print("get parts")
    parts_output_dir = f"{output_dir}All/PCA/Global/parts/"
    plot_pca_parts(pc_per_replica, global_xlim, global_ylim, parts_output_dir, parts=parts, time_step=0.1)
    print("get video")
    if parts != video_frames:
        frames_dir = f"{output_dir}All/PCA/Global/video_frames"
        plot_pca_parts(pc_per_replica, global_xlim, global_ylim, frames_dir, parts=video_frames,
                              time_step=time_step)
    else:
        frames_dir = parts_output_dir
    video_output_dir = f"{output_dir}All/PCA/Global/"
    plot_pca_video(frames_dir, video_output_dir, video_frames=video_frames)
    return global_data_dict


def pca_local(global_data_dict, input_dir, output_dir, samples_replicas_dict, reference_sample, labels_dict,
              position_dict, time_step, parts=20, video_frames=100):
    """

    :param global_data_dict: A dictionary with all data resulting from global analysis on all replicas of all samples
    :param input_dir: Directory in which the input folders are located
    :param output_dir: Directory in which the results will be saved
    :param samples_replicas_dict: Dictionary with all samples and their respective replicas
    :param reference_sample: Sample to which all the others are compared (usually WT)
    :param labels_dict: Dictionary with the label for each sample
    :param position_dict: Dictionary with the mutation position for each sample
    :param time_step: Time step between frames of the input .xtc trajectory (in ns)
    :param parts: Number of parts in which to divide the trajectory for the parts plot
    :param video_frames: The number of frames that are generated for making the video
    :return:
    """
    for sample in global_data_dict:
        if sample != reference_sample and position_dict[sample]:
            print(f"Get local data {labels_dict[sample]}")
            local_data_dict = {sample: global_data_dict[sample], reference_sample: global_data_dict[reference_sample]}
            new_local_data = local_multianalysis(input_dir, output_dir, reference_sample, sample,
                                                 samples_replicas_dict[reference_sample], samples_replicas_dict[sample],
                                                 labels_dict, position_dict[sample])
            for item in new_local_data:
                for replica in new_local_data[item]:
                    local_data_dict[item][replica].update(new_local_data[item][replica])
            print(f"get local pca model {labels_dict[sample]}")
            output_pca_data = f"{output_dir}{labels_dict[reference_sample]}_vs_{labels_dict[sample]}/PCA/Local/"
            pc_per_replica, local_xlim, local_ylim = get_pca_model(local_data_dict, samples_replicas_dict,
                                                                   output_pca_data, include_local=True)
            print(f"get full local PCA plot {labels_dict[sample]}")
            all_local_output_dir = (f"{output_dir}{labels_dict[reference_sample]}_vs_{labels_dict[sample]}/PCA/Local/"
                                    f"all/")
            plot_pca_all(pc_per_replica, local_xlim, local_ylim, all_local_output_dir, time_step=0.1)
            print(f"get local parts {labels_dict[sample]}")
            local_parts_output_dir = (f"{output_dir}{labels_dict[reference_sample]}_vs_{labels_dict[sample]}/PCA/Local/"
                                      f"parts/")
            plot_pca_parts(pc_per_replica, local_xlim, local_ylim, local_parts_output_dir, parts=parts, time_step=0.1)
            print(f"get local video {labels_dict[sample]}")
            if parts != video_frames:
                local_frames_dir = (f"{output_dir}{labels_dict[reference_sample]}_vs_{labels_dict[sample]}/PCA/Local/"
                                    f"video_frames")
                plot_pca_parts(pc_per_replica, local_xlim, local_ylim, local_frames_dir, parts=video_frames,
                               time_step=time_step)
            else:
                local_frames_dir = local_parts_output_dir
            video_output_dir = f"{output_dir}{labels_dict[reference_sample]}_vs_{labels_dict[sample]}/PCA/Local/"
            plot_pca_video(local_frames_dir, video_output_dir, video_frames=video_frames)
