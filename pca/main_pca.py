from pca.pca_functions import pca_local, pca_global


def main_pca(input_dir, output_dir, samples_replicas_dict, reference_sample, labels_dict, position_dict, time_step,
             final, starting, protein_type, parts=20, video_frames=100):
    """
    Makes all calculations and plots of the PCA module

    :param input_dir: Directory in which the input folders are located
    :param output_dir: Directory in which the results will be saved
    :param samples_replicas_dict: Dictionary with all samples and their respective replicas
    :param reference_sample: Sample to which all the others are compared (usually WT)
    :param labels_dict: Dictionary with the label for each sample
    :param position_dict: Dictionary with the mutation position for each sample
    :param time_step: Time step between frames of the input .xtc trajectory (in ns)
    :param final: Ending time of the simulation (in ns)
    :param starting: Starting time of the simulations (in ns, usually 0)
    :param protein_type: Folding according to secondary structure. Expected inputs are a/alpha, b/beta or a+b/alpha+beta
    :param parts: Number of parts in which to divide the trajectory for the parts plot
    :param video_frames: The number of frames that are generated for making the video

    :return: None
    """
    global_data = pca_global(input_dir, output_dir, samples_replicas_dict, reference_sample, labels_dict, time_step,
                             final, starting, protein_type, parts=parts, video_frames=video_frames)
    pca_local(global_data, input_dir, output_dir, samples_replicas_dict, reference_sample, labels_dict, position_dict,
              time_step, parts=parts, video_frames=video_frames)
