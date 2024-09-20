from matplotlib import pyplot as plt
from Multianalysis.Plotting.colors import colors_replicas

energetics_dict = {"proper": {"title": "Proper dihedrals", "y_label": "(kJ/mol)", "x_label": "Time (ns)"},
                   "LJ_pp": {"title": "L-J Prot-Prot", "y_label": "(kJ/mol)", "x_label": "Time (ns)"},
                   "coul_pp": {"title": "Coul Prot-Prot", "y_label": "(kJ/mol)", "x_label": "Time (ns)"},
                   "improper": {"title": "Improper dihedrals", "y_label": "(kJ/mol)", "x_label": "Time (ns)"},
                   "LJ_pn": {"title": "L-J Prot-non-Prot", "y_label": "(kJ/mol)", "x_label": "Time (ns)"},
                   "coul_pn": {"title": "Coul Prot-non-Prot", "y_label": "(kJ/mol)", "x_label": "Time (ns)"}}


def plot_energetics_simple(sample_data, sample_label, output_dir, out_name):
    """
    Plots a given number of thematically related global metrics of a single sample.

    :param energetics_dict: Dictionary with the related metrics, each with their title and x and y labels
    :param sample_data: Data for all metrics for a single sample
    :param sample_label: Label for the sample
    :param output_dir: Directory in which the results will be saved
    :param protein_type: Folding according to secondary structure. Expected inputs are a/alpha, b/beta or a+b/alpha+beta
    :param out_name: Name for the thematically related plot (e.g., Interactions, Energetics...)
    :return: Saves the desired plot in the output folder
    """
    plots = len(energetics_dict)
    rows = int(plots / 3)
    fig = plt.figure(figsize=(plots / 2 * 2.5 + 2.5, rows * 2.5 + 1))
    axes = []
    sample_replicas = [x for x in sample_data]
    i = 0
    for item in energetics_dict:
        n_replicas = len(sample_replicas)
        colors = colors_replicas(n_replicas)
        ax = fig.add_subplot(rows, int(plots / 2), i + 1)
        j = 0
        for replica in sample_replicas:
            plot_data = sample_data[replica][item]
            ax.plot(plot_data[0], plot_data[1], label=replica, c=next(colors), alpha=1 - j * (0.5 / n_replicas))
            j += 1
        ax.set_title(energetics_dict[item]["title"], fontsize=18)
        ax.set_xlabel(energetics_dict[item]["x_label"])
        ax.set_ylabel(sample_label + " " + energetics_dict[item]["y_label"])
        ax.legend(loc="best")
        ax.margins(0.1)
        axes.append(ax)
        i += 1
    fig.tight_layout()
    fig.savefig(output_dir + "%s_%s.png" % (out_name, sample_label))
    plt.close(fig)


def plot_energetics_compared(all_data_dict, labels_dict, samples, output_dir, out_name):
    """
    Plots a given number of thematically related global metrics of a sample (or samples) compared to the reference one.

    :param energetics_dict: Dictionary with the related metrics, each with their title and x and y labels
    :param all_data_dict: Data for all metrics for all samples
    :param labels_dict: Dictionary with the label for each sample
    :param samples: Samples to be compared (usually reference and mutant)
    :param output_dir: Directory in which the results will be saved
    :param protein_type: Folding according to secondary structure. Expected inputs are a/alpha, b/beta or a+b/alpha+beta
    :param out_name: Name for the thematically related plot (e.g., Interactions, Energetics...)
    :return: Saves the desired plot in the output folder
    """
    sample_compared_label = "_vs_".join(samples)
    plot_out_dir = output_dir + sample_compared_label + "/Energetics/Plots/"
    plots = len(energetics_dict)
    rows = int(len(samples) * plots / 3)
    fig = plt.figure(figsize=(plots/2 * 2.5 + 2.5, rows * 2.5 + 1))
    i = 0
    axes = []
    for sample in samples:
        sample_data = all_data_dict[sample]
        sample_label = labels_dict[sample]
        sample_replicas = [x for x in sample_data]
        for item in energetics_dict:
            n_replicas = len(sample_replicas)
            colors = colors_replicas(n_replicas)
            ax = fig.add_subplot(rows, int(plots / 2), i + 1)
            j = 0
            for replica in sample_replicas:
                plot_data = sample_data[replica][item]
                ax.plot(plot_data[0], plot_data[1], label=replica, c=next(colors), alpha=1 - j * (0.5 / n_replicas))
                j += 1
            ax.set_xlabel(energetics_dict[item]["x_label"])
            ax.set_title(energetics_dict[item]["title"], fontsize=18)
            ax.set_ylabel(sample_label + " " + energetics_dict[item]["y_label"])
            ax.legend(loc="best")
            ax.margins(0.1)
            axes.append(ax)
            i += 1
    for i in range(plots):
        xmin = 9999999999999
        xmax = -9999999999999
        ymin = 9999999999999
        ymax = -9999999999999
        for j in range(len(samples)):
            xlim = axes[i + j * plots].get_xlim()
            xmin = min(xlim[0], xmin)
            xmax = max(xlim[1], xmax)
            ylim = axes[i + j * plots].get_ylim()
            ymin = min(ylim[0], ymin)
            ymax = max(ylim[1], ymax)
        for j in range(len(samples)):
            axes[i + j * plots].set_xlim(xmin, xmax)
            axes[i + j * plots].set_ylim(ymin, ymax)
            axes[i + j * plots].legend(loc="best")
    fig.tight_layout()
    fig.savefig(plot_out_dir + "%s_%s.png" % (out_name, sample_compared_label))
    plt.close(fig)
