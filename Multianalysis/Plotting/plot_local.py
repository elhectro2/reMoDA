from matplotlib import pyplot as plt
from Multianalysis.Plotting.colors import colors_replicas

interactions_dict = {"hbonds_inter_l": {"title": "Intra protein H-bonds", "y_label": "H-bonds", "x_label": "Time (ns)"},
                     "hbonds_intra_l": {"title": "Protein-water H-bonds", "y_label": "H-bonds", "x_label": "Time (ns)"},
                     "native_contacts_l": {"title": "Native contacts", "y_label": "Fraction of native contacts", "x_label": "Time (ns)"},
                     "sasa_l": {"title": "SASA", "y_label": "SASA(A²)", "x_label": "Time (ns)"}}


def plot_property_simple(property_dict, sample_data, sample_label, output_dir, out_name, residue):
    """
    Plots a given number of thematically related local metrics of a single sample.

    :param property_dict: Dictionary with the related metrics, each with their title and x and y labels
    :param sample_data: Data for all metrics for a single sample
    :param sample_label: Label for the sample
    :param output_dir: Directory in which the results will be saved
    :param out_name: Name for the thematically related plot (e.g., Interactions, Energetics...)
    :param residue: Residue of interest (usually the mutated one)
    :return: Saves the desired plot in the output folder
    """
    plots = len(property_dict)
    fig = plt.figure(figsize=(plots * 2.5 + 2.5, 3.5))
    axes = []
    sample_replicas = [x for x in sample_data]
    i = 0
    for item in property_dict:
        n_replicas = len(sample_replicas)
        colors = colors_replicas(n_replicas)
        ax = fig.add_subplot(1, plots, i + 1)
        j = 0
        for replica in sample_replicas:
            plot_data = sample_data[replica][item]
            ax.plot(plot_data[0], plot_data[1], label=replica, c=next(colors), alpha=1 - j * (0.5 / n_replicas))
            j += 1
        ax.set_title(property_dict[item]["title"], fontsize=18)
        ax.set_xlabel(property_dict[item]["x_label"])
        ax.set_ylabel(sample_label + " " + property_dict[item]["y_label"])
        ax.legend(loc="best")
        ax.margins(0.1)
        axes.append(ax)
        i += 1
    fig.tight_layout()
    fig.savefig(output_dir + "%s_%s_resid%s.png" % (out_name, sample_label, str(residue)))
    plt.close(fig)


def plot_property_compared(property_dict, all_data_dict, labels_dict, samples, output_dir, out_name, residue):
    """
    Plots a given number of thematically related global metrics of a sample (or samples) compared to the reference one.

    :param property_dict: Dictionary with the related metrics, each with their title and x and y labels
    :param all_data_dict: Data for all metrics for all samples
    :param labels_dict: Dictionary with the label for each sample
    :param samples: Samples to be compared (usually reference and mutant)
    :param output_dir: Directory in which the results will be saved
    :param out_name: Name for the thematically related plot (e.g., Interactions, Energetics...)
    :param residue: Residue of interest (usually the mutated one)
    :return: Saves the desired plot in the output folder
    """
    sample_compared_label = "_vs_".join(samples)
    plot_out_dir = output_dir + sample_compared_label + "/Multianalysis/Plots/"
    plots = len(property_dict)
    rows = len(samples)
    fig = plt.figure(figsize=(plots * 2.5 + 2.5, rows * 2.5 + 1))
    i = 0
    axes = []
    for sample in samples:
        sample_data = all_data_dict[sample]
        sample_label = labels_dict[sample]
        sample_replicas = [x for x in sample_data]
        for item in property_dict:
            n_replicas = len(sample_replicas)
            colors = colors_replicas(n_replicas)
            ax = fig.add_subplot(rows, plots, i + 1)
            j = 0
            for replica in sample_replicas:
                plot_data = sample_data[replica][item]
                ax.plot(plot_data[0], plot_data[1], label=replica, c=next(colors), alpha=1 - j * (0.5 / n_replicas))
                j += 1
            ax.set_title(property_dict[item]["title"], fontsize=18)
            ax.set_xlabel(property_dict[item]["x_label"])
            ax.set_ylabel(sample_label + " " + property_dict[item]["y_label"])
            ax.legend(loc="best")
            ax.margins(0.1)
            axes.append(ax)
            i += 1
    for i in range(plots):
        xmin = 9999999999999
        xmax = -9999999999999
        ymin = 9999999999999
        ymax = -9999999999999
        for j in range(rows):
            xlim = axes[i + j * plots].get_xlim()
            xmin = min(xlim[0], xmin)
            xmax = max(xlim[1], xmax)
            ylim = axes[i + j * plots].get_ylim()
            ymin = min(ylim[0], ymin)
            ymax = max(ylim[1], ymax)
        for j in range(rows):
            axes[i + j * plots].set_xlim(xmin, xmax)
            axes[i + j * plots].set_ylim(ymin, ymax)
            axes[i + j * plots].legend(loc="best")
    fig.tight_layout()
    fig.savefig(plot_out_dir + "%s_%s_resid%s.png" % (out_name, sample_compared_label, str(residue)))
    plt.close(fig)


def plot_local_simple(sample_data, sample_label, output_dir, residue):
    """
    Plots all local metrics, grouped thematically, for a single sample

    :param sample_data: Data for all metrics for a single sample
    :param sample_label: Label for the sample
    :param output_dir: Directory in which the results will be saved
    :param residue: Residue of interest (usually the mutated one)
    :return: Saves all metrics' plots in the output folder for the single sample
    """
    plot_property_simple(interactions_dict, sample_data, sample_label, output_dir, "Local_interactions", residue)


def plot_local_compared(all_data_dict, labels_dict, samples, output_dir, residue):
    """
    Plots all local metrics, grouped thematically, of a sample (or samples) compared to the reference one.

    :param all_data_dict: Data for all metrics for all samples
    :param labels_dict: Dictionary with the label for each sample
    :param samples: Samples to be compared (usually reference and mutant)
    :param output_dir: Directory in which the results will be saved
    :param residue: Residue of interest (usually the mutated one)
    :return: Saves all metrics' plots in the output folder for the compared samples
    """
    plot_property_compared(interactions_dict, all_data_dict, labels_dict, samples, output_dir, "Local_interactions",
                           residue)
