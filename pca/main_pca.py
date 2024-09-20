from matplotlib import pyplot as plt
from matplotlib import colormaps
import matplotlib
from Multianalysis.multianalysis_functions import global_multianalysis, local_multianalysis
from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
global_pca_analysis_functions = ["gyration", "hbonds_inter", "hbonds_intra", "native_contacts", "rmsd",
                                 "rmsdist", "sasa", "ss", "coil", "tm_score"]
local_pca_analysis_functions = ["hbonds_inter_l", "hbonds_intra_l", "native_contacts_l", "sasa_l"]
cmaps = [colormaps[i] for i in ["Blues", "Oranges", "Greens", "Purples", "Greys", "Reds", "BuPu", "OrRd", "YlGn", "PuRd"]]


def main_pca(input_dir, output_dir, samples_replicas_dict, reference_sample, labels_dict, position_dict, time_step,
             final, starting, protein_type, out_pca, parts=20):
    all_data_dict = {}
    ref_label = labels_dict[reference_sample]
    all_data_dict[reference_sample] = global_multianalysis(input_dir, output_dir, reference_sample,
                                                           samples_replicas_dict[reference_sample],
                                                           time_step, final, starting, protein_type, ref_label)
    for sample in samples_replicas_dict:
        if sample != reference_sample:
            sample_label = labels_dict[sample]
            all_data_dict[sample] = global_multianalysis(input_dir, output_dir, sample, samples_replicas_dict[sample],
                                                         time_step, final, starting, protein_type, sample_label)
            if position_dict[sample]:
                reference_replicas = samples_replicas_dict[reference_sample]
                replicas = samples_replicas_dict[sample]
                residue = position_dict[sample]
                local_multianalysis(input_dir, output_dir, reference_sample, sample, reference_replicas, replicas,
                                    labels_dict, residue)
    canonic_x = []
    y_dict = {}
    for sample in samples_replicas_dict:
        if sample not in y_dict:
            y_dict[sample] = {}
        for replica in all_data_dict[sample]:
            if replica not in y_dict[sample]:
                y_dict[sample][replica] = {}
            for function in global_pca_analysis_functions:
                if not canonic_x:
                    canonic_x = [float(round(i, 2)) for i in all_data_dict[sample][replica][function][0]]
                new_y_data = ["NA" for i in canonic_x]
                for i in range(len(all_data_dict[sample][replica][function][0])):
                    time = float(round(all_data_dict[sample][replica][function][0][i], 2))
                    data_point = all_data_dict[sample][replica][function][1][i]
                    if time in canonic_x:
                        list_index = canonic_x.index(time)
                        new_y_data[list_index] = data_point
                y_dict[sample][replica][function] = new_y_data
    y_per_replica = {}
    y = []
    for sample in y_dict:
        for replica in y_dict[sample]:
            name = sample + "_" + replica
            y_replica = []
            for i in range(len(canonic_x)):
                y_time = []
                for function in global_pca_analysis_functions:
                    y_time.append(y_dict[sample][replica][function][i])
                y.append(y_time)
                y_replica.append(y_time)
            y_per_replica[name] = y_replica
    header = global_pca_analysis_functions
    pca_pipe = make_pipeline(StandardScaler(), PCA())
    pca_pipe.fit(y)
    pca_model = pca_pipe.named_steps['pca']
    # Getting PCA weights and explained variance
    file_out = open(out_pca + "All_PCA_global/PC_stats.txt", "w")
    file_out.write("PC  %variance " + " ".join(header) + "\n")
    comps = pca_model.components_
    variances = pca_model.explained_variance_ratio_
    for i in range(len(variances)):
        pc = "PC%i " % (i + 1)
        variance = str(round(variances[i] * 100, 2)) + "% "
        comp = [str(j) for j in comps[i]]
        file_out.write(pc + variance + " ".join(comp) + "\n")
    file_out.close()
    # Getting the PC values
    x = canonic_x * len(y_per_replica)
    pc_data = pca_pipe.transform(y)
    file_out = open(out_pca + "All_PCA_global/PC_values.txt", "w")
    file_out.write("IID " + " ".join("PC%i" % (j + 1) for j in range(len(pc_data[0]))) + "\n")
    for i in range(len(pc_data)):
        file_out.write("%s %s\n" % (x[i], " ".join([str(j) for j in pc_data[i]])))
    file_out.close()
    for name in y_per_replica:
        pc_data = pca_pipe.transform(y_per_replica[name])
        plt.plot(pc_data[:, 0], pc_data[:, 1], "o", label=name)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend(loc="best")
    plt.savefig("%sAll_PCA_global/PC_compare.jpg" % out_pca)
    plt.close()
    for i in range(parts):
        j = 0
        # General plot
        for name in y_per_replica:
            start = int(round(i * len(y_per_replica[name]) / parts, 0))
            end = int(round((i + 1) * len(y_per_replica[name]) / parts, 0) + 1)
            new_data = y_per_replica[name][start:end]
            scale = range(start, start + len(new_data))
            pc_data = pca_pipe.transform(new_data)
            plt.scatter(pc_data[:, 0], pc_data[:, 1], cmap=cmaps[j % 10], c=scale)
            plt.scatter([], [], color=cmaps[j % 10](0.5), label=name)
            j += 1
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        norm = matplotlib.colors.Normalize(vmin=min(scale), vmax=max(scale))
        bar = plt.colorbar(matplotlib.cm.ScalarMappable(cmap=colormaps["Greys"], norm=norm))
        bar.ax.set_ylabel("Time (ns)")
        xlim = plt.xlim()
        ylim = plt.ylim()
        plt.legend(loc="best")
        plt.savefig("%sAll_PCA_global/part%i_PC_compare.jpg" % (out_pca, i))
        plt.close()
        # Replica plots:
        j = 0
        for name in y_per_replica:
            start = int(round(i * len(y_per_replica[name]) / parts, 0))
            end = int(round((i + 1) * len(y_per_replica[name]) / parts, 0) + 1)
            new_data = y_per_replica[name][start:end]
            scale = range(start, start + len(new_data))
            pc_data = pca_pipe.transform(new_data)
            plt.scatter(pc_data[:, 0], pc_data[:, 1], cmap=cmaps[j % 10], c=scale)
            plt.scatter([], [], color=cmaps[j % 10](0.5), label=name)
            j += 1
            plt.xlabel("PC1")
            plt.ylabel("PC2")
            norm = matplotlib.colors.Normalize(vmin=min(scale), vmax=max(scale))
            bar = plt.colorbar(matplotlib.cm.ScalarMappable(cmap=colormaps["Greys"], norm=norm))
            bar.ax.set_ylabel("Time (ns)")
            name_split = name.split("_")
            mut = name_split[0]
            rep = "_".join(name_split[1:])
            plt.xlim(xlim)
            plt.ylim(ylim)
            plt.legend(loc="best")
            # plt.savefig("%s%s/PCA_global/%s/part%i_PC_compare.jpg" % (out_pca, mut, rep, i))
            plt.savefig("%sAll_PCA_global/%s_%s/part%i_PC_compare.jpg" % (out_pca, mut, rep, i))
            plt.close()
        # Compared plot
        for sample in samples_replicas_dict:
            j = 0
            if sample != reference_sample:
                namelist = ["%s_%s" % (reference_sample, x) for x in samples_replicas_dict[reference_sample]]
                namelist += ["%s_%s" % (sample, x) for x in samples_replicas_dict[sample]]
                for name in y_per_replica:
                    start = int(round(i * len(y_per_replica[name]) / parts, 0))
                    end = int(round((i + 1) * len(y_per_replica[name]) / parts, 0) + 1)
                    new_data = y_per_replica[name][start:end]
                    scale = range(start, start + len(new_data))
                    pc_data = pca_pipe.transform(new_data)
                    plt.scatter(pc_data[:, 0], pc_data[:, 1], cmap=cmaps[j % 10], c=scale)
                    plt.scatter([], [], color=cmaps[j % 10](0.5), label=name)
                    j += 1
                plt.xlabel("PC1")
                plt.ylabel("PC2")
                norm = matplotlib.colors.Normalize(vmin=min(scale), vmax=max(scale))
                bar = plt.colorbar(matplotlib.cm.ScalarMappable(cmap=colormaps["Greys"], norm=norm))
                bar.ax.set_ylabel("Time (ns)")
                plt.legend(loc="best")
                plt.savefig("%s%s_vs_%s/PCA/Global/All/part%i_PC_compare.jpg" % (out_pca, reference_sample, sample, i))
                plt.close()
                j = 0
                for name in y_per_replica:
                    print(name, sample, namelist)
                    if name in namelist:
                        start = int(round(i * len(y_per_replica[name]) / parts, 0))
                        end = int(round((i + 1) * len(y_per_replica[name]) / parts, 0) + 1)
                        new_data = y_per_replica[name][start:end]
                        scale = range(start, start + len(new_data))
                        pc_data = pca_pipe.transform(new_data)
                        plt.scatter(pc_data[:, 0], pc_data[:, 1], cmap=cmaps[j % 10], c=scale)
                        plt.scatter([], [], color=cmaps[j % 10](0.5), label=name)
                        j += 1
                        plt.xlabel("PC1")
                        plt.ylabel("PC2")
                        norm = matplotlib.colors.Normalize(vmin=min(scale), vmax=max(scale))
                        bar = plt.colorbar(matplotlib.cm.ScalarMappable(cmap=colormaps["Greys"], norm=norm))
                        bar.ax.set_ylabel("Time (ns)")
                        plt.legend(loc="best")
                        plt.savefig("%s%s_vs_%s/PCA/Global/%s/part%i_PC_compare.jpg" % (out_pca, reference_sample, sample,
                                                                                        name, i))
                        plt.close()
    for sample in samples_replicas_dict:
        if sample != reference_sample and position_dict[sample]:
            print(sample, reference_sample)
            y_dict = {}
            y_dict[reference_sample] = {}
            y_dict[sample] = {}
            reference_replicas = samples_replicas_dict[reference_sample]
            replicas = samples_replicas_dict[sample]
            residue = position_dict[sample]
            all_data_dict_local = local_multianalysis(input_dir, output_dir, reference_sample, sample,
                                                      reference_replicas, replicas, labels_dict, residue)
            for item in y_dict:
                for replica in all_data_dict_local[item]:
                    if replica not in y_dict[item]:
                        y_dict[item][replica] = {}
                    for function in global_pca_analysis_functions:
                        if not canonic_x:
                            canonic_x = [float(round(i, 2)) for i in all_data_dict[item][replica][function][0]]
                        new_y_data = ["NA" for i in canonic_x]
                        for i in range(len(all_data_dict[item][replica][function][0])):
                            time = float(round(all_data_dict[item][replica][function][0][i], 2))
                            data_point = all_data_dict[item][replica][function][1][i]
                            if time in canonic_x:
                                list_index = canonic_x.index(time)
                                new_y_data[list_index] = data_point
                        y_dict[item][replica][function] = new_y_data
                for replica in all_data_dict_local[item]:
                    for function in local_pca_analysis_functions:
                        if not canonic_x:
                            canonic_x = [float(round(i, 2)) for i in all_data_dict_local[item][replica][function][0]]
                        new_y_data = ["NA" for i in canonic_x]
                        for i in range(len(all_data_dict_local[item][replica][function][0])):
                            time = float(round(all_data_dict_local[item][replica][function][0][i], 2))
                            data_point = all_data_dict_local[item][replica][function][1][i]
                            if time in canonic_x:
                                list_index = canonic_x.index(time)
                                new_y_data[list_index] = data_point
                        y_dict[item][replica][function] = new_y_data
                        print(item, replica, function, y_dict[item][replica][function])
            y = []
            y_per_replica = {}
            for y_sample in y_dict:
                for replica in y_dict[y_sample]:
                    name = y_sample + "_" + replica
                    y_replica = []
                    for i in range(len(canonic_x)):
                        y_time = []
                        for function in global_pca_analysis_functions + local_pca_analysis_functions:
                            y_time.append(y_dict[y_sample][replica][function][i])
                        y.append(y_time)
                        y_replica.append(y_time)
                    y_per_replica[name] = y_replica
            # Obtaining PCA local
            header = global_pca_analysis_functions + local_pca_analysis_functions
            pca_pipe = make_pipeline(StandardScaler(), PCA())
            pca_pipe.fit(y)
            pca_model = pca_pipe.named_steps['pca']
            # Getting PCA weights and explained variance
            file_out = open("%s%s_vs_%s/PCA/Local/All/PC_stats.txt" % (out_pca, reference_sample, sample), "w")
            file_out.write("PC  %variance " + " ".join(header) + "\n")
            comps = pca_model.components_
            variances = pca_model.explained_variance_ratio_
            for i in range(len(variances)):
                pc = "PC%i " % (i + 1)
                variance = str(round(variances[i] * 100, 2)) + "% "
                comp = [str(j) for j in comps[i]]
                file_out.write(pc + variance + " ".join(comp) + "\n")
            file_out.close()
            # Getting the PC values
            x = canonic_x * len(y_per_replica)
            pc_data = pca_pipe.transform(y)
            file_out = open("%s%s_vs_%s/PCA/Local/All/PC_values.txt" % (out_pca, reference_sample, sample), "w")
            file_out.write("IID " + " ".join("PC%i" % (j + 1) for j in range(len(pc_data[0]))) + "\n")
            for i in range(len(pc_data)):
                file_out.write("%s %s\n" % (x[i], " ".join([str(j) for j in pc_data[i]])))
            file_out.close()
            for i in range(parts):
                # General plot
                j = 0
                for name in y_per_replica:
                    start = int(round(i * len(y_per_replica[name]) / parts, 0))
                    end = int(round((i + 1) * len(y_per_replica[name]) / parts, 0) + 1)
                    new_data = y_per_replica[name][start:end]
                    scale = range(start, start + len(new_data))
                    pc_data = pca_pipe.transform(new_data)
                    plt.scatter(pc_data[:, 0], pc_data[:, 1], cmap=cmaps[j % 10], c=scale)
                    plt.scatter([], [], color=cmaps[j % 10](0.5), label=name)
                    j += 1
                plt.xlabel("PC1")
                plt.ylabel("PC2")
                norm = matplotlib.colors.Normalize(vmin=min(scale), vmax=max(scale))
                bar = plt.colorbar(matplotlib.cm.ScalarMappable(cmap=colormaps["Greys"], norm=norm))
                bar.ax.set_ylabel("Time (ns)")
                xlim = plt.xlim()
                ylim = plt.ylim()
                plt.legend(loc="best")
                plt.savefig("%s%s_vs_%s/PCA/Local/All/part%i_PC_compare.jpg" % (out_pca, reference_sample, sample, i))
                plt.close()
                # Replica plots:
                j = 0
                for name in y_per_replica:
                    start = int(round(i * len(y_per_replica[name]) / parts, 0))
                    end = int(round((i + 1) * len(y_per_replica[name]) / parts, 0) + 1)
                    new_data = y_per_replica[name][start:end]
                    scale = range(start, start + len(new_data))
                    print(name, "window", new_data)
                    pc_data = pca_pipe.transform(new_data)
                    plt.scatter(pc_data[:, 0], pc_data[:, 1], cmap=cmaps[j % 10], c=scale)
                    plt.scatter([], [], color=cmaps[j % 10](0.5), label=name)
                    j += 1
                    plt.xlabel("PC1")
                    plt.ylabel("PC2")
                    norm = matplotlib.colors.Normalize(vmin=min(scale), vmax=max(scale))
                    bar = plt.colorbar(matplotlib.cm.ScalarMappable(cmap=colormaps["Greys"], norm=norm))
                    bar.ax.set_ylabel("Time (ns)")
                    name_split = name.split("_")
                    mut = name_split[0]
                    rep = "_".join(name_split[1:])
                    plt.xlim(xlim)
                    plt.ylim(ylim)
                    plt.legend(loc="best")
                    plt.savefig("%s%s_vs_%s/PCA/Local/%s_%s/part%i_PC_compare.jpg" % (out_pca, reference_sample, sample,
                                                                                      mut, rep, i))
                    plt.close()
