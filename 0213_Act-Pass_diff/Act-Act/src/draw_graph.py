import matplotlib.pyplot as plt
from matplotlib import ticker
import os
import numpy as np

COLOR_LIST = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
STYLE_LIST = ['-', '--', '-.', ':', '-', '--', '-.', ':']

def draw_prob_cumprob_by_simulation(all_data):
    plot_data = all_data.plot_data
    X = [0]
    Y1 = [0]
    Y2 = [0]
    X.extend([x[1] + 1 for x in plot_data.plot_range]) # Step count
    Y1.extend(plot_data.prob) # Probability
    Y2.extend(plot_data.cum_prob) # Cumulative Probability

    fig, ax1 = plt.subplots()
    ln1 = ax1.plot(X, Y1, color=COLOR_LIST[0], label="Probability")
    ax2 = ax1.twinx()
    ln2 = ax2.plot(X, Y2, color=COLOR_LIST[1], label="Cumulative Probability", linestyle=STYLE_LIST[0])

    ax1.set_xlabel('RTT')
    ax1.set_ylabel('Probability of RTT')
    ax2.set_ylabel('Cumulative Probability of RTT')

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2, loc='right')

    plt.grid(True)

    plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    plt.gca().ticklabel_format(style="sci",  axis="x",scilimits=(0,0))

    fig_name = "./result/" + all_data.file_name + "_prob_cumprob.png"
    plt.savefig(fig_name, dpi=90, bbox_inches="tight", pad_inches=0.0)
    plt.close('all')
    all_data.set_prob_cumprob_fig_name(fig_name)

def draw_mean_by_simulation(all_data):
    output_data = all_data.output_data
    mean = output_data.mean
    analytical_mean = output_data.analytical_model.rtt
    X = []
    X_label = []

    if all_data.is_ptime:
        txrx_mean = output_data.txrx_mean
        X = range(3)
        X_label = ["Mean", "TxRx Mean", "Analytical"]
        Y = [float(mean), float(txrx_mean), float(analytical_mean)]
    else:
        X = range(2)
        X_label = ["Mean", "Analytical"]
        Y = [float(mean), float(analytical_mean)]

    for x, y in zip(X, Y):
        y = round(y, 1)
        plt.text(x, y, y, ha='center', va='bottom')

    plt.ylabel("Mean of RTT")
    plt.bar(X, Y, color=COLOR_LIST, tick_label=X_label, width=0.5)

    fig_name = "./result/" + all_data.file_name + "_mean.png"
    plt.savefig(fig_name, dpi=90, bbox_inches="tight", pad_inches=0.0)
    plt.close('all')

    all_data.set_mean_fig_name(fig_name)

def draw_by_simulation(data_dict, classify_dict):
    for file_name in classify_dict["all"]:
        draw_prob_cumprob_by_simulation(data_dict[file_name])
        draw_mean_by_simulation(data_dict[file_name])

def check_directory(path):
    if not os.path.isdir(path):
        os.makedirs(path)

def draw_mean_by_distance(data_dict, classify_dict):
    dir_path = "./compare_mean_by_distance/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]
    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        is_ptime = info[0]
        fig_name = dir_path + "ARQ{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])
        X = []
        X_labels = ["Mean", "Analytical"]
        if is_ptime:
            X_labels.append("TxRx Mean")
        Y = [[] for i in range(len(X_labels))]

        finish_arr = []
        for file_name in file_list:
            datas = data_dict[file_name].get_mean_by_distance()
            X.append(datas[0])
            Y[0].append(datas[1])
            Y[1].append(datas[2])
            if is_ptime:
                Y[2].append(datas[3])

        for i in range(len(Y)):
            plt.plot(X, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5")

        plt.xlabel('Distance from Tx to Rx')
        plt.ylabel('Mean of RTT')
        plt.ylim(ymin=0)
        # plt.xticks(X, [str(x) for x in X])
        plt.grid(True)
        plt.legend(loc='upper left')

        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def calc_max_range(data_dict, file_list):
    result_range = []
    for file_name in file_list:
        data = data_dict[file_name]
        if len(result_range) < len(data.plot_data.plot_range_by_duplication):
            result_range = data.plot_data.plot_range_by_duplication

    return result_range

def draw_cumprob_each_duplication(data_dict, classify_dict):
    dir_path = "./compare_cumprob_each_duplication/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["tr"]

    X_arr = calc_max_range(data_dict, classify_dict["all"])
    X = [0]
    for x in X_arr:
        X.append(x[1])

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        is_ptime = info[0]
        fig_name = dir_path + "TxRx{}_{}-{}.png".format(info[1], info[4], info[5])
        Y = []
        X_labels = []

        for file_name in file_list:
            datas = data_dict[file_name].get_cumprob_each_duplication()
            X_labels.append("SW-ARQ{}_{}".format(datas[0], datas[1]))
            Y.append(datas[2])

        for i in range(len(Y)):
            Y[i].insert(0, 0)
            while len(X) != len(Y[i]):
                Y[i].append(100)
            plt.plot(X, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('RTT')
        plt.ylabel('Cumulative Probability of RTT')
        plt.ylim([0, 100])
        plt.grid(True)
        plt.legend(loc='lower right')

        plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="x",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def search_by_duplication(data_dict, file_list):
    X = []
    X_labels = []
    Y_analytical = []
    for file_name in file_list:
        datas = data_dict[file_name].get_median_by_distance_each_duplication()
        if not datas[0] in X:
            X.append(datas[0])
            Y_analytical.append(datas[4])
        label = "SW-ARQ{}-{}".format(datas[1], datas[2])
        if not label in X_labels:
            X_labels.append(label)
    X_labels.append("Analytical")
    return [X, X_labels, Y_analytical]

def draw_median_by_distance_each_duplication(data_dict, classify_dict):
    dir_path = "./compare_median_by_distance_each_duplication/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["t"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}.png".format(info[4], info[5])

        X, X_labels, Y_analytical = search_by_duplication(data_dict, file_list)
        Y = [[] for i in range(len(X_labels))]
        Y[-1] = Y_analytical

        for file_name in file_list:
            datas = data_dict[file_name].get_median_by_distance_each_duplication()
            for i in range(len(X_labels)):
                if "{}-{}".format(datas[1], datas[2]) in X_labels[i]:
                    Y[i].append(datas[3])
                    break

        for i in range(len(Y)):
            plt.plot(X, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('Distance from Tx to Rx')
        plt.ylabel('Median of RTT')
        # plt.xticks(X, [str(x) for x in X])
        plt.grid(True)
        plt.legend(loc='upper left')

        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def draw_mean_by_distance_each_duplication(data_dict, classify_dict):
    dir_path = "./compare_mean_by_distance_each_duplication/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["t"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}.png".format(info[4], info[5])

        X, X_labels, Y_analytical = search_by_duplication(data_dict, file_list)
        Y = [[] for i in range(len(X_labels))]
        Y[-1] = Y_analytical

        for file_name in file_list:
            datas = data_dict[file_name].get_mean_by_distance_each_duplication()
            for i in range(len(X_labels)):
                if "{}-{}".format(datas[1], datas[2]) in X_labels[i]:
                    Y[i].append(datas[3])
                    break

        for i in range(len(Y)):
            plt.plot(X, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('Distance from Tx to Rx')
        plt.ylabel('Mean of RTT')
        # plt.xticks(X, [str(x) for x in X])
        plt.grid(True)
        plt.legend(loc='upper left')

        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def search_by_distance(data_dict, file_list):
    X = []
    X_labels = []
    for file_name in file_list:
        datas = data_dict[file_name].get_mean_by_duplication_each_distance()
        label = "d={}".format(datas[0])
        if not datas[1] in X:
            X.append(datas[1])
        if not label in X_labels:
            X_labels.append(label)

    return [X, X_labels]

def draw_mean_by_duplication_each_distance(data_dict, classify_dict):
    dir_path = "./compare_mean_by_duplication_each_distance/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["t"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}.png".format(info[4], info[5])

        X, X_labels = search_by_distance(data_dict, file_list)
        Y = [[] for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_mean_by_duplication_each_distance()
            for i in range(len(X_labels)):
                if "d={}".format(datas[0]) in X_labels[i]:
                    Y[i].append(datas[2])
                    break

        for i in range(len(Y)):
            plt.plot(X, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('Duplication')
        plt.ylabel('Mean of RTT')
        # plt.xticks(X, [str(x) for x in X])
        plt.grid(True)
        plt.legend(loc='upper left')

        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def draw_txrx_mean_by_distance_each_duplication(data_dict, classify_dict):
    dir_path = "./compare_txrx_mean_by_distance_each_duplication/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["t"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        if not info[0]:
            break
        fig_name = dir_path + "{}-{}.png".format(info[4], info[5])

        X, X_labels, Y_analytical = search_by_duplication(data_dict, file_list)
        Y = [[] for i in range(len(X_labels))]
        Y[-1] = Y_analytical

        for file_name in file_list:
            datas = data_dict[file_name].get_txrx_mean_by_distance_each_duplication()
            for i in range(len(X_labels)):
                if "{}-{}".format(datas[1], datas[2]) in X_labels[i]:
                    Y[i].append(datas[3])
                    break

        for i in range(len(Y)):
            plt.plot(X, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('Distance from Tx to Rx')
        plt.ylabel('Mean of RTT')
        # plt.xticks(X, [str(x) for x in X])
        plt.grid(True)
        plt.legend(loc='upper left')

        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def draw_jitter_by_duplication_each_distance(data_dict, classify_dict):
    dir_path = "./compare_jitter_by_duplication_each_distance/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["t"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}.png".format(info[4], info[5])

        X, X_labels = search_by_distance(data_dict, file_list)
        Y = [[] for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_std_by_duplication_each_distance()
            for i in range(len(X_labels)):
                if "d={}".format(datas[0]) in X_labels[i]:
                    Y[i].append(datas[2])
                    break

        for i in range(len(Y)):
            plt.plot(X, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('Duplication')
        plt.ylabel('Jitter of RTT')
        # plt.xticks(X, [str(x) for x in X])
        plt.grid(True)
        plt.legend(loc='upper left')

        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def draw_regression_median_by_duplication_each_distance(data_dict, classify_dict):
    dir_path = "./regression_median_by_duplication_each_distance/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["t"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}.png".format(info[4], info[5])

        X, X_labels = search_by_distance(data_dict, file_list)
        Y = [[] for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_median_by_duplication_each_distance()
            for i in range(len(X_labels)):
                if "d={}".format(datas[0]) in X_labels[i]:
                    Y[i].append(datas[2])
                    break

        for i in range(len(Y)):
            y = np.poly1d(np.polyfit(X, Y[i], len(X) - 1))(X)
            y1 = np.poly1d(np.polyfit(X, Y[i], len(X)-1))

            f_c = ": $"
            for j in range(len(y1.c)):
                num_str = ""

                if len(y1.c) - j == 2:
                    num_str = "{0:+g}x".format(round(y1.c[j], 1))
                elif len(y1.c) - j == 1:
                    num_str = "{0:+g}".format(round(y1.c[j], 1))
                else:
                    num_str = "{0:+g}x^{1}".format(round(y1.c[j], 1), len(y1.c) - 1 - j)

                f_c += num_str

            f_c += '$'

            X_labels[i] += f_c
            X_labels[i] = repr(X_labels[i])
            plt.plot(X, Y[i], COLOR_LIST[i] + "o")
            plt.plot(X, y, color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('Duplication')
        plt.ylabel('Median of RTT')
        # plt.xticks(X, [str(x) for x in X])
        plt.grid(True)
        plt.legend(loc='upper right', fontsize=8)
        # plt.legend(bbox_to_anchor=(0, 1.2), loc=2, borderaxespad=0)

        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def draw_regression_jitter_by_duplication_each_distance(data_dict, classify_dict):
    dir_path = "./regression_jitter_by_duplication_each_distance/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["t"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}.png".format(info[4], info[5])

        X, X_labels = search_by_distance(data_dict, file_list)
        Y = [[] for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_std_by_duplication_each_distance()
            for i in range(len(X_labels)):
                if "d={}".format(datas[0]) in X_labels[i]:
                    Y[i].append(datas[2])
                    break

        for i in range(len(Y)):
            y = np.poly1d(np.polyfit(X, Y[i], len(X) - 1))(X)
            y1 = np.poly1d(np.polyfit(X, Y[i], len(X)-1))

            f_c = ": $"
            for j in range(len(y1.c)):
                num_str = ""

                if len(y1.c) - j == 2:
                    num_str = "{0:+g}x".format(round(y1.c[j], 1))
                elif len(y1.c) - j == 1:
                    num_str = "{0:+g}".format(round(y1.c[j], 1))
                else:
                    num_str = "{0:+g}x^{1}".format(round(y1.c[j], 1), len(y1.c) - 1 - j)

                f_c += num_str

            f_c += '$'

            X_labels[i] += f_c
            X_labels[i] = repr(X_labels[i])


            plt.plot(X, Y[i], COLOR_LIST[i] + "o")
            plt.plot(X, y, color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('Duplication')
        plt.ylabel('Jitter of RTT')
        # plt.xticks(X, [str(x) for x in X])
        plt.grid(True)
        plt.legend(loc='upper right', fontsize=8)

        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def draw_mean_by_diffusioncoefficient(data_dict, classify_dict):
    dir_path = "./compare_mean_by_diffusioncoefficient/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]
    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "ARQ{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])
        X = []
        X_labels = ["Measured Value"]
        Y = [[] for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_mean_by_diffusioncoefficient()
            X.append(datas[0])
            Y[0].append(datas[1])

        for i in range(len(Y)):
            plt.plot(X, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5")

        plt.xlabel('Diffusion Coefficient')
        plt.ylabel('Mean of RTT')
        plt.ylim(ymin=0)
        # plt.xticks(X, [str(x) for x in X])
        plt.grid(True)
        plt.legend(loc='upper right')

        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')


def draw_graph(data_dict, classify_dict):
    draw_by_simulation(data_dict, classify_dict)
    draw_mean_by_distance(data_dict, classify_dict)
    draw_cumprob_each_duplication(data_dict, classify_dict)
    draw_median_by_distance_each_duplication(data_dict, classify_dict)
    draw_mean_by_distance_each_duplication(data_dict, classify_dict)
    draw_mean_by_duplication_each_distance(data_dict, classify_dict)
    draw_txrx_mean_by_distance_each_duplication(data_dict, classify_dict)
    draw_jitter_by_duplication_each_distance(data_dict, classify_dict)
    draw_regression_median_by_duplication_each_distance(data_dict, classify_dict)
    draw_regression_jitter_by_duplication_each_distance(data_dict, classify_dict)
    draw_mean_by_diffusioncoefficient(data_dict, classify_dict)
