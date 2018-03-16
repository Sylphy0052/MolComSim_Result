import matplotlib.pyplot as plt
from matplotlib import ticker
import os
import numpy as np

COLOR_LIST = ['r', 'g', 'b', 'm', 'c', 'y', 'k']
STYLE_LIST = ['-', '--', '-.', ':', '-', '--', '-.', ':']

def check_directory(path):
    if not os.path.isdir(path):
        os.makedirs(path)

def calc_max_range(data_dict, file_list):
    result_range = []
    for file_name in file_list:
        data = data_dict[file_name]
        if len(result_range) < len(data.plot_data.plot_range_by_duplication):
            result_range = data.plot_data.plot_range_by_duplication

    return result_range


def draw_two_datas_graph(X, Y1, Y2, X_label, Y_labels, location):
    fig, ax1 = plt.subplots()
    ln1 = ax1.plot(X, Y1, color=COLOR_LIST[0], label=Y_labels[0])
    ax2 = ax1.twinx()
    ln2 = ax2.plot(X, Y2, color=COLOR_LIST[1], label=Y_labels[1], linestyle=STYLE_LIST[0])

    ax1.set_xlabel(X_label)
    ax1.set_ylabel(Y_labels[2])
    ax2.set_ylabel(Y_labels[3])

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    plt.grid(True)
    plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    plt.gca().ticklabel_format(style="sci",  axis="x",scilimits=(0,0))

    plt.savefig(fig_name, dpi=90, bbox_inches="tight", pad_inches=0.0)
    plt.close('all')
    all_data.set_prob_cumprob_fig_name(fig_name)

def draw_bar_graph(X, Y, X_label, Y_label, fig_name):
    for x, y in zip(X, Y):
        y = round(y, 1)
        plt.text(x, y, y, ha='center', va='bottom')

    plt.ylabel(Y_label)
    plt.bar(X, Y, color=COLOR_LIST, tick_label=X_label, width=0.5)

    plt.savefig(fig_name, dpi=90, bbox_inches="tight", pad_inches=0.0)
    plt.close('all')

def draw_line_graph(X, Y, X_label, Y_label, location, fig_name):
    for i in range(len(Y)):
        plt.plot(X, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5")

    plt.xlabel(X_label)
    plt.ylabel(Y_label)
    plt.ylim(ymin=0)

    plt.grid(True)
    plt.legend(loc=location)

    plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

    plt.savefig(fig_name)
    plt.close('all')

def draw_by_simulation(data_dict, classify_dict):
    for file_name in classify_dict["all"]:
        draw_prob_cumprob_by_simulation(data_dict[file_name])
        draw_mean_by_simulation(data_dict[file_name])

def draw_prob_cumprob_by_simulation(all_data):
    plot_data = all_data.plot_data
    fig_name = "./result/" + all_data.file_name + "_prob_cumprob.png

    X = [0]
    Y1 = [0]
    Y2 = [0]

    Y_labels = ["Probability", "Cumulative Probability", "Probability of RTT", "Cumulative Probability of RTT"]

    X.extend([x1 + 1 for x in plot_data.plot_range])
    Y1.extend(plot_data.prob)
    Y2.extend(plot_data.cum_prob)

    draw_two_datas_graph(X, Y1, Y2, "RTT", "Probability", "Cumulative Probability", "right")

def draw_mean_by_simulation(all_data):
    fig_name = "./result/" + all_data.file_name + "_mean.png"
    X = []
    X_label = []

    output_data = all_data.output_data
    mean = output_data.mean
    analytical_mean = output_data.analytical_model.rtt

    if all_data.is_ptime:
        txrx_mean = output_data.txrx_mean
        X = range(3)
        X_label = ["Mean", "TxRx Mean", "Analytical"]
        Y = [float(mean), float(txrx_mean), float(analytical_mean)]
    else:
        X = range(2)
        X_label = ["Mean", "Analytical"]
        Y = [float(mean), float(analytical_mean)]

    draw_bar_graph(X, Y, X_label, Y_label, fig_name)

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

        draw_line_graph(X, Y, "Distance from Tx to Rx", "Mean of RTT", X_labels, "upper left", fig_name)

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

        draw_line_graph(X, Y, "RTT", "Cumulative Probability of RTT", X_labels, "lower right", fig_name)
