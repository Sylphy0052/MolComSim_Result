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

def search_by_diffusioncoefficient(data_dict, file_list):
    X = []
    X_labels = []
    for file_name in file_list:
        datas = data_dict[file_name].get_median_by_diffusioncoefficient_each_distance()
        label = "d={}".format(datas[1])
        if not datas[0] in X:
            X.append(datas[0])
        if not label in X_labels:
            X_labels.append(label)

    return [X, X_labels]

def search_by_distance_each_rto(data_dict, file_list):
    X = []
    X_labels = []
    for file_name in file_list:
        datas = data_dict[file_name].get_median_by_distance_each_rto()
        label = "RTO{}".format(datas[1])
        if label == "RTO0":
            label = "NoRetransmit"
        elif label == "RTO1":
            label = "RTO1 = 2 * Median"
        elif label == "RTO2":
            label = "RTO2 = Median + 1/3 * stdev"
        elif label == "RTO3":
            label = "RTO3 = Median + 1/2 * stdev"
        elif label == "RTO4":
            label = "RTO4 = Median + stdev"

        if not datas[0] in X:
            X.append(datas[0])
        if not label in X_labels:
            X_labels.append(label)

    return [X, X_labels]

def search_regression_by_distance_each_rto(data_dict, file_list):
    X = []
    X_labels = []
    for file_name in file_list:
        datas = data_dict[file_name].get_median_by_distance_each_rto()
        label = "RTO{}".format(datas[1])
        if label == "RTO0":
            label = "NoRetransmit"

        if not datas[0] in X:
            X.append(datas[0])
        if not label in X_labels:
            X_labels.append(label)

    return [X, X_labels]

def search_by_rto_each_distance(data_dict, file_list):
    X = []
    X_labels = []
    for file_name in file_list:
        datas = data_dict[file_name].get_by_rto_each_distance()
        label = "r={}".format(datas[0])
        x = "RTO{}".format(datas[1])

        if not label in X_labels:
            X_labels.append(label)
        if not x in X:
            X.append(x)

    return [X, X_labels]

def search_each_distance(data_dict, file_list):
    X = []
    X_labels = []
    for file_name in file_list:
        datas = data_dict[file_name].get_median_by_diffusioncoefficient_each_distance()
        label = "r={}".format(datas[1])

        # if not datas[0] in X:
        #     X.append(datas[0])
        if not label in X_labels:
            X_labels.append(label)

    # return [X, X_labels]
    return X_labels

def search_by_distance_each_decomposing(data_dict, file_list):
    X = []
    X_labels = []
    for file_name in file_list:
        data = data_dict[file_name]
        if data.config.config_dict["decomposing"] == 1:
            label = "Decomposing"
        else:
            label = "Non Decomposing"

        if not data.r in X:
            X.append(data.r)
        if not label in X_labels:
            X_labels.append(label)

    return [X, X_labels]

def draw_regression_graph(XY, X_labels, xName, yName, fig_name):
    for i in range(len(X_labels)):
        X = []
        Y = []
        for k, v in sorted(XY[i].items()):
            X.append(k)
            Y.append(v)

        y = np.poly1d(np.polyfit(X, Y, len(X) - 1))(X)
        y1 = np.poly1d(np.polyfit(X, Y, len(X)-1))

        f_c = ": $"
        for j in range(len(y1.c)):
            num_str = ""

            if len(y1.c) - j == 2:
                num_str = "{0:+g}x".format(y1.c[j], 3)
            elif len(y1.c) - j == 1:
                num_str = "{0:+g}".format(y1.c[j], 3)
            else:
                num_str = "{0:+g}x^{1}".format(y1.c[j], len(y1.c) - 1 - j)

            f_c += num_str

        f_c += '$'

        X_labels[i] += f_c
        X_labels[i] = repr(X_labels[i])
        plt.plot(X, Y, COLOR_LIST[i] + "o")
        plt.plot(X, y, color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

    plt.xlabel(xName)
    plt.ylabel(yName)
    plt.grid(True)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)

    plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

    plt.savefig(fig_name, bbox_inches='tight')
    plt.close('all')

def draw_by_simulation(data_dict, classify_dict):
    for file_name in classify_dict["all"]:
        draw_prob_cumprob_by_simulation(data_dict[file_name])
        draw_mean_by_simulation(data_dict[file_name])

def draw_prob_cumprob_by_simulation(all_data):
    plot_data = all_data.plot_data
    fig_name = "./result/" + all_data.file_name + "_prob_cumprob.png"

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


    plt.savefig(fig_name, dpi=90, bbox_inches="tight", pad_inches=0.0)
    plt.close('all')
    all_data.set_prob_cumprob_fig_name(fig_name)

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

    for x, y in zip(X, Y):
        y = round(y, 1)
        plt.text(x, y, y, ha='center', va='bottom')

    plt.ylabel("Mean of RTT")
    plt.bar(X, Y, color=COLOR_LIST, tick_label=X_label, width=0.5)

    plt.savefig(fig_name, dpi=90, bbox_inches="tight", pad_inches=0.0)
    plt.close('all')

    all_data.set_mean_fig_name(fig_name)

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
        plt.legend(loc='upper right')

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
        plt.legend(loc='upper right')

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

def draw_median_by_diffusioncoefficient(data_dict, classify_dict):
    dir_path = "./compare_median_by_diffusioncoefficient/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]
    for file_list in file_list_by_conditions:
        sorted(file_list, reverse=True)
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "ARQ{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])
        X, X_labels = search_by_diffusioncoefficient(data_dict, file_list)
        Y = [[] for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_median_by_diffusioncoefficient_each_distance()
            for i in range(len(X_labels)):
                if "d={}".format(datas[1]) in X_labels[i]:
                    Y[i].append(datas[2])
                    break

        for i in range(len(Y)):
            plt.plot(X, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5")

        plt.xlabel('Diffusion Coefficient')
        plt.ylabel('Median of RTT')
        plt.ylim(ymin=0)
        # plt.xticks(X, [str(x) for x in X])
        plt.grid(True)
        plt.legend(loc='upper right')

        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def draw_regression_median_by_diffusioncoefficient_each_distance(data_dict, classify_dict):
    dir_path = "./regression_median_by_diffusioncoefficient_each_distance/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "ARQ{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        X, X_labels = search_by_diffusioncoefficient(data_dict, file_list)
        Y = [[] for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_median_by_diffusioncoefficient_each_distance()
            for i in range(len(X_labels)):
                if "d={}".format(datas[1]) in X_labels[i]:
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

def draw_median_by_distance_each_rto(data_dict, classify_dict):
    dir_path = "./compare_median_by_distance_each_rto/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        X, X_labels = search_by_distance_each_rto(data_dict, file_list)
        Y = [[] for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_median_by_distance_each_rto()
            for i in range(len(X_labels)):
                if "RTO{}".format(datas[1]) in X_labels[i]:
                    Y[i].append(datas[2])
                    break
                elif str(datas[1]) == "0":
                    Y[i].append(datas[2])
                    break

        for i in range(len(Y)):
            plt.plot(X, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('Distance from Tx to Rx')
        plt.ylabel('Median of RTT')
        plt.grid(True)
        plt.legend(loc='upper left')

        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def draw_regression_median_by_distance_each_rto(data_dict, classify_dict):
    dir_path = "./compare_regression_median_by_distance_each_rto/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        _, X_labels = search_regression_by_distance_each_rto(data_dict, file_list)
        XY = [{} for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_median_by_distance_each_rto()
            for i in range(len(X_labels)):
                if "RTO{}".format(datas[1]) in X_labels[i]:
                    XY[i][datas[0]] = datas[2]
                    # Y[i].append(datas[2])
                    break
                elif str(datas[1]) == "0":
                    XY[i][datas[0]] = datas[2]
                    # Y[i].append(datas[2])
                    break

        draw_regression_graph(XY, X_labels, 'Distance from Tx to Rx', 'Median of RTT', fig_name)

def draw_median_by_rto_each_distance(data_dict, classify_dict):
    dir_path = "./compare_median_by_rto_each_distance/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        X, X_labels = search_by_rto_each_distance(data_dict, file_list)
        x = range(len(X))
        Y = [[] for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_median_by_rto_each_distance()
            for i in range(len(X_labels)):
                if "r={}".format(datas[0]) in X_labels[i]:
                    Y[i].append(datas[2])
                    break

        for i in range(len(Y)):
            plt.plot(x, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('RTO Type')
        plt.ylabel('Median of RTT')
        plt.grid(True)
        plt.legend(loc='upper left')
        plt.xticks(x, X)
        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def draw_jitter_by_rto_each_distance(data_dict, classify_dict):
    dir_path = "./compare_jitter_by_rto_each_distance/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        X, X_labels = search_by_rto_each_distance(data_dict, file_list)
        x = range(len(X))
        Y = [[] for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_jitter_by_rto_each_distance()
            for i in range(len(X_labels)):
                if "r={}".format(datas[0]) in X_labels[i]:
                    Y[i].append(datas[2])
                    break

        for i in range(len(Y)):
            plt.plot(x, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('RTO Type')
        plt.ylabel('Jitter of RTT')
        plt.grid(True)
        plt.legend(loc='upper left')
        plt.xticks(x, X)
        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

###
def draw_mean_by_distance_each_rto(data_dict, classify_dict):
    dir_path = "./compare_mean_by_distance_each_rto/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        X, X_labels = search_by_distance_each_rto(data_dict, file_list)
        Y = [[] for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_mean_by_distance_each_rto()
            for i in range(len(X_labels)):
                if "RTO{}".format(datas[1]) in X_labels[i]:
                    Y[i].append(datas[2])
                    break
                elif str(datas[1]) == "0":
                    Y[i].append(datas[2])
                    break

        for i in range(len(Y)):
            plt.plot(X, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('Distance from Tx to Rx')
        plt.ylabel('Mean of RTT')
        plt.grid(True)
        plt.legend(loc='upper left')

        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def draw_regression_mean_by_distance_each_rto(data_dict, classify_dict):
    dir_path = "./compare_regression_mean_by_distance_each_rto/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        _, X_labels = search_regression_by_distance_each_rto(data_dict, file_list)
        # Y = [[] for i in range(len(X_labels))]
        XY = [{} for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_mean_by_distance_each_rto()
            for i in range(len(X_labels)):
                if "RTO{}".format(datas[1]) in X_labels[i]:
                    XY[i][datas[0]] = datas[2]
                    # Y[i].append(datas[2])
                    break
                elif str(datas[1]) == "0":
                    XY[i][datas[0]] = datas[2]
                    # Y[i].append(datas[2])
                    break

        draw_regression_graph(XY, X_labels, 'Distance from Tx to Rx', 'Mean of RTT', fig_name)

def draw_mean_by_rto_each_distance(data_dict, classify_dict):
    dir_path = "./compare_mean_by_rto_each_distance/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        X, X_labels = search_by_rto_each_distance(data_dict, file_list)
        x = range(len(X))
        Y = [[] for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_mean_by_rto_each_distance()
            for i in range(len(X_labels)):
                if "r={}".format(datas[0]) in X_labels[i]:
                    Y[i].append(datas[2])
                    break

        for i in range(len(Y)):
            plt.plot(x, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('RTO Type')
        plt.ylabel('Mean of RTT')
        plt.grid(True)
        plt.legend(loc='upper left')
        plt.xticks(x, X)
        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def draw_regression_jitter_by_distance_each_rto(data_dict, classify_dict):
    dir_path = "./compare_regression_jitter_by_distance_each_rto/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        _, X_labels = search_regression_by_distance_each_rto(data_dict, file_list)
        XY = [{} for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_jitter_by_distance_each_rto()
            for i in range(len(X_labels)):
                if "RTO{}".format(datas[1]) in X_labels[i]:
                    XY[i][datas[0]] = datas[2]
                    # Y[i].append(datas[2])
                    break
                elif str(datas[1]) == "0":
                    XY[i][datas[0]] = datas[2]
                    # Y[i].append(datas[2])
                    break

        draw_regression_graph(XY, X_labels, 'Distance from Tx to Rx', 'Jitter of RTT', fig_name)

def draw_jitter_by_distance_each_rto(data_dict, classify_dict):
    dir_path = "./compare_jitter_by_distance_each_rto/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        X, X_labels = search_by_distance_each_rto(data_dict, file_list)
        Y = [[] for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_jitter_by_distance_each_rto()
            for i in range(len(X_labels)):
                if "RTO{}".format(datas[1]) in X_labels[i]:
                    Y[i].append(datas[2])
                    break
                elif str(datas[1]) == "0":
                    Y[i].append(datas[2])
                    break

        for i in range(len(Y)):
            plt.plot(X, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('Distance from Tx to Rx')
        plt.ylabel('Jitter of RTT')
        plt.grid(True)
        plt.legend(loc='upper left')

        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def draw_retransmission_wait_time_by_rto_each_distance(data_dict, classify_dict):
    dir_path = "./compare_retransmission_wait_time_by_rto_each_distance/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        X, X_labels = search_by_rto_each_distance(data_dict, file_list)
        x = range(len(X))
        Y = [[] for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_retransmission_wait_time_by_rto_each_distance()
            for i in range(len(X_labels)):
                if "r={}".format(datas[0]) in X_labels[i]:
                    Y[i].append(datas[2])
                    break

        for i in range(len(Y)):
            plt.plot(x, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('RTO Type')
        plt.ylabel('Retransmission Wait Time')
        plt.grid(True)
        plt.legend(loc='upper left')
        plt.xticks(x, X)
        # plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        # plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def draw_median_by_retransmission_wait_time_each_distance(data_dict, classify_dict):
    dir_path = "./compare_median_by_retransmission_wait_time_each_distance/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        # X, X_labels = search_by_diffusioncoefficient_each_distance(data_dict, file_list)
        X_labels = search_each_distance(data_dict, file_list)
        # X = [[] for i in range(len(X_labels))]
        # Y = [[] for i in range(len(X_labels))]

        XY = [{} for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_median_by_retransmission_wait_time_each_distance()
            for i in range(len(X_labels)):
                if "r={}".format(datas[1]) == X_labels[i]:
                    # X[i].append(datas[0])
                    # Y[i].append(datas[2])
                    XY[i][datas[0]] = datas[2]

        for i in range(len(X_labels)):
            x = []
            y = []
            for k, v in sorted(XY[i].items()):
                x.append(k)
                y.append(v)
            plt.plot(x, y, color=COLOR_LIST[i], label=X_labels[i], marker='o', markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('Diffusion Coefficient')
        plt.ylabel('Median of RTT')
        plt.grid(True)
        plt.legend(loc='upper left')

        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def draw_mean_by_retransmission_wait_time_each_distance(data_dict, classify_dict):
    dir_path = "./compare_mean_by_retransmission_wait_time_each_distance/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        X_labels = search_each_distance(data_dict, file_list)

        XY = [{} for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_mean_by_retransmission_wait_time_each_distance()
            for i in range(len(X_labels)):
                if "r={}".format(datas[1]) == X_labels[i]:
                    XY[i][datas[0]] = datas[2]

        for i in range(len(X_labels)):
            x = []
            y = []
            for k, v in sorted(XY[i].items()):
                x.append(k)
                y.append(v)
            plt.plot(x, y, color=COLOR_LIST[i], label=X_labels[i], marker='o', markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('Diffusion Coefficient')
        plt.ylabel('Mean of RTT')
        plt.grid(True)
        plt.legend(loc='upper left')

        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def draw_jitter_by_retransmission_wait_time_each_distance(data_dict, classify_dict):
    dir_path = "./compare_jitter_by_retransmission_wait_time_each_distance/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        X_labels = search_each_distance(data_dict, file_list)

        XY = [{} for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_jitter_by_retransmission_wait_time_each_distance()
            for i in range(len(X_labels)):
                if "r={}".format(datas[1]) == X_labels[i]:
                    XY[i][datas[0]] = datas[2]

        for i in range(len(X_labels)):
            x = []
            y = []
            for k, v in sorted(XY[i].items()):
                x.append(k)
                y.append(v)
            plt.plot(x, y, color=COLOR_LIST[i], label=X_labels[i], marker='o', markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('Diffusion Coefficient')
        plt.ylabel('Jitter of RTT')
        plt.grid(True)
        plt.legend(loc='upper left')

        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def draw_regression_median_by_retransmission_wait_time_each_distance(data_dict, classify_dict):
    dir_path = "./compare_regression_median_by_retransmission_wait_time_each_distance/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        X_labels = search_each_distance(data_dict, file_list)

        XY = [{} for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_median_by_retransmission_wait_time_each_distance()
            for i in range(len(X_labels)):
                if "r={}".format(datas[1]) == X_labels[i]:
                    XY[i][datas[0]] = datas[2]

        draw_regression_graph(XY, X_labels, 'Diffusion Coefficient', 'Median of RTT', fig_name)

def draw_regression_mean_by_retransmission_wait_time_each_distance(data_dict, classify_dict):
    dir_path = "./compare_regression_mean_by_retransmission_wait_time_each_distance/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        X_labels = search_each_distance(data_dict, file_list)

        XY = [{} for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_mean_by_retransmission_wait_time_each_distance()
            for i in range(len(X_labels)):
                if "r={}".format(datas[1]) == X_labels[i]:
                    XY[i][datas[0]] = datas[2]

        draw_regression_graph(XY, X_labels, 'Diffusion Coefficient', 'Mean of RTT', fig_name)

def draw_regression_jitter_by_retransmission_wait_time_each_distance(data_dict, classify_dict):
    dir_path = "./compare_regression_jitter_by_retransmission_wait_time_each_distance/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        X_labels = search_each_distance(data_dict, file_list)

        XY = [{} for i in range(len(X_labels))]

        for file_name in file_list:
            datas = data_dict[file_name].get_jitter_by_retransmission_wait_time_each_distance()
            for i in range(len(X_labels)):
                if "r={}".format(datas[1]) == X_labels[i]:
                    XY[i][datas[0]] = datas[2]

        draw_regression_graph(XY, X_labels, 'Diffusion Coefficient', 'Jitter of RTT', fig_name)

# decomposing
def classifyDecomposingDict(data_dict, classify_dict):
    decomposing_file = []
    non_decomposing_file = []

    for file_name in classify_dict['all']:
        data = data_dict[file_name]
        if data.config.config_dict["decomposing"] == 1:
            decomposing_file.append(file_name)
        else:
            non_decomposing_file.append(file_name)

    return [decomposing_file, non_decomposing_file]

def draw_mean_by_distance_each_decomposing(data_dict, classify_dict):
    dir_path = "./compare_mean_by_distance_each_decomposing/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        X, X_labels = search_by_distance_each_decomposing(data_dict, file_list)
        Y = [[] for i in range(len(X_labels))]

        for file_name in file_list:
            data = data_dict[file_name]
            for i in range(len(X_labels)):
                if data.config.config_dict["decomposing"] == 1 and "Decomposing" == X_labels[i]:
                    Y[i].append(data.mean)
                    break
                elif data.config.config_dict["decomposing"] == 0 and "Non Decomposing" == X_labels[i]:
                    Y[i].append(data.mean)
                    break

        for i in range(len(Y)):
            plt.plot(X, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('Distance from Tx to Rx')
        plt.ylabel('Mean of RTT')
        plt.grid(True)
        plt.legend(loc='upper left')

        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def draw_median_by_distance_each_decomposing(data_dict, classify_dict):
    dir_path = "./compare_median_by_distance_each_decomposing/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        X, X_labels = search_by_distance_each_decomposing(data_dict, file_list)
        Y = [[] for i in range(len(X_labels))]

        for file_name in file_list:
            data = data_dict[file_name]
            for i in range(len(X_labels)):
                if data.config.config_dict["decomposing"] == 1 and "Decomposing" == X_labels[i]:
                    Y[i].append(data.med)
                    break
                elif data.config.config_dict["decomposing"] == 0 and "Non Decomposing" == X_labels[i]:
                    Y[i].append(data.med)
                    break

        for i in range(len(Y)):
            plt.plot(X, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('Distance from Tx to Rx')
        plt.ylabel('Median of RTT')
        plt.grid(True)
        plt.legend(loc='upper left')

        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')

def draw_jitter_by_distance_each_decomposing(data_dict, classify_dict):
    dir_path = "./compare_jitter_by_distance_each_decomposing/"
    check_directory(dir_path)
    file_list_by_conditions = classify_dict["dt"]

    for file_list in file_list_by_conditions:
        info = data_dict[file_list[0]].get_info()
        fig_name = dir_path + "{}-{}_{}-{}.png".format(info[2], info[3], info[4], info[5])

        X, X_labels = search_by_distance_each_decomposing(data_dict, file_list)
        Y = [[] for i in range(len(X_labels))]

        for file_name in file_list:
            data = data_dict[file_name]
            for i in range(len(X_labels)):
                if data.config.config_dict["decomposing"] == 1 and "Decomposing" == X_labels[i]:
                    Y[i].append(data.std)
                    break
                elif data.config.config_dict["decomposing"] == 0 and "Non Decomposing" == X_labels[i]:
                    Y[i].append(data.std)
                    break

        for i in range(len(Y)):
            plt.plot(X, Y[i], color=COLOR_LIST[i], label=X_labels[i], markersize="5", linestyle=STYLE_LIST[i])

        plt.xlabel('Distance from Tx to Rx')
        plt.ylabel('Jitter of RTT')
        plt.grid(True)
        plt.legend(loc='upper left')

        plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
        plt.gca().ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        plt.savefig(fig_name)
        plt.close('all')


def draw_regression_mean_by_distance_each_decomposing(data_dict, classify_dict):
    pass

def draw_regression_median_by_distance_each_decomposing(data_dict, classify_dict):
    pass

def draw_regression_jitter_by_distance_each_decomposing(data_dict, classify_dict):
    pass

# Adjust
def draw_adjust_graph(data_dict, classify_dict):
    check_directory("./adjust_graph")
    for file_name in classify_dict["all"]:
        steps = []
        info_adjust_num = []
        ack_adjust_num = []

        adjust_file_name = data_dict[file_name].adjust_file_name
        fig_name = "./adjust_graph/" + adjust_file_name.split(".txt")[0].split("/")[2] + ".png"
        with open(adjust_file_name, 'r') as f:
            for line in f:
                datas = line.split(',')
                for data in datas:
                    info, ack = data.split('/')
                    info_adjust_num.append(info)
                    ack_adjust_num.append(ack)
        for i in range(len(info_adjust_num) - 1):
            steps.append(i * 100)
        steps.append(data_dict[file_name].input_data.steps[-1])

        plt.plot(steps, info_adjust_num, color=COLOR_LIST[0], label="INFO Num", markersize="5")
        plt.plot(steps, ack_adjust_num, color=COLOR_LIST[1], label="ACK Num", markersize="5")

        plt.xlabel('Steps')
        plt.ylabel('Molecule Number')

        plt.grid(True)
        plt.legend(loc='upper left')

        plt.savefig(fig_name)
        plt.close('all')

def draw_graph(data_dict, classify_dict):
    # Each simulation graph
    draw_by_simulation(data_dict, classify_dict)

    # draw_mean_by_distance(data_dict, classify_dict)
    # draw_cumprob_each_duplication(data_dict, classify_dict)
    # draw_median_by_distance_each_duplication(data_dict, classify_dict)
    # draw_mean_by_distance_each_duplication(data_dict, classify_dict)
    # draw_mean_by_duplication_each_distance(data_dict, classify_dict)
    # draw_txrx_mean_by_distance_each_duplication(data_dict, classify_dict)
    # draw_jitter_by_duplication_each_distance(data_dict, classify_dict)
    # draw_regression_median_by_duplication_each_distance(data_dict, classify_dict)
    # draw_regression_jitter_by_duplication_each_distance(data_dict, classify_dict)

    # Diffusion Coefficient
    # draw_median_by_diffusioncoefficient(data_dict, classify_dict)
    # draw_regression_median_by_diffusioncoefficient_each_distance(data_dict, classify_dict)

    # RTO
    # draw_retransmission_wait_time_by_rto_each_distance(data_dict, classify_dict)
    # draw_regression_median_by_distance_each_rto(data_dict, classify_dict)
    # draw_median_by_rto_each_distance(data_dict, classify_dict)
    # draw_jitter_by_rto_each_distance(data_dict, classify_dict)
    # draw_regression_mean_by_distance_each_rto(data_dict, classify_dict)
    # draw_mean_by_rto_each_distance(data_dict, classify_dict)
    # draw_regression_jitter_by_distance_each_rto(data_dict, classify_dict)
    # draw_median_by_retransmission_wait_time_each_distance(data_dict, classify_dict)
    # draw_mean_by_retransmission_wait_time_each_distance(data_dict, classify_dict)
    # draw_jitter_by_retransmission_wait_time_each_distance(data_dict, classify_dict)
    # draw_regression_median_by_retransmission_wait_time_each_distance(data_dict, classify_dict)
    # draw_regression_mean_by_retransmission_wait_time_each_distance(data_dict, classify_dict)
    # draw_regression_jitter_by_retransmission_wait_time_each_distance(data_dict, classify_dict)

    # decomposing
    # draw_mean_by_distance_each_decomposing(data_dict, classify_dict)
    # draw_median_by_distance_each_decomposing(data_dict, classify_dict)
    # draw_jitter_by_distance_each_decomposing(data_dict, classify_dict)
    #
    # draw_regression_mean_by_distance_each_decomposing(data_dict, classify_dict)
    # draw_regression_median_by_distance_each_decomposing(data_dict, classify_dict)
    # draw_regression_jitter_by_distance_each_decomposing(data_dict, classify_dict)

    # Adjust
    draw_adjust_graph(data_dict, classify_dict)
