import glob
from natsort import natsorted
import os
import sys
from statistics import mean
from matplotlib import ticker
import numpy as np
import pickle

HEADER = """
<!DOCTYPE html>
<html lang=\"ja\">
<head>
<meta charset=\"UTF-8\">
<title>Result</title>
"""

STYLE = """
<style>
.tabbox input { display: none; }

.tab {
    display: inline-block;
    border: thin solid black;
    /* border-bottom: 1px solid #5ab4bd; */
    /* box-shadow: 0 0 10px rgba(0, 0, 0, 0.2); */
    width: 30em;
    /* padding: 0.50em 0.75em; */
    padding: 0.25em 0.20em;
    /* color: #565656; */
    color: #000000;
    background-color: #d9d9d9;
    font-weight: bold;
    font-size: 10px;
}

.tab:hover {
    opacity: 0.75;
    cursor: pointer;
}

input:checked + .tab {
    color: #ffffff;
    background-color: #5ab4bd;
    position: relative;
    z-index: 10;
}

.tabcontent {
    display: none;
    /* border: 1px solid black; */
    margin-top: -1px;
    padding: 1em;
    position: relative;
    z-index: 0;
    /* background-color: #ffffcc; */
}

table { border-collapse: collapse; }
td {
    border: solid 1px;
    padding: 0.5em;
}
thead.scrollHead,tbody.scrollBody { display:block; }
tbody.scrollBody {
    overflow-y:scroll;
    height:750px;
}

td,th {
    table-layout:fixed;
}

"""

TAB_STYLE = """
#tabcheck{}:checked ~ #tabcontent{}{{ display: block; }}
"""

BODY = """
<input type=\"radio\" name=\"tabset\" id=\"tabcheck{}\">
<label for=\"tabcheck{}\" class=\"tab\">{}</label>
"""

BODY_DETAIL = """
<div class=\"tabcontent\" id=\"tabcontent{}\">
{}
</div>
"""

TAB1 = """
<table style=\"font-size:12px;position: absolute; left: 10px; top: 10px;\">
{}
</table>
"""

TAB1_1 = """
<tr>
<td>{}</td>
<td>{}</td>
</tr>
"""

TAB2 = """
<table align=\"left\" style=\"font-size:12px;position: absolute; left: 360px; top: 10px;\">
{}
</table>
"""

TAB2_1 = """
<tr>
<td>{}</td>
<td>{}</td>
</tr>
"""

TAB2_HEADER = ["Count", "Min", "Max", "Range", "VAR", "STD", "MED", "MEAN", "TxRxAVG"]

TAB3 = """
<table align=\"left\" style=\"font-size:12px;position: absolute; left: 360px; top: 360px;\">
{}
</table>
"""

TAB3_1 = """
<tr>
<td>{}</td>
<td>{}</td>
</tr>
"""

TAB3_HEADER = ["Dinf", "r", "L", "R", "RTT"]

TAB4 = """
<table align=\"left\" style=\"font-size:12px;position: absolute; left: 700px; top: 10px;\">
<thead class=\"scrollHead\">
<tr>
<th>Min</th>
<th>Max</th>
<th>Num</th>
<th>Prob</th>
<th>CumProb</th>
</tr>
</thead>
<tbody class=\"scrollBody\">
{}
</tbody>
</table>
"""

TAB4_1 = """
<tr>
<td>{}</td>
<td>{}</td>
<td>{}</td>
<td>{}</td>
<td>{}</td>
</tr>
"""

TAB5 = """
<table align=\"left\" style=\"font-size:12px;position: absolute; left: 360px; top: 560px;\">
{}
</table>
"""

TAB5_1 = """
<tr>
<td>{}</td>
<td>{}</td>
</tr>
"""

TAB5_HEADER = ["coll-a/a", "coll-a/i", "coll-a/n", "coll-i/i", "coll-i/n"]

FIG = """
<img src=\"{}\" align=\"right\" style=\"position: absolute; right: 50px; top: {}px;\">
"""

LINE_COLOR = ['r', 'g', 'b', 'm', 'c', 'y', 'k']

LINE_STYLE = ['-', '--', '-.', ':', '-', '--', '-.', ':']

ARQ_ARRAY = ["ARQ1-1_", "ARQ10-10_", "ARQ20-20_", "ARQ30-30_", "ARQ50-50_", "ARQ100-100_"]
ARQ_ARRAY_COLLISION = ["ARQ10-10_", "ARQ20-20_", "ARQ30-30_", "ARQ50-50_", "ARQ100-100_"]
DISTANCE_ARRAY = ["TxRx30_", "TxRx50_", "TxRx70_", "TxRx90_"]
NOISE_ARRAY = ["0NoiseNoCollisions.", "0Noise."]

class HtmlWriter:
    def __init__(self, data_dict, file_list):
        self.data_dict = data_dict
        self.file_list = file_list
        self.file_name = "index.html"

        if os.path.exists(self.file_name):
            os.remove(self.file_name)

        self.write_html()

    def write_html(self):
        self.write_header()
        self.write_body()

    def write_header(self):
        header = HEADER
        header += self.write_style()
        header += """
</style>
</head>
"""
        self.write_text(header)

    def write_style(self):
        style = STYLE
        index = 0
        for file_name in self.file_list:
            index += 1
            style += TAB_STYLE.format(str(index), str(index))

        return style


    def write_body(self):
        body = "<body><div class=\"tabbox\">"
        index = 0
        for config_file_name in self.file_list:
            file_name = self.data_dict[config_file_name].input_data.input_file_name
            index += 1
            tab_name = file_name.split(".txt")[0].split("batch_")[1]
            body += BODY.format(str(index), str(index), tab_name)

        # 表や図を入れる
        index = 0
        for file_name in self.file_list:
            index += 1
            body += BODY_DETAIL.format(str(index), self.write_detail_body(file_name))

        body += """
</div>
</body>
</html>
"""
        self.write_text(body)

    def write_detail_body(self, file_name):
        all_data = self.data_dict[file_name]
        is_coll = all_data.is_coll
        text = ""
        text += self.create_tab1(all_data.config)
        text += self.create_tab2(all_data.output_data)
        text += self.create_tab3(all_data.output_data.analytical_model)
        # if is_coll:
        #     text += self.create_tab5(all_data.input_data)
        text += self.create_tab4(all_data.plot_data)
        text += self.create_fig_text(all_data.mean_fig_name, 10)
        text += self.create_fig_text(all_data.prob_cumprob_fig_name, 430)

        return text

    def check_config(self, config_dict, key):
        ret_val = ""

        if key in ["transmitter", "receiver"]:
            v = config_dict[key].toString()
            ret_val = TAB1_1.format(key, v)

        elif key in "intermediateNode":
            v = config_dict[key].toString()
            ret_val = TAB1_1.format(key, v)

        elif key in "moleculeParams":
            for val in config_dict[key]:
                v = val.toString()
                ret_val += TAB1_1.format(key, v)

        elif key in "microtubuleParams":
            for val in config_dict[key]:
                v = val.toString()
                ret_val += TAB1_1.format(key, v)

        else:
            ret_val = TAB1_1.format(key, config_dict[key])

        return ret_val

    def create_tab1(self, config):
        tab1 = ""
        for k in config.keys():
            if k == "outputFile":
                continue
            tab1 += self.check_config(config, k)
        return TAB1.format(tab1)

    def create_tab2(self, output_data):
        output_data_arr = output_data.toArray()
        is_ptime = output_data.is_ptime
        tab2 = ""
        for i in range(len(TAB2_HEADER)):
            if i == 8 and not is_ptime:
                continue
            tab2 += TAB2_1.format(TAB2_HEADER[i], round(output_data_arr[i], 1))

        tab2 += TAB2_1.format("Analytical Mean", round(output_data.analytical_model.rtt, 1))

        return TAB2.format(tab2)

    def create_tab3(self, analytical_model):
        tab3 = ""

        model_arr = analytical_model.toArray()
        tab = ""
        for i in range(len(TAB3_HEADER)):
            tab += TAB3_1.format(TAB3_HEADER[i], round(model_arr[i], 2))
        tab3 += TAB3.format(tab)

        return tab3

    def create_tab4(self, plot_data):
        tab4 = ""
        for i in range(len(plot_data.plot_range)):
            tab4 += TAB4_1.format(plot_data.plot_range[i][0], plot_data.plot_range[i][1], len(plot_data.plot_data[i]), round(plot_data.prob[i], 1), round(plot_data.cum_prob[i], 1))

        return TAB4.format(tab4)

    def create_fig_text(self, file_name, height):
        return FIG.format(file_name, height)

    def write_text(self, text):
        with open(self.file_name, 'a') as f:
            f.write(text)

    def create_tab5(self, input_data):
        tab5 = ""
        Y = []
        Y.append(mean(input_data.coll_aa))
        Y.append(mean(input_data.coll_ai))
        Y.append(mean(input_data.coll_an))
        Y.append(mean(input_data.coll_ii))
        Y.append(mean(input_data.coll_in))

        for i in range(len(Y)):
            tab5 += TAB5_1.format(TAB5_HEADER[i], round(Y[i], 1))

        return TAB5.format(tab5)
