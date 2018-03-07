import glob
from natsort import natsorted
import os
import pickle
from data import SimData
from draw_graph import drawGraph
from make_html import HtmlWriter

def createDataDict(file_list):
    data_dict = {}
    if os.path.isfile('pickle_data.txt'):
        with open('pickle_data.txt', 'rb') as f:
            data_dict = pickle.load(f)
    else:
        count = 1
        for file_name in file_list:
            print("{}/{} - {} file reading...".format(count, len(file_list), file_name.split("/")[2]))
            data_dict[file_name] = SimData(file_name)
            count += 1
        with open('pickle_data.txt', 'wb') as f:
            pickle.dump(data_dict, f)

    return data_dict

def main():
    all_file_list = natsorted(glob.glob("./dat/*.dat"))
    data_dict = createDataDict(all_file_list)
    drawGraph(data_dict, all_file_list)
    HtmlWriter(data_dict, all_file_list)

if __name__ == '__main__':
    main()
