import glob
from natsort import natsorted
import os
from data import AllData
import pickle
from draw_graph import draw_graph
from make_html import HtmlWriter
from make_csv import make_csv

class FileType:
    def __init__(self, d, t, r, file_name):
        self.d = d
        self.t = t
        self.r = r
        self.file_name = file_name

    def get_info(self):
        return [self.d, self.t, self.r, self.file_name]

def classify_file_by_condition(all_file_type, conditions, d_arr, t_arr, r_arr):
    ret_arr = []
    use_arr1, use_arr2, use_arr3 = [[], [], []]
    index1, index2, index3 = [0, 0, 0]

    if conditions == "d":
        use_arr1 = d_arr
        index1 = 0
    elif conditions == "t":
        use_arr1 = t_arr
        index1 = 1
    elif conditions == "r":
        use_arr1 = r_arr
        index1 = 2
    elif conditions == "dt":
        use_arr1 = d_arr
        use_arr2 = t_arr
        index1 = 0
        index2 = 1
    elif conditions == "dr":
        use_arr1 = d_arr
        use_arr2 = r_arr
        index1 = 0
        index2 = 2
    elif conditions == "tr":
        use_arr1 = t_arr
        use_arr2 = r_arr
        index1 = 1
        index2 = 2
    else:
        use_arr1 = d_arr
        use_arr2 = t_arr
        use_arr3 = r_arr
        index1 = 0
        index2 = 1
        index3 = 2

    if len(conditions) == 1:
        for i in use_arr1:
            temp_arr = []
            for file_type in all_file_type:
                info = file_type.get_info()
                if i == info[index1]:
                    temp_arr.append(info[3])
            ret_arr.append(temp_arr)
    elif len(conditions) == 2:
        for i in use_arr1:
            for j in use_arr2:
                temp_arr = []
                for file_type in all_file_type:
                    info = file_type.get_info()
                    if i == info[index1] and j == info[index2]:
                        temp_arr.append(info[3])
                ret_arr.append(temp_arr)
    else:
        for i in use_arr1:
            for j in use_arr2:
                for k in use_arr3:
                    temp_arr = []
                    for file_type in all_file_type:
                        info = file_type.get_info()
                        if i == info[index1] and j == info[index2] and k == info[index3]:
                            temp_arr.append(info[3])
                    ret_arr.append(temp_arr)

    return ret_arr

# duplication:d / movement_type:t / dictance:r
def classify_file(all_file_list):
    classify_dict = {}
    header = ["d", "t", "r", "dt", "dr", "tr", "dtr"]
    all_file_type = []
    d_arr = []
    t_arr = []
    r_arr = []

    for f in all_file_list:
        file_name = f.split('/')[2].split('.')[0].split('_')
        all_file_type.append(FileType(file_name[1], file_name[2], file_name[0], f))
        if not file_name[0] in r_arr:
            r_arr.append(file_name[0])
        if not file_name[1] in d_arr:
            d_arr.append(file_name[1])
        if not file_name[2] in t_arr:
            t_arr.append(file_name[2])

    for i in header:
        classify_dict[i] = classify_file_by_condition(all_file_type, i, d_arr, t_arr, r_arr)

    classify_dict["all"] = all_file_list

    return classify_dict

def search_input_file():
    file_list = glob.glob("./dat/*.dat")
    file_list = natsorted(file_list) # 自然順ソート
    return file_list

def create_pickle_file(all_file_list):
    data_dict = {}
    if os.path.isfile('pickle_file_data.txt'):
        with open('pickle_file_data.txt', 'rb') as f:
            data_dict = pickle.load(f)
    else:
        count = 1
        for file_name in all_file_list:
            print("{}/{} - {} file reading...".format(count, len(all_file_list), file_name.split("/")[2]))
            data_dict[file_name] = AllData(file_name)
            count += 1
        with open('pickle_file_data.txt', 'wb') as f:
            pickle.dump(data_dict, f)

    classify_dict = {}
    if os.path.isfile('pickle_file_list.txt'):
        with open('pickle_file_list.txt', 'rb') as f:
            classify_dict = pickle.load(f)
    else:
        count = 1
        for file_name in all_file_list:
            classify_dict = classify_file(all_file_list)
        with open('pickle_file_list.txt', 'wb') as f:
            pickle.dump(classify_dict, f)

    return [data_dict, classify_dict]

def main():
    data_dict = {}
    all_file_list = search_input_file()

    data_dict, classify_dict = create_pickle_file(all_file_list)
    draw_graph(data_dict, classify_dict)
    # HtmlWriter(data_dict, classify_dict["all"])
    # make_csv(data_dict, classify_dict["all"])

if __name__ == '__main__':
    main()
