import pandas as pd
from data import MovementType

def is_same_conditions(search_data, target_data):
    # 距離と重複度
    if search_data.info_arq != target_data.info_arq:
        return False
    if search_data.ack_arq != target_data.ack_arq:
        return False
    if search_data.r != target_data.r:
        return False
    return True


def search_same_conditions(data_dict, search_file_name, file_list, already_array):
    ret_array = [search_file_name]
    search_data = data_dict[search_file_name]
    tmp_arr = []
    for file_name in file_list:
        if file_name in already_array:
            continue
        target_data = data_dict[file_name]
        if is_same_conditions(search_data, target_data):
            tmp_arr.append(file_name)
            ret_array.append(file_name)

    return tmp_arr, ret_array

def to_csv_string(data_dict, file_list):
    for file_name in file_list:
        if data_dict[file_name].info_type == MovementType["PASSIVE"]:
            pass_data = data_dict[file_name].get_csv_data()
        elif data_dict[file_name].info_type == MovementType["ACTIVE"]:
            act_data = data_dict[file_name].get_csv_data()
    header = "{}_{}-{}".format(pass_data[0], pass_data[1], pass_data[2])
    return "{}, {}, {}, {}, {}, {}, {}\n".format(header, round(pass_data[3], 3), round(act_data[3], 3), round(pass_data[4], 3), round(act_data[4], 3), round(pass_data[5], 3), round(act_data[5], 3))

def write_header():
    return ",PASSIVE_MEAN,ACTIVE_MEAN,PASSIVE_MEDIAN,ACTIVE_MEDIAN,PASSIVE_TxRxMEAN,ACTIVE_TxRxMEAN\n"

def make_csv(data_dict, file_list):
    already_array = []
    file_list_by_conditions = []
    for file_name in file_list:
        if file_name in already_array:
            continue
        already_array.append(file_name)
        tmp_arr, same_arr = search_same_conditions(data_dict, file_name, file_list, already_array)
        already_array.extend(tmp_arr)
        if len(same_arr) != 0:
            file_list_by_conditions.append(same_arr)

    ret_str = ""
    for file_arr in file_list_by_conditions:
        ret_str += to_csv_string(data_dict, file_arr)

    with open("csv_data.txt", "w") as f:
        f.write(write_header())
        f.write(ret_str)

    data = pd.read_csv("csv_data.txt")
    data.to_csv("data.csv", index=False)
