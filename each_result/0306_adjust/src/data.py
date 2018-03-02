import re
from enum import Enum
import linecache
import sys
import os
import numpy as np
from statistics import mean, variance, stdev, median
import math

DNA_DICT = {
    11.31: {"diameter": 0.01, "dc": 63.944},
    9.19: {"diameter": 0.02, "dc": 42.187},
    7.46: {"diameter": 0.02, "dc": 27.833},
    6.06: {"diameter": 0.04, "dc": 18.363},
    4.92: {"diameter": 0.06, "dc": 12.115},
    4.00: {"diameter": 0.09, "dc": 7.993},
    3.25: {"diameter": 0.13, "dc": 5.273},
    2.64: {"diameter": 0.20, "dc": 3.479},
    2.14: {"diameter": 0.30, "dc": 2.295},
    1.18: {"diameter": 0.98, "dc": 0.695},
    0.67: {"diameter": 3.00, "dc": 0.227},
    2.20: {"diameter": 0.28, "dc": 2.416},
    1.42: {"diameter": 0.67, "dc": 1.013},
    10.58: {"diameter": 0.01, "dc": 55.931},
}

class MoleculeType(Enum):
    INFO = 0
    ACK = 1
    NOISE = 2

class MovementType(Enum):
    PASSIVE = 0
    ACTIVE = 1

class AllData:
    def __init__(self, config_file_name):
        self.config = Config(config_file_name)
        self.input_data = InputData(self.config)
        self.output_data = OutputData(self.config, self.input_data)
        self.plot_data = PlotData(self.output_data)
        self.is_ptime = self.input_data.is_ptime
        self.is_coll = self.input_data.is_coll
        self.define_file_name()

    def define_file_name(self):
        self.r = int(self.output_data.analytical_model.r)
        molecule_params = self.config.config_dict["moleculeParams"]
        for mp in molecule_params:
            if mp.type_of_molecule == MoleculeType["INFO"]:
                self.info_arq = mp.num_of_molecules
                self.info_type = mp.type_of_movement
            elif mp.type_of_molecule == MoleculeType["ACK"]:
                self.ack_arq = mp.num_of_molecules
                self.ack_type = mp.type_of_movement
        self.step_length = self.config.config_dict["stepLengthX"]
        self.define_compare_data()

        if "RTO" in self.config.config_file_name:
            self.file_name = "TxRx{}_ARQ{}-{}_{}-{}_RTO{}".format(self.r, self.info_arq, self.ack_arq, self.info_type.name, self.ack_type.name, self.rto_type)
        else:
            self.file_name = "TxRx{}_ARQ{}-{}_{}-{}".format(self.r, self.info_arq, self.ack_arq, self.info_type.name, self.ack_type.name)

        if self.config.config_dict["decomposing"] == 1:
            self.file_name += "_decomposing"

        if self.config.isAdjust:
            self.file_name += "_adjust{}".format(self.config.adjust_num)
            self.adjust_file_name = "./result/adjust_batch_" + self.config.config_dict["outputFile"]

    def define_compare_data(self):
        self.mean = self.output_data.mean
        self.analytical_rtt = self.output_data.analytical_model.rtt
        self.med = self.output_data.med
        self.std = self.output_data.std
        self.diffusion_coefficient = self.config.config_dict["stepLengthX"]
        self.retransmitWaitTime = self.config.config_dict["retransmitWaitTime"]
        self.txrx_mean = self.output_data.txrx_mean
        if "RTO" in self.config.config_file_name:
            self.rto_type = self.config.config_file_name.split("RTO")[1][0]

    def get_compare_info(self):
        return [self.r,
                self.info_arq,
                self.ack_arq,
                self.info_type.name,
                self.ack_type.name,
                "{03.2f}".format(self.diffusion_coefficient).replace('.', ''),
                "{:05d}".format(self.retransmitWaitTime)]

    def set_mean_fig_name(self, fig_name):
        self.mean_fig_name = fig_name

    def get_mean_fig_name(self):
        return self.mean_fig_name

    def set_prob_cumprob_fig_name(self, fig_name):
        self.prob_cumprob_fig_name = fig_name

    def get_prob_cumprob_fig_name(self):
        return self.prob_cumprob_fig_name

    def get_info(self):
        return [self.is_ptime, self.r, self.info_arq, self.ack_arq, self.info_type.name, self.ack_type.name]

    def get_mean_by_distance(self):
        ret_arr = [self.r, self.mean, self.analytical_rtt]
        if self.is_ptime:
            ret_arr.append(self.txrx_mean)
        return ret_arr

    def get_cumprob_each_duplication(self):
        return [self.info_arq, self.ack_arq, self.plot_data.cum_prob_by_duplication]

    def get_median_by_distance_each_duplication(self):
        return [self.r, self.info_arq, self.ack_arq, self.med, self.analytical_rtt]

    def get_mean_by_distance_each_duplication(self):
        return [self.r, self.info_arq, self.ack_arq, self.med, self.analytical_rtt]

    def get_mean_by_duplication_each_distance(self):
        return [self.r, self.info_arq, self.mean]

    def get_txrx_mean_by_distance_each_duplication(self):
        return [self.r, self.info_arq, self.ack_arq, self.txrx_mean, self.analytical_rtt]

    def get_std_by_duplication_each_distance(self):
        return [self.r, self.info_arq, self.std]

    def get_median_by_duplication_each_distance(self):
        return [self.r, self.info_arq, self.med]

    def get_median_by_diffusioncoefficient(self):
        return [self.diffusion_coefficient, self.med]

    def get_csv_data(self):
        return [self.r, self.info_arq, self.ack_arq, self.mean, self.med, self.txrx_mean]

    def get_median_by_diffusioncoefficient_each_distance(self):
        return [self.diffusion_coefficient, self.r, self.med]

    def get_median_by_distance_each_rto(self):
        return [self.r, self.rto_type, self.med]

    def get_mean_by_distance_each_rto(self):
        return [self.r, self.rto_type, self.mean]

    def get_jitter_by_distance_each_rto(self):
        return [self.r, self.rto_type, self.std]

    def get_by_rto_each_distance(self):
        return [self.r, self.rto_type]

    def get_median_by_rto_each_distance(self):
        return [self.r, self.rto_type, self.med]

    def get_jitter_by_rto_each_distance(self):
        return [self.r, self.rto_type, self.std]

    def get_mean_by_rto_each_distance(self):
        return [self.r, self.rto_type, self.mean]

    def get_retransmission_wait_time_by_rto_each_distance(self):
        return [self.r, self.rto_type, self.retransmitWaitTime]

    def get_median_by_retransmission_wait_time_each_distance(self):
        return [self.retransmitWaitTime, self.r, self.med]

    def get_mean_by_retransmission_wait_time_each_distance(self):
        return [self.retransmitWaitTime, self.r, self.mean]

    def get_jitter_by_retransmission_wait_time_each_distance(self):
        return [self.retransmitWaitTime, self.r, self.std]


# 設定(datファイル)の値
class Config:
    def __init__(self, config_file_name):
        self.config_dict = {}
        self.config_header = []
        self.config_file_name = config_file_name
        self.isAdjust = True

        self.init_config_dict()
        self.parse_config_file(config_file_name)

    def init_config_dict(self):
        # 複数ありえるもの
        self.config_dict["moleculeParams"] = []
        self.config_dict["microtubuleParams"] = []

    def parse_config_file(self, config_file_name):
        with open(config_file_name, 'r') as f:
            for line in f:
                # コメントと空行読み飛ばし
                if line[0] == '*' or line[0] == '\n' in line:
                    continue
                # 改行取り除き
                line = line.rstrip()
                # 一番左のスペースで区切る
                key, val = line.split(" ", 1)

                # switchがないのでifでなんとかする
                # transmitterかreceiverの時
                if key in ["transmitter", "receiver"]:
                    self.config_dict[key] = NanoMachine(val)
                    self.config_header.append(key)
                elif key in "intermediateNode":
                    self.config_dict[key] = IntermediateNode(val)
                    self.config_header.append(key)
                elif key in "moleculeParams":
                    molParams = MoleculeParams(val)
                    self.config_dict[key].append(molParams)
                    if not key in self.config_header:
                        self.config_header.append(key)
                    if not molParams.type_of_molecule == MoleculeType.NOISE:
                        if molParams.adaptive_change_number == 0:
                            self.isAdjust = False
                        else:
                            self.adjust_num = molParams.adaptive_change_number
                elif key in "microtubuleParams":
                    self.config_dict[key].append(MicrotubuleParams(val))
                    if not key in self.config_header:
                        self.config_header.append(key)
                elif key in ["probDRail", "stepLengthX", "stepLengthY", "stepLengthZ"]:
                    self.config_dict[key] = float(val)
                    self.config_header.append(key)
                elif key in "outputFile":
                    self.config_dict[key] = val
                    self.config_header.append(key)
                else:
                    self.config_dict[key] = int(val)
                    self.config_header.append(key)

# 入力値
# シミュレーション結果
class InputData:
    def __init__(self, config):
        self.input_file_name = "./result/batch_" + config.config_dict["outputFile"]
        self.steps = []
        self.is_ptime = False
        self.is_coll = config.config_dict["useCollisions"] == 1

        # use if is_ptime is True
        self.info_time = []
        self.info_num = []
        self.ack_time = []
        self.ack_num = []

        self.info_sum_time = 0
        self.info_sum_num = 0
        self.ack_sum_time = 0
        self.ack_sum_num = 0
        self.txrx_mean = 0

        # use if is_coll is True
        self.coll_aa = []
        self.coll_ai = []
        self.coll_an = []
        self.coll_ii = []
        self.coll_in = []

        self.sum_coll_aa = 0
        self.sum_coll_ai = 0
        self.sum_coll_an = 0
        self.sum_coll_ii = 0
        self.sum_coll_in = 0
        self.sum_coll = 0

        self.parse_data()

        if self.is_ptime:
            self.calc_avg()
        if self.is_coll:
            self.calc_coll()

    def parse_data(self):
        self.check_ptime(linecache.getline(self.input_file_name, 1).split(','))
        with open(self.input_file_name, 'r') as f:
            if self.is_ptime and not self.is_coll:
                for line in f:
                    data = line.split(',')
                    self.steps.append(int(data[0]))
                    self.info_time.append(int(data[1]))
                    self.info_num.append(int(data[2]))
                    self.ack_time.append(int(data[3]))
                    self.ack_num.append(int(data[4]))

            elif not self.is_ptime and self.is_coll:
                for line in f:
                    data = line.split(',')
                    self.steps.append(int(data[0]))
                    self.coll_aa.append(int(data[1]))
                    self.coll_ai.append(int(data[2]))
                    self.coll_an.append(int(data[3]))
                    self.coll_ii.append(int(data[4]))
                    self.coll_in.append(int(data[5]))

            elif self.is_ptime and self.is_coll:
                for line in f:
                    data = line.split(',')
                    self.steps.append(int(data[0]))
                    self.info_time.append(int(data[1]))
                    self.info_num.append(int(data[2]))
                    self.ack_time.append(int(data[3]))
                    self.ack_num.append(int(data[4]))
                    self.coll_aa.append(int(data[5]))
                    self.coll_ai.append(int(data[6]))
                    self.coll_an.append(int(data[7]))
                    self.coll_ii.append(int(data[8]))
                    self.coll_in.append(int(data[9]))

            else:
                for line in f:
                    self.steps.append(int(line))

    def check_ptime(self, datas):
        num_of_line = len(datas)

        if self.is_coll:
            num_of_line -= 5

        if num_of_line == 5:
            self.is_ptime = True

    def calc_avg(self):
        self.info_sum_time = sum(self.info_time)
        self.info_sum_num = sum(self.info_num)
        self.ack_sum_time = sum(self.ack_time)
        self.ack_sum_num = sum(self.ack_num)
        self.txrx_mean= float(self.info_sum_time) / self.info_sum_num + float(self.ack_sum_time) / self.ack_sum_num

    def calc_coll(self):
        self.sum_coll_aa = sum(self.coll_aa)
        self.sum_coll_ai = sum(self.coll_ai)
        self.sum_coll_an = sum(self.coll_an)
        self.sum_coll_ii = sum(self.coll_ii)
        self.sum_coll_in = sum(self.coll_in)
        self.sum_coll = self.sum_coll_aa + self.sum_coll_ai + self.sum_coll_an + self.sum_coll_ii + self.sum_coll_in

# 出力値
# シミュレーション結果の統計
class OutputData:
    def __init__(self, config, input_data):
        self.steps = np.sort(input_data.steps)
        self.count = len(self.steps)
        self.minimum = int(self.steps[0])
        self.maximum = int(self.steps[-1])
        self.range_num = int(self.maximum / 1000) * 10
        if self.range_num == 0:
            self.range_num = int(self.maximum / 100)
        self.var = np.var(self.steps)
        self.std = np.std(self.steps)
        self.med = median(self.steps)
        self.mean = mean(self.steps)

        self.is_ptime = input_data.is_ptime

        self.txrx_mean = 0
        if self.is_ptime:
            self.txrx_mean = input_data.txrx_mean

        self.analytical_model = AnalyticalModel(config)

    def to_array(self):
        if self.is_ptime:
            return [self.count, self.minimum, self.maximum, self.range_num,
            self.var, self.std, self.med, self.mean, self.txrx_mean]
        else:
            return [self.count, self.minimum, self.maximum, self.range_num,
            self.var, self.std, self.med, self.mean]

# グラフを書くためのデータ
class PlotData:
    def __init__(self, output_data):
        self.range_num = output_data.range_num
        self.plot_range = []
        self.prob = []
        self.cum_prob = []

        self.range_num_by_duplication = 10000
        self.plot_range_by_duplication = []
        self.cum_prob_by_duplication = []

        self.plt_data = self.make_plt_data(output_data.steps)
        self.plt_data_by_duplication = self.make_plt_data_by_duplication(output_data.steps)
        self.calc_prob(len(output_data.steps))
        self.calc_prob_by_duplication(len(output_data.steps))

    def make_plt_data(self, datas):
        plt_data = []
        data_arr = []
        head = 0
        tail = self.range_num
        self.plot_range.append([head, tail - 1])
        i = 0
        while True:
            if datas[i] < tail:
                data_arr.append(datas[i])
                i += 1
                if i == len(datas):
                    break
            else:
                plt_data.append(data_arr)
                head += self.range_num
                tail += self.range_num
                self.plot_range.append([head, tail - 1])
                data_arr = []

        plt_data.append(data_arr)
        return plt_data

    def make_plt_data_by_duplication(self, datas):
        plt_data = []
        data_arr = []
        head = 0
        tail = self.range_num_by_duplication
        self.plot_range_by_duplication.append([head, tail - 1])
        i = 0
        while True:
            if datas[i] < tail:
                data_arr.append(datas[i])
                i += 1
                if i == len(datas):
                    break
            else:
                plt_data.append(data_arr)
                head += self.range_num_by_duplication
                tail += self.range_num_by_duplication
                self.plot_range_by_duplication.append([head, tail - 1])
                data_arr = []

        plt_data.append(data_arr)
        return plt_data

    def calc_prob(self, step_num):
        prob = 0.0
        cum_prob = 0.0
        for i in range(len(self.plt_data)):
            head, tail = self.plot_range[i]
            prob = float(len(self.plt_data[i])) / step_num * 100.0
            cum_prob += prob
            self.prob.append(prob)
            self.cum_prob.append(cum_prob)

    def calc_prob_by_duplication(self, step_num):
        prob = 0.0
        cum_prob = 0.0
        for i in range(len(self.plt_data_by_duplication)):
            head, tail = self.plot_range_by_duplication[i]
            prob = float(len(self.plt_data_by_duplication[i])) / step_num * 100.0
            cum_prob += prob
            self.cum_prob_by_duplication.append(cum_prob)

# 推論モデルのデータ
class AnalyticalModel:
    def __init__(self, config):
        rx = config.config_dict["receiver"].center_position
        tx = config.config_dict["transmitter"].center_position

        step_length = config.config_dict["stepLengthX"]
        if step_length in DNA_DICT.keys():
            self.D = DNA_DICT[step_length]["dc"]
        else:
            self.D = 0.5

        self.r = self.calc_distance(tx, rx)
        self.L = int(config.config_dict["mediumDimensionX"] / 2)
        self.l = (int(config.config_dict["receiver"].size) * 2 - 1) / 2

        self.info_winf = 0.0
        self.ack_winf = 0.0

        molecule_params = config.config_dict["moleculeParams"]
        for molecule_param in molecule_params:
            # from transmitter to receiver
            if molecule_param.type_of_molecule == MoleculeType["INFO"]:
                if molecule_param.type_of_movement == MovementType["PASSIVE"]:
                    self.info_winf = self.calc_passive_rtt()
                else:
                    passive_winf = self.calc_passive_rtt()
                    self.info_winf = self.calc_active_rtt(passive_winf)

            # from receiver to transmitter
            elif molecule_param.type_of_molecule == MoleculeType["ACK"]:
                if molecule_param.type_of_movement == MovementType["PASSIVE"]:
                    self.ack_winf = self.calc_passive_rtt()
                else:
                    passive_winf = self.calc_passive_rtt()
                    self.ack_winf = self.calc_active_rtt(passive_winf)

        self.rtt = self.info_winf + self.ack_winf

    def calc_distance(self, tx, rx):
        x = (tx.x - rx.x) ** 2
        y = (tx.y - rx.y) ** 2
        z = (tx.z - rx.z) ** 2
        return math.sqrt(x + y + z)

    def calc_passive_rtt(self):
        winf = (self.r - self.l) * (2 * self.L ** 3 - self.l * self.r ** 2 - self.l ** 2 * self.r) / (2 * self.D * self.l * self.r)
        return winf

    def calc_active_rtt(self, passive_winf):
        V1 = 2 * 2 * self.r
        V = 4.0 / 3.0 * math.pi * self.L ** 3
        p = V1 / V
        va = 1
        vp = self.r / passive_winf
        ve = p * va + (1.0 - p) * vp
        winf = self.r / ve
        return winf

    def to_array(self):
        return [self.D, self.r, self.L, self.l, self.rtt]

# 3次元位置情報
class Position:
    def __init__(self, args):
        self.x = int(args[0])
        self.y = int(args[1])
        self.z = int(args[2])

    def to_string(self):
        return "({}, {}, {})".format(self.x, self.y, self.z)

# transmitterとreceiverのクラス
class NanoMachine:
    def __init__(self, val):
        args = self.parse_val(val)
        self.center_position = Position(args[0:3])
        self.size = int(args[3])
        self.release_position = Position(args[4:7])

    def parse_val(self, val):
        # めんどくさい文字列なので非英数文字で区切って数値だけの配列に
        val = [i for i in re.split(r"[,( )]", val) if i != '']
        return val

    def to_string(self):
        return "{} {} {}".format(self.center_position.to_string(), str(self.size), self.release_position.to_string())

# 中間ノードクラス
class IntermediateNode:
    def __init__(self, val):
        args = self.parse_val(val)
        self.center_position = Position(args[0:3])
        self.size = int(args[3])
        self.info_release_position = Position(args[4:7])
        self.ack_release_position = Position(args[7:10])

    def parse_val(self, val):
        val = [i for i in re.split(r"[,( )]", val) if i != '']
        return val

    def to_string(self):
        return "{} {} {} {}".format(self.center_position.to_string(), str(self.size), self.info_release_position.to_string(), self.ack_release_position.to_string())

# 分子の情報のクラス
class MoleculeParams:
    def __init__(self, val):
        args = self.parse_val(val)
        self.num_of_molecules = int(args[0])
        self.type_of_molecule = MoleculeType[args[1]]
        if self.type_of_molecule != MoleculeType["NOISE"]:
            self.type_of_movement = MovementType[args[2]]
            self.adaptive_change_number = int(args[3])
            self.size = float(args[4])
        else:
            self.size = float(args[2])

    def parse_val(self, val):
        val = [i for i in re.split(r"[ ]", val) if i != '']
        return val

    def to_string(self):
        if self.type_of_molecule != MoleculeType["NOISE"]:
            return "{} {} {} {}".format(str(self.num_of_molecules), self.type_of_molecule.name, self.type_of_movement.name, str(self.adaptive_change_number))
        else:
            return "{} {}".format(str(self.num_of_molecules), self.type_of_molecule.name)

# マイクロチューブのクラス
class MicrotubuleParams:
    def __init__(self, val):
        args = self.parse_val(val)
        self.start_position = Position(args[0:3])
        self.end_position = Position(args[3:6])

    def parse_val(self, val):
        val = [i for i in re.split(r"[,( )]", val) if i != '']
        return val

    def to_string(self):
        return "{} {}".format(self.start_position.to_string(), self.end_position.to_string())
