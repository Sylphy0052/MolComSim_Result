DD = 24 * 60 * 60
HH = 60 * 60
MM = 60
SS = 1

class DataTime:
    def __init__(self, file_name, seconds, count):
        self.file_name = file_name
        self.seconds = int(seconds)
        self.count = int(count)
        self.second_by_one = float(self.seconds) / self.count

def show_one_sim_time(arr):
    for data in arr:
        print("{}: {}".format(data.file_name, data.second_by_one))

def show_all_time(arr):
    all_time = 0.0
    for data in arr:
        all_time += data.second_by_one

    all_time *= 1000

    dd = int(all_time / DD)
    all_time %= DD
    hh = int(all_time / HH)
    all_time %= HH
    mm = int(all_time / MM)
    all_time %= MM
    ss = int(all_time)

    print("All Time: {}d:{}h:{}m:{}s".format(dd,hh,mm,ss))

def main():
    data_arr = []
    with open("./result/time.txt", 'r') as f:
        for line in f:
            if len(line.split(',')) != 3:
                continue
            file_name, seconds, count = line.split(',')
            data_arr.append(DataTime(file_name, seconds, count))

    show_one_sim_time(data_arr)
    show_all_time(data_arr)

if __name__ == '__main__':
    main()
