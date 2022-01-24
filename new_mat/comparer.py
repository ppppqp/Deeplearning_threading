import sys

from progress import ProgressBar
if __name__ == "__main__":

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    f1_seq = []
    f2_seq = []
    with open(file1) as f1:
        line = f1.readline()
        f1_seq = line.split(' ')
    with open(file2) as f2:
        line = f2.readline()
        f2_seq = line.split(' ')
    print("file 1 length:", len(f1_seq))
    print("file 2 length:", len(f2_seq))

    # progress = ProgressBar(min(len(f1_seq), len(f2_seq)),"comparing:")
    for idx in range(min(len(f1_seq), len(f2_seq))):
        # progress.current += 1
        # progress()
        if ':' in f1_seq[idx] and ':' in f2_seq[idx]:
            f1_num = float(f1_seq[idx].split(':')[1])
            f2_num = float(f2_seq[idx].split(':')[1])
            if not abs(f1_num - f2_num) < 0.000001:
                print("DIFFERENT: idx: ", idx, " f1:", f1_num, " f2:", f2_num)
            else:
                print("SAME: idx: ", idx, " f1:", f1_num, " f2:", f2_num)
        else:
            print("f1",f1_seq[idx])
            print("f2",f2_seq[idx])
    # progress.done()
    print("Comparision finished")