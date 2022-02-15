import numpy as np
import sys

def GetRSA(AA, ASA):
    look_up = {
        'A': 110.2,
        'D': 144.1,
        'C': 140.4,
        'E': 174.7,
        'F': 200.7,
        'G': 78.7,
        'H': 181.9,
        'I': 185.0,
        'K': 205.7,
        'L': 183.1,
        'M': 200.1,
        'N': 146.4,
        'P': 141.9,
        'Q': 178.6,
        'R': 229.0,
        'S': 117.2,
        'T': 138.7,
        'V': 153.7,
        'W': 240.5,
        'Y': 213.7
    }
    if AA not in look_up:
        AA = 'G'
    return min(1.0, ASA/look_up[AA])

def extract(input, src, dest):
    BG_RSA = 0.3213596246201607
    with open(input) as pdb_list:
        pdb = pdb_list.readline().rstrip()
        while(pdb):
        # 3672 * atom_number 
            print(pdb)
            with open(src + pdb) as raw_label:
                label = []
                count = 1
                line = raw_label.readline().lstrip()
                while not line[0] == '#':
                    line = raw_label.readline().strip()
                line = raw_label.readline()
                while line:
                    # items = [l for l in line.split() if not l == ' ']
                    # print(items)
                    AA = line[13]
                    if not AA == '!':
                        idx = int(line[5:10])
                        if not count == idx:
                            print(count, idx)
                        while count < idx:
                            label.append(BG_RSA)
                            count += 1
                        ASA = int(line[35:38])
                        RSA = GetRSA(AA,ASA)
                        count += 1
                        label.append(RSA)
                    else:
                        pass
                        # label.append(BG_RSA)
                    # count += 1
                    line = raw_label.readline()

                np.save(dest + pdb, np.array(label))
            pdb = pdb_list.readline().rstrip()
    # return BG_RSA, count

mode = sys.argv[1]
if mode == "train":
    input = "/nfs/amino-home/panqp/protein_stru/repo/EAGERER/datasets/tr3672.txt"
    src = "/nfs/amino-home/panqp/protein_stru/repo/EAGERER/data/train_label_raw/"
    dest = "/nfs/amino-home/panqp/protein_stru/repo/EAGERER/data/train_label/"
elif mode == "validation":
    input = "/nfs/amino-home/panqp/protein_stru/repo/EAGERER/datasets/val918.txt"
    src = "/nfs/amino-home/panqp/protein_stru/repo/EAGERER/data/validation_label_raw/"
    dest = "/nfs/amino-home/panqp/protein_stru/repo/EAGERER/data/validation_label/"
else:
    input = "/nfs/amino-home/panqp/protein_stru/repo/EAGERER/datasets/ts1199.txt"
    src = "/nfs/amino-home/panqp/protein_stru/repo/EAGERER/data/test_label_raw/"
    dest = "/nfs/amino-home/panqp/protein_stru/repo/EAGERER/data/test_label/"

extract(input, src, dest)

 


