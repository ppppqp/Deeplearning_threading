#!/nfs/amino-home/zhng/local_library/anaconda3/bin/python
import getopt
import os
from progress import ProgressBar
import sys
from typing import Sequence
import math
from utils import get_stru_files
ResName = "ARNDCQEGHILKMFPSTWYV"
factor1 = [ -0.591, 1.538, 0.945, 1.05, -1.343, 0.931, 1.357, -0.384, 0.336, -1.239, -1.019, 1.831, -0.663, -1.006, 0.189, -0.228, -0.032, -0.595, 0.26, -1.337]
factor2 = [ -1.302, -0.055, 0.828, 0.302, 0.465, -0.179, -1.453, 1.652, -0.417, -0.547, -0.987, -0.561, -1.524, -0.59, 2.081, 1.399, 0.326, 0.009, 0.83, -0.279 ]
factor3 = [ -0.733, 1.502, 1.299, -3.656, -0.862, -3.005, 1.477, 1.33, -1.673, 2.131, -1.505, 0.533, 2.219, 1.891, -1.628, -4.76, 2.213, 0.672, 3.097, -0.544 ]
factor4 = [1.57, 0.44, -0.169, -0.259, -1.02, -0.503, 0.113, 1.045, -1.474, 0.393, 1.266, -0.277, -1.005, -0.397, 0.421, 0.67, 0.908, -2.128, -0.838, 1.242]
factor5 = [ -0.146, 2.897, 0.933, -3.242, -0.255, -1.853, -0.837, 2.064, -0.078, 0.816, -0.912, 1.648, 1.212, 0.412, -1.392, -2.647, 1.313, -0.184, 1.512, -1.262 ]
BLOSUM62aa = "ARNDCQEGHILKMFPSTWYVBZX*"
BLOSUM62 = [[4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,-1,0,-4], 
            [-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,-1,0,-1,-4],
            [-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,3,0,-1,-4],
            [-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,4,1,-1,-4],
            [0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4],
            [-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0,3,-1,-4],
            [-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4],
            [0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,-1,-2,-1,-4],
            [-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,0,0,-1,-4],
            [-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,-3,-3,-1,-4],
            [-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,-4,-3,-1,-4],
            [-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,0,1,-1,-4],
            [-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,-3,-1,-1,-4],
            [-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,-3,-3,-1,-4],
            [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,-2,-1,-2,-4],
            [1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,0,0,0,-4],
            [0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,-1,-1,0,-4],
            [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,-4,-3,-2,-4],
            [-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,-3,-2,-1,-4],
            [0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4,-3,-2,-1,-4],
            [-2,-1,3,4,-3,0,1,-1,0,-3,-4,0,-3,-3,-2,0,-1,-4,-3,-3,4,1,-1,-4],
            [-1,0,0,1,-3,3,4,-2,0,-3,-3,1,-1,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4],
            [0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,0,0,-2,-1,-1,-1,-1,-1,-4],
            [-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,1]]
class Residue:
    def __init__(self, protein_name, residue_name, position, segment):
        self.protein_name = protein_name
        self.residue_name = residue_name
        self.position = position
        self.segment = segment
        self.PSSM_con = []
        self.PSSM_prob = []
        self.AA_indexs = []
        self.second_structure = ' '
        self.second_stru_prob = []
        self.HH_list = []
        self.disorder = 0
        self.ASA = 0
        self.RSA = 0
        self.Phi = 0
        self.Psi = 0 
        self.factor1 = 0
        self.factor2 = 0
        self.factor3 = 0
        self.factor4 = 0
        self.factor5 = 0
        self.label = 0
        self.alpha_HSE_d = 0
        self.alpha_HSE_u = 0
        self.beta_HSE_d = 0
        self.beta_HSE_u = 0




class Feature: 
    def __init__(self) -> None:
        self.terminal_factor = []
        self.PWAA = []
        self.EBGW = []
        self.CKSAAP = []
        self.AAindex_mean = []
        self.AAindex_pair_mean = []
        self.factors = []
        self.disorder = []
        self.disorderCTD = []

        self.ss_win = []
        self.ss_prob = []
        self.ss_percent = []
        self.ss_seg_cnt = []
        self.ss_max_seg_cnt = []
        self.RSA_aa_mean = []
        self.RSA_ss_mean = []
        self.RSA_win_mean = 0
        self.RSA_win = []
        self.Phi_aa_mean = []
        self.Phi_ss_mean = []
        self.Phi_win_mean = 0
        self.Phi_win = []
        self.Psi_aa_mean = []
        self.Psi_ss_mean = []
        self.Psi_win_mean = 0
        self.Psi_win = []
        self.PSSM_con_win = []
        self.PSSM_prob_win = []
        self.PSSM_con_AA_mean = []
        self.PSSM_prob_AA_mean = []

        self.PSSM_con_row_max = []
        self.PSSM_con_row_min = []
        self.PSSM_con_row_mean = []

        self.PSSM_con_col_max = []
        self.PSSM_con_col_min = []
        self.PSSM_con_col_mean = [] 

        self.PSSM_prob_row_max = []
        self.PSSM_prob_col_max = []
        self.PSSM_prob_col_mean = []

        self.PSSM_entropy = []
        self.KNN = []

        self.HH_frequency = []
        self.HH_entropy = []
        self.alpha_HSE_d = []
        self.alpha_HSE_u = []
        self.beta_HSE_d = []
        self.beta_HSE_u = []


class KNN_node:
    def __init__(self):
        self.label = 0
        self.segment = ' '
        self.distance = 0.0
        self.index = 0

def KNN_compare(node):
    return node.distance

def KNN_distance(seg1, seg2):
    dist = 0
    blosum_max = 11
    blosum_min = -4
    # for i in range(21):
    #     for j in range(21):
    #         blosum_max = max(BLOSUM62[i][j], blosum_max)
    #         blosum_min = min(BLOSUM62[i][j], blosum_min)
    # print(blosum_max)
    # print(blosum_min)
    if len(seg1) == len(seg2):
        sim = []
        for i in range(len(seg1)):
            i1 = i2 = 0
            if seg1[i] in BLOSUM62aa:
                i1 = BLOSUM62aa.index(seg1[i])
            else:
                i1 = 20    
            if seg2[i] in BLOSUM62aa:
                i2 = BLOSUM62aa.index(seg2[i])
            else:
                i2 = 20
            sim.append((BLOSUM62[i1][i2] - blosum_min) / float(blosum_max - blosum_min))
        dist = 1 - sum(sim) / len(seg1)
    else:
        dist = -100000
    return dist

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
    return min(1.0, ASA/look_up[AA])

def format_num(num):
    return str(num)[0:17]

def read_feature_files(fasta_file, AA_file, feature_dir, HH_dir, hse_dir, window_length):

    res_seq = []
    with open(fasta_file) as f:
        lines = f.readlines()
        name = lines[0][:-1] # get rid of the last \n 
        seq = lines[1][:-1]
        print(seq)
        label = ""
        label_file = False
        f.seek(0)
        all_text = f.read().strip()
        if all_text[-1].isdigit():
            label = lines[2]
            label_file = True
        format_seq = seq
        for ii in range((window_length-1) // 2):
            format_seq = '*' + format_seq + '*'
        for ii in range(len(seq)):
            half_win_len = (window_length - 1)//2
            res_seq.append(Residue(name, seq[ii], ii+1, format_seq[ii:ii+window_length]))
            if label_file:
                res_seq[-1].label = int(label[ii])
        

    fasta_file = fasta_file.split('/')[-1].split('.')[0]
    # read fasta_file and construct res_seq
    p = 0
    with open(os.path.join(feature_dir, fasta_file + ".diso")) as f:
        # read diso
        for ii in range(5): f.readline()

        line = f.readline()
        line = line.rstrip()
        while line:
            if len(line) > 10:
                disoname = line[6]
                if disoname == res_seq[p].residue_name:
                    if line[8] == '*':
                        res_seq[p].disorder = 1
                    else:
                        res_seq[p].disorder = 0
                else:
                    print(disoname, res_seq[p].residue_name)
                    raise Exception("disorder: not matching")
                p += 1
            line = f.readline() # read next line
            if line: line = line.rstrip()
    p = 0

    with open(os.path.join(feature_dir, fasta_file + ".ss2")) as f:
        # read ss
        for ii in range(2): f.readline()
        line = f.readline()
        line = line.rstrip()
        while line:
            if len(line) > 10:
                ss_name = line[5]
                if ss_name == res_seq[p].residue_name:
                    temp_ss_val = line.split()
                    res_seq[p].second_structure = temp_ss_val[2]
                    res_seq[p].second_stru_prob.extend([float(i) for i in temp_ss_val[3:6]])
                else:
                    raise Exception('ss: not matching')
                p += 1
            line = f.readline() # read next line
            if line: line = line.rstrip()
    p = 0
    with open(os.path.join(feature_dir, fasta_file + ".spXout")) as f:
        # read spineX
        f.readline()
        line = f.readline()
        line = line.rstrip()
        while line:
            if len(line) > 10:
                spineX_name = line[13]
                if spineX_name == res_seq[p].residue_name:
                    spineX_str = line.split()
                    res_seq[p].Phi = float(spineX_str[3])
                    res_seq[p].Psi = float(spineX_str[4])
                    res_seq[p].ASA = float(spineX_str[10])
                    res_seq[p].RSA = GetRSA(res_seq[p].residue_name, res_seq[p].ASA)
                else:
                    raise Exception('spineX: not matching')
                p += 1
            line = f.readline() # read next line
            if line: line = line.rstrip()
    p = 0
    with open(os.path.join(feature_dir, fasta_file + ".mat")) as f:
        # read psiblast
        for ii in range(3): f.readline()
        line = f.readline()
        line = line.rstrip()
        while line:
            if len(line) > 10 and line.lstrip()[0].isdigit():
                psi_name = line[6]
                if psi_name == res_seq[p].residue_name:
                    line = line[8:]
                    items = line.split()
                    for k1 in range(20):
                        res_seq[p].PSSM_con.append(int(items[k1]))
                    for k2 in range(20,40):
                        res_seq[p].PSSM_prob.append(int(items[k2]) / 100.0)
                else:
                    raise Exception('psiblast: not matching')
                p += 1
            line = f.readline() # read next line
            if line: line = line.rstrip()
    p = 0        
    with open(os.path.join(HH_dir, fasta_file + ".aln")) as f:
        # read hhaln
        line = f.readline()
        line = line.rstrip()
        while line:
            for ii in range(len(line)):
                res_seq[ii].HH_list.append(line[ii])
            line = f.readline() # read next line
            if line: line = line.rstrip()
    p = 0        
    with open(os.path.join(hse_dir, fasta_file + ".hsa2")) as f:
        # read hsa
        line = f.readline()
        line = f.readline()
        line = line.rstrip()
        while line:
            items = line.split('\t')
            if items[1] == res_seq[p].residue_name: 
                res_seq[p].alpha_HSE_d = float(items[3])
                res_seq[p].alpha_HSE_u = float(items[4])
            else:
                print(items[1], res_seq[p].residue_name)
                raise Exception('hsa: not matching')
            line = f.readline() # read next line
            if line: line = line.rstrip()
            p += 1
    p = 0
    with open(os.path.join(hse_dir, fasta_file + ".hsb2")) as f:
        # hsb
        f.readline()
        line = f.readline()
        line = line.rstrip()
        while line:
            items = line.split('\t')
            if res_seq[p].residue_name == items[1]:
                res_seq[p].beta_HSE_d = items[3]
                res_seq[p].beta_HSE_u = items[4]
            else:
                raise Exception('hsb: not match')
            line = f.readline() # read next line
            if line: line = line.rstrip()
            p += 1
    p = 0
    with open(AA_file) as f:
        AAdatabase = []
        f.readline()
        line = f.readline()
        line = line.rstrip()
        while line:
            if len(line) > 0:
                aastr1 = line.split('\t')[1]
                aastr2 = aastr1.split(',')
                AA = [float(str) for str in aastr2]
                AAdatabase.append(AA)
            line = f.readline() # read next line
            if line: line = line.rstrip()
        for node in res_seq:
            if node.residue_name in ResName:
                aa_position = ResName.index(node.residue_name)
                for aalist in AAdatabase:
                    node.AA_indexs.append(aalist[aa_position])
                node.factor1 = factor1[aa_position]
                node.factor2 = factor2[aa_position]
                node.factor3 = factor3[aa_position]
                node.factor4 = factor4[aa_position]
                node.factor5 = factor5[aa_position]
    return res_seq
def compute_feature(res_seq, window_length, sample_file):
    feature_seq = []
    half_win_len = (window_length - 1) // 2
    res_seq_len = len(res_seq)
    for ii in range(res_seq_len):
        feature_seq.append(Feature())
    # terminal factor
    for i in range(res_seq_len):
        if i - half_win_len < 0:
            feature_seq[i].terminal_factor.extend([1,0,0])
        elif i + half_win_len >= res_seq_len:
            feature_seq[i].terminal_factor.extend([0,0,1])
        else:
            feature_seq[i].terminal_factor.extend([0,1,0])
    
    # position weight amino acide composition
    for i in range(res_seq_len):
        sub_seq = []
        for j in range(i - half_win_len, i + half_win_len + 1):
            if j < 0 or j >= res_seq_len:
                sub_seq.append('*')
            else:
                sub_seq.append(res_seq[j].residue_name)
    
        # get feature
        for j in range(0, 20):
            c = 0
            for k in range(len(sub_seq)):
                if sub_seq[k] == ResName[j]:
                    temp = k - half_win_len
                    c += (temp + abs(temp)/half_win_len)
            initialPWAA = c / (half_win_len * (half_win_len + 1))
            feature_seq[i].PWAA.append(initialPWAA + 1.0 /2.0 - 1.0 / (2.0 * window_length))

    hydrophobic_group = "AFGILMPVW"
    polar_group = "CNQSTY"
    pos_charged_group = "HKR"
    neg_charged_group = "DE"

    J = 5
    L = []
    for i in range(J):
        L.append(int(round((i+1.0)*window_length/J)))
        #17/5
        #17*2/5
        #17*3/5
        #17*4/5
        #17*5/5
    for i in range(res_seq_len):
        aa = ''
        H1 = []
        H2 = []
        H3 = []
        for j in range(i - half_win_len, i + half_win_len + 1):
            if j < 0 or j >= res_seq_len:
                aa = '*'
                H1.append(0)
                H2.append(0)
                H3.append(0)
            else:
                aa = res_seq[j].residue_name
                if aa in hydrophobic_group:
                    H1.append(1)
                    H2.append(1)
                    H3.append(1)             
                elif aa in polar_group:
                    H1.append(1)
                    H2.append(0)
                    H3.append(0)                    
                elif aa in pos_charged_group:
                    H1.append(0)
                    H2.append(1)
                    H3.append(0)
                elif aa in neg_charged_group:
                    H1.append(0)
                    H2.append(0)
                    H3.append(1)
                else:
                    H1.append(0)
                    H2.append(0)
                    H3.append(0)
        for j in range(J):
            temp = len([p for p in H1[0:L[j]] if p == 1])
            feature_seq[i].EBGW.append(float(temp)/L[j])
        for j in range(J):
            temp = len([p for p in H2[0:L[j]] if p == 1])
            feature_seq[i].EBGW.append(float(temp)/L[j])
        for j in range(J):
            temp = len([p for p in H3[0:L[j]] if p == 1])
            feature_seq[i].EBGW.append(float(temp)/L[j])
    
    # considering residue pairs CKSAAP encoding scheme
    for i in range(res_seq_len):
        sub_seq = []
        for j in range(i - half_win_len, i + half_win_len + 1):
            if j >= 0 and j < res_seq_len:
                sub_seq.append(res_seq[j].residue_name)
    
        # get paris number
        pairs_num = []
        for j in range(400):
            pairs_num.append(0)
        for j in range(len(sub_seq)-1):
            if sub_seq[j] in ResName and sub_seq[j+1] in ResName:
                pairs_1_index = ResName.index(sub_seq[j])
                pairs_2_index = ResName.index(sub_seq[j+1])
                pairs_num[pairs_1_index * 20 + pairs_2_index] += 1

        # get feature   
        for j in range(400):
            feature_seq[i].CKSAAP.append(pairs_num[j]/ (window_length - 1.0))

    # aa index mean
    for i in range(res_seq_len):
        for k in range(len(res_seq[0].AA_indexs)):
            AAindex_sum = 0.0
            AAindex_N = 0
            for j in range(i - half_win_len, i + half_win_len + 1):
                if j >= 0 and j < res_seq_len and len(res_seq[j].AA_indexs) > 0:
                    AAindex_sum += + res_seq[j].AA_indexs[k]
                    AAindex_N += 1
            feature_seq[i].AAindex_mean.append(AAindex_sum / AAindex_N)

    # AAindexPairMean
    for i in range(res_seq_len):
        for k in range(len(res_seq[0].AA_indexs)):
            AAindex_pair_sum = 0.0
            AAindex_N = 0
            for j in range(i - half_win_len, i + half_win_len):
                  if j >= 0 and j < res_seq_len - 1 and len(res_seq[j].AA_indexs) > 0 and len(res_seq[j+1].AA_indexs) > 0:
                    temp = res_seq[j].AA_indexs[k] * res_seq[j+1].AA_indexs[k]
                    AAindex_pair_sum += temp
                    AAindex_N += 1
            feature_seq[i].AAindex_pair_mean.append(AAindex_pair_sum / AAindex_N)

    # Factors
    for i in range(res_seq_len):
        for j in range(i - half_win_len, i + half_win_len + 1):
            if j < 0 or j >= res_seq_len:
                feature_seq[i].factors.extend([0,0,0,0,0])
            else:
                feature_seq[i].factors.append(res_seq[j].factor1)
                feature_seq[i].factors.append(res_seq[j].factor2)
                feature_seq[i].factors.append(res_seq[j].factor3)
                feature_seq[i].factors.append(res_seq[j].factor4)
                feature_seq[i].factors.append(res_seq[j].factor5)
    # Disorder
    for i in range(res_seq_len):
        for j in range(i - half_win_len, i + half_win_len + 1):
            if j < 0 or j >= res_seq_len:
                feature_seq[i].disorder.append(0)
            else:
                feature_seq[i].disorder.append(res_seq[j].disorder)
    
    # DisorderCTD
    for i in range(res_seq_len):
        # get subsequence
        sub_seq = []
        for j in range(i - half_win_len, i + half_win_len + 1):
            if j >= 0 and j < res_seq_len:
                sub_seq.append(res_seq[j].disorder)
        
        # get C - composition
        disorderN = len([p for p in sub_seq if p == 1])
        orderN = len([p for p in sub_seq if p == 0])
        feature_seq[i].disorderCTD .append(disorderN/float(len(sub_seq)))
        feature_seq[i].disorderCTD .append(orderN/float(len(sub_seq)))

        # get T - transition
        transitionN = 0
        for j in range(len(sub_seq) - 1):
            if sub_seq[j] != sub_seq[j+1]:
                transitionN += 1
        feature_seq[i].disorderCTD.append(transitionN/(len(sub_seq)-1.0))

        # get D - distribution
        tempDisorder = 0
        d1 = 0
        d2 = 0
        d3 = 0
        d4 = 0
        d5 = 0
        for j in range(len(sub_seq)):
            if sub_seq[j] == 1:
                tempDisorder += 1
            else:
                continue

            if tempDisorder == 1:
                d1 = j + 1.0
            if tempDisorder == int(disorderN * 0.25 + 0.5):
                d2 = j + 1.0
            if tempDisorder == int(disorderN * 0.5 + 0.5):
                d3 = j + 1.0            
            if tempDisorder == int(disorderN * 0.75 + 0.5):
                d4 = j + 1.0
            if tempDisorder == disorderN:
                d5 = j + 1.0

        feature_seq[i].disorderCTD.extend([
            d1 / len(sub_seq),
            d2 / len(sub_seq),
            d3 / len(sub_seq),
            d4 / len(sub_seq),
            d5 / len(sub_seq)
        ])

        tempOrder = 0
        d6 = 0
        d7 = 0
        d8 = 0
        d9 = 0
        d10 = 0
        for j in range(len(sub_seq)):
            if sub_seq[j] == 0:
                tempOrder += 1
            else:
                continue

            if tempOrder == 1:
                d6 = j + 1.0
            if tempOrder == int(orderN * 0.25 + 0.5):
                d7 = j + 1.0
            if tempOrder == int(orderN * 0.5 + 0.5):
                d8 = j + 1.0            
            if tempOrder == int(orderN * 0.75 + 0.5):
                d9 = j + 1.0
            if tempOrder == orderN:
                d10 = j + 1.0

        feature_seq[i].disorderCTD.extend([
            d6 / len(sub_seq),
            d7 / len(sub_seq),
            d8 / len(sub_seq),
            d9 / len(sub_seq),
            d10 / len(sub_seq)
        ])
    # SSwindow
    for i in range(res_seq_len):
        for j in range(i - half_win_len, i + half_win_len + 1):
            if j < 0 or j >= res_seq_len:
                feature_seq[i].ss_win.extend([0,0,0])
            else:
                if res_seq[j].second_structure == 'H':
                    feature_seq[i].ss_win.extend([1,0,0])
                elif res_seq[j].second_structure == 'C':
                    feature_seq[i].ss_win.extend([0,1,0])
                elif res_seq[j].second_structure == 'E':
                    feature_seq[i].ss_win.extend([0,0,1])
                else:
                    feature_seq[i].ss_win.extend([0,0,0])
    
    # SSprob, SSpercent
    for i in range(res_seq_len):
        sub_char = []
        sub_Cprob = []
        sub_Hprob = []
        sub_Eprob = []
        for j in range(i - half_win_len, i + half_win_len + 1):
            if j >= 0 and j < len(res_seq):
                sub_char.append(res_seq[j].second_structure)
                sub_Cprob.append(res_seq[j].second_stru_prob[0])
                sub_Hprob.append(res_seq[j].second_stru_prob[1])
                sub_Eprob.append(res_seq[j].second_stru_prob[2])

        # SSprob
        feature_seq[i].ss_prob.append(sum(sub_Cprob)/len(sub_Cprob))
        feature_seq[i].ss_prob.append(sum(sub_Hprob)/len(sub_Hprob))
        feature_seq[i].ss_prob.append(sum(sub_Eprob)/len(sub_Eprob))
        # SSpercent
        feature_seq[i].ss_percent.append(len([p for p in sub_char if p == 'C'])/len(sub_char))
        feature_seq[i].ss_percent.append(len([p for p in sub_char if p == 'H'])/len(sub_char))
        feature_seq[i].ss_percent.append(len([p for p in sub_char if p == 'E'])/len(sub_char))

        # get seg
        sub_char.append('*')
        segC = []
        segH = []
        segE = []
        count = 1
        temp = sub_char[0]
        for j in range(1, len(sub_char)):
            if sub_char[j] == temp:
                count += 1
            elif sub_char[j] != temp and count > 1:
                if temp == 'H':
                    segH.append(count)
                elif temp == 'E':
                    segE.append(count)
                else:
                    segC.append(count)
                count = 1
            temp = sub_char[j]
        
        # ss segcount
        for j in range(2, window_length + 1):
            a1 = len([p for p in segH if p >= j])
            b1 = len(segH) + len(segE) + len(segC)
            feature_seq[i].ss_seg_cnt.append(a1 / float(b1))

            a2 = len([p for p in segE if p >= j])
            b2 = len(segH) + len(segE) + len(segC)
            feature_seq[i].ss_seg_cnt.append(a2 / float(b2))
            
            a3 = len([p for p in segC if p >= j])
            b3 = len(segH) + len(segE) + len(segC)
            feature_seq[i].ss_seg_cnt.append(a3 / float(b3))
        
        # ss max seg cnt
        if len(segH) > 0:
            feature_seq[i].ss_max_seg_cnt.append(max(segH)/float(window_length))
        else:
            feature_seq[i].ss_max_seg_cnt.append(0)

        if len(segE) > 0:
            feature_seq[i].ss_max_seg_cnt.append(max(segE)/float(window_length))
        else:
            feature_seq[i].ss_max_seg_cnt.append(0)

        if len(segC) > 0:
            feature_seq[i].ss_max_seg_cnt.append(max(segC)/float(window_length))
        else:
            feature_seq[i].ss_max_seg_cnt.append(0)

    # RSA, PHI, PSI
    for i in range(res_seq_len):
        sub_ss = []
        sub_seq = []
        sub_rsa = []
        sub_phi = []
        sub_psi = []
        # RSAwindow, PhiWindow, PsiWindow
        for j in range(i - half_win_len, i + half_win_len + 1):
            if j >= 0 and j < res_seq_len:
                sub_ss.append(res_seq[j].second_structure)
                sub_seq.append(res_seq[j].residue_name)
                sub_rsa.append(res_seq[j].RSA)
                sub_phi.append(res_seq[j].Phi)
                sub_psi.append(res_seq[j].Psi)

                feature_seq[i].RSA_win.append(res_seq[j].RSA)
                feature_seq[i].Phi_win.append(res_seq[j].Phi)
                feature_seq[i].Psi_win.append(res_seq[j].Psi)
            else:
                feature_seq[i].RSA_win.append(0)
                feature_seq[i].Phi_win.append(0)
                feature_seq[i].Psi_win.append(0)
        # RSAaamean Phiaaean Psi aamean
        for j in range(20):
            Len = 0
            sumRSA = 0
            sumPhi = 0
            sumPsi = 0
            for k in range(len(sub_seq)):
                if sub_seq[k] == ResName[j]:
                    sumRSA += sub_rsa[k]
                    sumPhi += sub_phi[k]
                    sumPsi += sub_psi[k]
                    Len += 1 
            feature_seq[i].RSA_aa_mean.append(0 if sumRSA == 0 else sumRSA / Len )
            feature_seq[i].Phi_aa_mean.append(0 if sumPhi == 0 else sumPhi / Len )
            feature_seq[i].Psi_aa_mean.append(0 if sumPsi == 0 else sumPsi / Len )

        # RSAssmean Phissmean
        ss = "HEC"
        for j in range(3):
            Len = 0
            sumRSA = 0
            sumPhi = 0
            sumPsi = 0
            for k in range(len(sub_ss)):
                if sub_ss[k] == ss[j]:
                    sumRSA += sub_rsa[k]
                    sumPhi += sub_phi[k]
                    sumPsi += sub_psi[k]
                    Len += 1

            feature_seq[i].RSA_ss_mean.append(0 if sumRSA == 0 else sumRSA / Len )
            feature_seq[i].Phi_ss_mean.append(0 if sumPhi == 0 else sumPhi / Len )
            feature_seq[i].Psi_ss_mean.append(0 if sumPsi == 0 else sumPsi / Len )

        feature_seq[i].RSA_win_mean = sum(sub_rsa) / len(sub_seq)
        feature_seq[i].Phi_win_mean = sum(sub_phi) / len(sub_seq)
        feature_seq[i].Psi_win_mean = sum(sub_psi) / len(sub_seq)
        
        # PSSM
    for i in range(res_seq_len):
        sub_seq = []
        sub_col_PSSM_prob = []
        sub_col_PSSM_con = []
        sub_row_PSSM_prob = []
        sub_row_PSSM_con = []
        temp_double = []
        temp_int = []
        for j in range(20):
            sub_col_PSSM_con.append([])
            sub_col_PSSM_prob.append([])
            temp_double.append(0)
            temp_int.append(0)
        for j in range(i - half_win_len, i + half_win_len + 1):
            if j >=0 and j < res_seq_len:
                sub_seq.append(res_seq[j].residue_name)
                sub_row_PSSM_con.append(res_seq[j].PSSM_con)
                sub_row_PSSM_prob.append(res_seq[j].PSSM_prob)
            else:
                sub_seq.append('*')
                sub_row_PSSM_con.append(temp_int)
                sub_row_PSSM_prob.append(temp_double)


            for k in range(20):
                if j >=0 and j < res_seq_len and len(res_seq[j].PSSM_con) > 0 and  len(res_seq[j].PSSM_prob) > 0:
                    sub_col_PSSM_con[k].append(res_seq[j].PSSM_con[k])
                    sub_col_PSSM_prob[k].append(res_seq[j].PSSM_prob[k])
                else:
                    sub_col_PSSM_con[k].append(0)
                    sub_col_PSSM_prob[k].append(0)
        
        # window
        for j in sub_col_PSSM_con:
            for k in j:
                feature_seq[i].PSSM_con_win.append(k)
        for j in sub_col_PSSM_prob:
            for k in j:
                feature_seq[i].PSSM_prob_win.append(k)
        
        # mean, max
        for j in range(20):
            feature_seq[i].PSSM_con_col_mean.append(0 if len(sub_col_PSSM_con[j]) == 0 else int(sum(sub_col_PSSM_con[j])/window_length))
            feature_seq[i].PSSM_con_col_max.append(0 if len(sub_col_PSSM_con[j]) == 0 else max(sub_col_PSSM_con[j]))
            feature_seq[i].PSSM_con_col_min.append(0 if len(sub_col_PSSM_con[j]) == 0 else min(sub_col_PSSM_con[j]))

            feature_seq[i].PSSM_prob_col_mean.append(0 if len(sub_col_PSSM_prob[j]) == 0 else sum(sub_col_PSSM_prob[j])/window_length)
            feature_seq[i].PSSM_prob_col_max.append(0 if len(sub_col_PSSM_prob[j]) == 0 else max(sub_col_PSSM_prob[j]))

        for j in range(len(sub_seq)):
            feature_seq[i].PSSM_con_row_mean.append(0 if len(sub_row_PSSM_con[j]) == 0 else int(sum(sub_row_PSSM_con[j])/window_length))
            feature_seq[i].PSSM_con_row_max.append(0 if len(sub_row_PSSM_con[j]) == 0 else max(sub_row_PSSM_con[j]))
            feature_seq[i].PSSM_con_row_min.append(0 if len(sub_row_PSSM_con[j]) == 0 else min(sub_row_PSSM_con[j]))
            
            feature_seq[i].PSSM_prob_row_max.append(0 if len(sub_row_PSSM_prob[j]) == 0 else max(sub_row_PSSM_prob[j]))
        
        # aamean
        for j in range(20):
            sumPSSMcom = 0
            sumPSSMprob = 0
            for k in range(len(sub_seq)):
                if sub_seq[k] == ResName[j]:
                    sumPSSMcom += sub_col_PSSM_con[j][k]
                    sumPSSMprob += sub_col_PSSM_prob[j][k]
                
            
            feature_seq[i].PSSM_con_AA_mean.append(sumPSSMcom/ window_length)
            feature_seq[i].PSSM_prob_AA_mean.append(sumPSSMprob/ window_length)
        # entropy
        for j in range(len(sub_seq)):
            sum_score = 0
            templist = [p for p in sub_row_PSSM_prob[j] if not p == 0]
            for tl in templist:
                sum_score -= math.log(tl, 2) * tl
            feature_seq[i].PSSM_entropy.append(sum_score)
    # KNN
    print(res_seq[3].segment)
    KNN_progress = ProgressBar(res_seq_len, description="KNN Progress:")
    for i in range(res_seq_len):
        KNN_progress.current += 1
        # KNN_progress()
        KNN_list = []
        with open(sample_file) as sr:
            line = sr.readline()
            line = line.rstrip()
            count = 0
            while line:
                items = [s for s in line.split('\t') if not s == '']
                node = KNN_node()
                node.label = int(items[0])
                node.index = count
                # print(node.label)
                node.segment = items[3][20-half_win_len : 20 + half_win_len+1]
                # print(len(node.segment), len(res_seq[i].segment))
                node.distance = KNN_distance(res_seq[i].segment, node.segment)
                KNN_list.append(node)
                line = sr.readline()
                line = line.rstrip()
                count += 1
        KNN_list.sort(key=KNN_compare, reverse=False) 
        # print(KNN_list[0].distance)

        for j in range(1,6):
            # print(KNN_list[j].distance)
            if i == 3:
                print(KNN_list[j].segment)
            K = len(KNN_list) // int(2**j)
            # print(len([p for p in KNN_list[0:K] if p.label==1]), K)
            perc = len([p for p in KNN_list[0:K] if p.label==1]) / float(K)
            feature_seq[i].KNN.append(perc)
        # if i == 3:
        #     for j in range(1,6):
        #     print(res_seq[i].residue_name)
        #     print(KNN_list[0].label)
        #     print(KNN_list[0].distance)
        #     print(KNN_list[-1].distance)
        #     print(len([p for p in KNN_list[0:K] if p.label==1]))
        #     print(K)
        #     print(len([p for p in KNN_list[0:K] if p.label==1]) / float(K))
            # exit(1)
    KNN_progress.done()
    #HHfrewuency, HHentropy

    # get all list

    HHfrelist = []
    for i in range(res_seq_len):
        temp_list = []
        for j in range(len(ResName)):
            num = len([p for p in res_seq[i].HH_list if p == ResName[j]])
            Len = len(res_seq[i].HH_list)
            temp_list.append(num/float(Len))
        HHfrelist.append(temp_list)
        
    
    # get HHfrequency
    for i in range(res_seq_len):
        for j in range(i - half_win_len, i + half_win_len+1):
            if j >= 0 and j < res_seq_len:
                for k in range(20):
                    feature_seq[i].HH_frequency.append(HHfrelist[j][k])
            else:
                for k in range(20):
                    feature_seq[i].HH_frequency.append(0)
    # get HHentropy
    for i in range(res_seq_len):
        for j in range(i-half_win_len, i + half_win_len + 1):
            if j >= 0 and j < res_seq_len:
                sum_score = 0
                temp_list = [p for p in HHfrelist[j] if not p == 0]
                for tl in temp_list:
                    sum_score -= math.log(tl, 2) * tl
                feature_seq[i].HH_entropy.append(sum_score)
            else:
                    feature_seq[i].HH_entropy.append(0)
    
    # hsa, hsb
    for i in range(res_seq_len):
        for j in range(i - half_win_len, i + half_win_len + 1):
            if j >= 0 and j < res_seq_len:
                feature_seq[i].alpha_HSE_d.append(res_seq[j].alpha_HSE_d)
                feature_seq[i].alpha_HSE_u.append(res_seq[j].alpha_HSE_u)
                feature_seq[i].beta_HSE_d.append(res_seq[j].beta_HSE_d)
                feature_seq[i].beta_HSE_u.append(res_seq[j].beta_HSE_u)
            else:
                feature_seq[i].alpha_HSE_d.append(0)
                feature_seq[i].alpha_HSE_u.append(0)
                feature_seq[i].beta_HSE_d.append(0)
                feature_seq[i].beta_HSE_u.append(0)
    return feature_seq
def write_feature(feature, p):
    output = ""
    for j in range(len(feature)):
        output += str(p) + ":" +format_num(feature[j]) + " "
        p += 1
    return output


def write_features(res_seq, feature_seq, output_index, file_name, add_stru):
    out_file = open(file_name, 'w')
    # write_progress = ProgressBar(len(output_index), description="Writing Progress: ")
    feature_list = [
            "terminal_factor",
            "PWAA",
            "EBGW",
            "CKSAAP",
            "AAindex_mean",
            "AAindex_pair_mean",
            "factors",
            "disorder",
            "disorderCTD",
            "ss_win",
            "ss_prob",
            "ss_percent",
            "ss_seg_cnt",
            "ss_max_seg_cnt",
            "RSA_aa_mean",
            "RSA_ss_mean",
            "RSA_win_mean",
            "RSA_win",
            "Phi_aa_mean",
            "Phi_ss_mean",
            "Phi_win_mean",
            "Phi_win",
            "Psi_aa_mean",
            "Psi_ss_mean",
            "Psi_win_mean",
            "Psi_win",
            "PSSM_con_win",
            "PSSM_prob_win",
            "PSSM_con_AA_mean",
            "PSSM_prob_AA_mean",
            "PSSM_con_row_max",
            "PSSM_con_row_min",
            "PSSM_con_row_mean",
            "PSSM_con_col_max",
            "PSSM_con_col_min",
            "PSSM_con_col_mean",
            "PSSM_prob_row_max",
            "PSSM_prob_col_max",
            "PSSM_prob_col_mean",
            "PSSM_entropy",
            "HH_frequency",
            "HH_entropy",
            "alpha_HSE_d",
            "alpha_HSE_u",
            "beta_HSE_d",
            "beta_HSE_u",
            "KNN"
        ]
    with open("feature_idx_map", 'w') as idx_map:
        count = 0
        for i in output_index:

            # write_progress()
            output = ""
            if res_seq[i].label == 1:
                output += "1 "
            elif res_seq[i].label == 0:
                output += "-1 "
            else:
                continue
            p = 1
            for feat_str in feature_list:
                idx_map.write("Feature: {} starting index:{:>5d}\n".format(feat_str, p))
                feat = getattr(feature_seq[i], feat_str)
                if(type(feat) == list):
                    # list, parse every child
                        output += write_feature(feat, p)
                        p += len(feat)
                else:
                    # number
                        output += str(p) + ":" + format_num(feat) + " "
                        p += 1
            if add_stru:
                output += get_stru_feature(res_seq,count,p)
            
            out_file.write(output + '\n')
            count += 1
        # write_progress.done()
        out_file.close()
def get_stru_feature(res_seq,index, feature_id):
    fasta_file = sys.argv[1]
    aa_file = sys.argv[2]
    feature_dir = sys.argv[3]
    sample_file = sys.argv[4]
    window = int(sys.argv[5])
    residue_type = sys.argv[6]
    out_file = sys.argv[7]
    add_stru = False if len(sys.argv) == 8 else True
    sitelist = [];
    for ii in range(0,len(res_seq)):
        # print(residue_type)
        # print(res_seq[ii].residue_name)
        if res_seq[ii].residue_name == residue_type:
            sitelist.append(ii+1);
    site = [sitelist[index]]
    feature_dir += '/'
    fasta_file = fasta_file.split('/')[-1].split('.')[0]
    feature_dir = '../met_predictor/Met-Predictor_12302019/example/example_features/P0CX53_features/'
    outfilelist = get_stru_files(fasta_file,residue_type,feature_dir,feature_dir,feature_dir,feature_dir,feature_dir,feature_dir,feature_dir,feature_dir,feature_dir,site,window,feature_id);
    return outfilelist[0].strip()


def run_extract():
    fasta_file = sys.argv[1]
    aa_file = sys.argv[2]
    feature_dir = sys.argv[3]
    sample_file = sys.argv[4]
    window = int(sys.argv[5])
    residue_type = sys.argv[6]
    out_file = sys.argv[7]
    add_stru = False if len(sys.argv) == 8 else True
    print("add_stru", add_stru)
    res_seq = read_feature_files(fasta_file, aa_file, feature_dir, feature_dir, feature_dir, window)
    feature_seq = compute_feature(res_seq, window, sample_file)
    output_index = []
    for i in range(len(res_seq)):
        if res_seq[i].residue_name == residue_type:
            output_index.append(i)
    write_features(res_seq, feature_seq, output_index, out_file, add_stru)

if __name__ == "__main__":
    run_extract()
    # run_stru_extract()