#!/usr/bin/env python
import torch
from dataloader import *
from model import *
import argparse
import os
import numpy as np  #from scipy.stats import pearsonr 
import pandas as pd 
from matplotlib import pyplot as plt
#torch.cuda.set_device(0)
#BATCHSIZE = 1
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
# device = torch.device('cpu')
parser = argparse.ArgumentParser()
parser.add_argument('--sequential_features', type=str, default='./features/d1', help='The full path of sequential features of test data.')
parser.add_argument('--pairwise_features', type=str, default='./features/d2', help='The full path of pairwise features of test data.')
parser.add_argument('--target_output', type=str, default='./data/train_label', help='The full path of the target outputs.')
parser.add_argument('--test_list', type=str, default='eg.list', help='The full path of the test list.')
parser.add_argument('--models_path', type=str, default='./models', help='path list of models.')
config = parser.parse_args()


def cov(a,b):
    n = a.shape[1]
    return np.sum((a - np.mean(a)) * (b - np.mean(b))) / (n-1)


def PCC(a,b):
    #Pearson's correlation coefficient = covariance(X, Y) / (stdv(X) * stdv(Y))
    return cov(a,b) / (np.std(a) * np.std(b))

def MAE(a, b):
    n = a.shape[1]
    return np.sum(np.absolute(a - b))/n


if __name__ == '__main__':
    myname = config.models_path
    pccs = []
    maes = []
    length = []
    nameList=os.listdir(myname)#models' names
    #fwn0="sprofpred"
    #print("#modelname overall_pcc average_pcc overall_mae average_mae")
    for modelname in nameList:
        net = SPROF(ResidualBlock).eval()
        net = net.to(device)
        net.load_state_dict(torch.load(f'{config.models_path}/{modelname}'))
        data_loader = get_loader(config.sequential_features, config.pairwise_features, config.target_output,
                                 config.test_list, batch_size=1, shuffle=False, num_workers=0)
        with torch.no_grad():
            for i, (test1d, test2d, targetdata, name) in enumerate(data_loader):
                print(i)
                print("pdb: %s" %name)
                try:
                    test1d = test1d.to(device)
                    test2d = test2d.to(device)
                    targetdata=targetdata.float()  #cliff 2020-05-31 
                    targetdata = targetdata.to(device)
                    output = net(test1d, test2d)  # model output
                    output = output.to(device)
                    pcc = PCC(targetdata.numpy(), output[0].numpy())
                    mae = MAE(targetdata.numpy(), output[0].numpy())
                    print(pcc)
                    print(mae)
                    length.append(targetdata.numpy().shape[1])
                    pccs.append(pcc)
                    maes.append(mae)

                except ValueError:
                    with open("/nfs/amino-home/panqp/protein_stru/repo/EAGERER/unsupported", 'a') as f:
                        f.write("%s\n" %name)
    print(np.mean(np.array(pccs)))
    print(np.mean(np.array(maes)))
    plt.plot(np.log10(np.array(length)), np.array(pccs))
    plt.savefig('pcc.png')
    print("Predition finshed")