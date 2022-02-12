#!/usr/bin/env python
import torch
from dataloader import *
from model import *
import argparse
import os
import numpy as np  #from scipy.stats import pearsonr 
import pandas as pd 
#torch.cuda.set_device(0)
#BATCHSIZE = 1
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
# device = torch.device('cpu')
parser = argparse.ArgumentParser()
parser.add_argument('--sequential_features', type=str, default='./features/d1', help='The full path of sequential features of test data.')
parser.add_argument('--pairwise_features', type=str, default='./features/d2', help='The full path of pairwise features of test data.')
parser.add_argument('--target_output', type=str, default='./features/lbs', help='The full path of the target outputs.')
parser.add_argument('--test_list', type=str, default='test_list', help='The full path of the test list.')
parser.add_argument('--models_path', type=str, default='./models', help='path list of models.')
config = parser.parse_args()

if __name__ == '__main__':
    myname = config.models_path
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
                    alb=targetdata[0].cpu().numpy()
                    apd=output[0].cpu().numpy()[0]
                    alb_pd=np.column_stack((alb,apd))
                    print(targetdata.shape)
                    df=pd.DataFrame(alb_pd,columns=["No","RSA"])
                    df["No"]=df["No"].astype(np.int64) #,dtype={"No":np.int64,"RSA":np.float64})
                    outfn="/nfs/amino-home/panqp/protein_stru/repo/EAGERER/outputs/%s.out" %os.path.splitext(name[0])[0].split('/')[-1]
                    df.to_csv(outfn,index=False,header=True,float_format="%10.6f")
                except ValueError:
                    with open("/nfs/amino-home/panqp/protein_stru/repo/EAGERER/unsupported", 'a') as f:
                        f.write("%s\n" %name)

    print("Predition finshed")