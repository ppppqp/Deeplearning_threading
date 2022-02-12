import torch
from dataloader import *
from model import *
import argparse
import os
import numpy as np  #from scipy.stats import pearsonr 
import pandas as pd 
import torch.optim as optim
from torch.autograd import Variable
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


if __name__ == '__main__':
    myname = config.models_path
    nameList=["met_pred"]#models' names

    for modelname in nameList:
        net = SPROF(ResidualBlock).eval()
        net = net.to(device)
        lr = 0.0005
        epoach_num = 30
        criterion = nn.MSELoss()
        optimizer = optim.Adam(net.parameters(), lr=lr)

        # net.load_state_dict(torch.load(f'{config.models_path}/{modelname}'))
        data_loader = get_loader(config.sequential_features, config.pairwise_features, config.target_output,
                                 config.test_list, batch_size=1, shuffle=False, num_workers=0)
        for step in range(epoach_num):
            for i, (test1d, test2d, targetdata, name) in enumerate(data_loader):
                print("pdb: %s" %name)
                print("step:", step)
                try:
                    optimizer.zero_grad()
                    test1d = test1d.to(device)
                    test2d = test2d.to(device)
                    targetdata=targetdata.float()  #cliff 2020-05-31 
                    targetdata = targetdata.to(device)
                    output = net(test1d, test2d)  # model output
                    output = output.to(device)
                    # alb=targetdata[0].cpu().numpy()
                    # print(targetdata.shape)
                    # print(targetdata)
                    # apd=output[0].cpu().numpy()[0]
                    print(output[0].shape, targetdata.shape)  
                    if output[0].shape != targetdata.shape:
                        raise Exception('label shape not valid')                 
                    loss = criterion(output[0], targetdata)
                    print(loss.data)
                    # loss = Variable(loss, requires_grad = True)
                    loss.backward()
                    optimizer.step()
                except Exception as e:
                    msg = e.args
                    if msg == 'label shape not valid':
                        with open("/nfs/amino-home/panqp/protein_stru/repo/EAGERER/label_invalid", 'a') as f:
                            f.write("%s\n" %name)
            if step % 10 == 9:
                torch.save(net.state_dict(), "./models/met_pred_"+str(step)+".pt")

   
    print("Train finshed")
