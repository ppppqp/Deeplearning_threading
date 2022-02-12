#!/usr/bin/env python 
import torch
import numpy as np
import torch.utils.data as data
import csv
import os 
#FEATURE = 125

class dataset(data.Dataset):
    def __init__(self, root1D, rootPair, labelRoot,  path):
        self.root1D = root1D
        self.rootPair = rootPair
        self.labelRoot = labelRoot
        self.nameCsv = csv.reader(open(path, 'r'))
        self.name = []
        for data in self.nameCsv:
            self.name.append(data[0]+'.npy')


    def __getitem__(self, index):
        ### cliff 2021-04-21
        apdbfn=os.path.basename(self.name[index])  #features file 
        # feature 1d  L* 21
        #feature1D = np.load(self.root1D + '/' + self.name[index])
        feature1D = np.load(self.root1D + '/' + apdbfn)
        # feature 2d  L* L 
        #featurePair = np.load(self.rootPair + '/' + self.name[index])
        featurePair = np.load(self.rootPair + '/' + apdbfn)
        # labels  1d  L* 1
        #label2D = np.load(self.labelRoot + '/' + self.name[index])
        label2D = np.load(self.labelRoot + '/' + apdbfn)


        feature1D = torch.FloatTensor(feature1D)

        feature2D = torch.FloatTensor(featurePair)
        feature2D = np.expand_dims(feature2D, axis=0)
        #targetTensor = torch.LongTensor(label2D)
        targetTensor = torch.FloatTensor(label2D)
        return feature1D, feature2D, targetTensor,self.name[index]

    def __len__(self):
        return len(self.name)

def get_loader(root1D, rootPair, labelRoot, path, batch_size, shuffle, num_workers):
    trainData = dataset(root1D=root1D, rootPair=rootPair,labelRoot = labelRoot, path =path)

    data_loader = torch.utils.data.DataLoader(dataset=trainData,
                                              batch_size=batch_size,
                                              shuffle=shuffle,
                                              num_workers=num_workers)
    return data_loader
