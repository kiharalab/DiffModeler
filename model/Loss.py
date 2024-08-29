import torch
import torch.nn as nn
import torch.nn.functional as F

class SoftDiceLoss(nn.Module):
    #do not use softmax for input, instead, using the relative weak sigmoid to train
    def __init__(self,weight=None):
        super(SoftDiceLoss, self).__init__()
        self.eps = 1e-6
        self.weight = weight
    def forward(self,input,target):
        if not input.shape == target.shape:
            raise ValueError("input and target shapes must be the same. Got: Input:{} Target:{}"
                             .format(input.shape, target.shape))
        if not input.device == target.device:
            raise ValueError(
                "input and target must be in the same device. Got: {}" .format(
                    input.device, target.device))
        input_soft = F.sigmoid(input)
        #if self.weight is None:
        self.weight = [1 for k in range(input_soft.shape[1])]
        assert input_soft.shape[1]==len(self.weight)
        dims = [k for k in range(0,len(input.shape))]
        loss = 0
        for j in range(len(self.weight)):
            intersection = torch.sum(input_soft[:,j] * target[:,j], dims)
            cardinality = torch.sum(input_soft[:,j] + target[:,j], dims)#take square is wrong here
            dice_score = 2. * intersection / (cardinality + self.eps)
            loss += self.weight[j]*torch.mean(1. - dice_score)
        loss /= len(self.weight)
        return loss

class BatchDiceLoss(nn.Module):
    #do not use softmax for input, instead, using the relative weak sigmoid to train
    def __init__(self,weight=None):
        super(BatchDiceLoss, self).__init__()
        self.eps = 1e-6
        self.weight = weight
    def forward(self,input,target):
        if not input.shape == target.shape:
            raise ValueError("input and target shapes must be the same. Got: Input:{} Target:{}"
                             .format(input.shape, target.shape))
        if not input.device == target.device:
            raise ValueError(
                "input and target must be in the same device. Got: {}" .format(
                    input.device, target.device))
        input_soft = F.sigmoid(input)
        #if self.weight is None:
        self.weight = [1 for k in range(input_soft.shape[1])]
        assert input_soft.shape[1]==len(self.weight)
        dims = [k for k in range(1,len(input.shape))]
        loss = 0
        for j in range(len(self.weight)):
            intersection = torch.sum(input_soft[:,j] * target[:,j], dims)
            cardinality = torch.sum(input_soft[:,j] + target[:,j], dims)#take square is wrong here
            dice_score = 2. * intersection / (cardinality + self.eps)
            loss += self.weight[j]*torch.mean(1. - dice_score)
        loss /= len(self.weight)
        loss = torch.mean(loss)
        return loss
