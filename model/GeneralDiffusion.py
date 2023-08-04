import math
import torch
from torch import device, nn, einsum
import torch.nn.functional as F
import torch.nn as nn
from inspect import isfunction
from functools import partial
import numpy as np
from tqdm import tqdm
from model.Loss import BatchDiceLoss

def gamma(t, ns=0.0002, ds=0.00025):
    # A scheduling function based on cosine function.
    return np.cos( ((t + ns) / (1 + ds)) * np.pi / 2)**2

def exists(x):
    return x is not None


def default(val, d):
    if exists(val):
        return val
    return d() if isfunction(d) else d

class GeneralDiffusion(nn.Module):
    def __init__(
        self,
        encode_fn,
        denoise_fn,
        box_size,
        opt=None,
        channels=1,
        loss_type="both",
        schedule_opt=None
    ):
        super().__init__()
        self.channels = channels
        self.encode_fn = encode_fn
        self.box_size = box_size
        self.denoise_fn = denoise_fn
        self.opt = opt
        self.loss_type = loss_type
        self.schedule_opt = schedule_opt
        self.num_timesteps = self.schedule_opt['n_timestep']
        self.clip_flag = self.schedule_opt['infer_clip']



    def set_loss(self):
        self.loss_func =BatchDiceLoss()

    def q_sample(self, x_start, gamma_t, noise=None):
        """
        :param x_start: density x
        :param continuous_sqrt_alpha_cumprod: [gamma_t-1, gamma_t]
        :param noise: random noise
        :return:
        """
        noise = default(noise, lambda: torch.randn_like(x_start))

        # random gama from sqrt(gamma_[t-1]) sqrt(gamma_[t])

        return (
            torch.sqrt(gamma_t) * x_start +
            torch.sqrt(1 -gamma_t) * noise
        )

    def p_losses(self, x_in, noise=None):
        x_start = x_in['backbone']
        [b, c, h, w, l ] = x_start.shape
        batch_size = b

        t = np.random.uniform()
        gamma_t = gamma(t)
        # t = torch.FloatTensor(
        #     [t]).repeat(batch_size, 1).to(x_start.device)
        gamma_t = torch.FloatTensor([gamma_t]).to(x_start.device)
        gamma_t_batch = gamma_t.repeat(batch_size, 1)

        print("sampled t:",t," gamma t:",gamma_t)

        noise = default(noise, lambda: torch.randn_like(x_start))
        input_decode = x_start
        x_noisy = self.q_sample(input_decode,gamma_t,noise=noise)
        density_condition = x_in['density']

        density_condition = self.encode_fn(density_condition)
        density_condition_norm = torch.sigmoid(density_condition)

        x_recon = self.denoise_fn(
                torch.cat([density_condition_norm, x_in['density'], x_noisy], dim=1), gamma_t_batch)


        loss = self.loss_func(x_recon,x_start)


        return loss,x_recon,x_start



    def forward(self, x, *args, **kwargs):
        return self.p_losses(x, *args, **kwargs)
    @torch.no_grad()
    def p_sample(self, x, t, clip_denoised=True, condition_x=None,td=1.0):
        batch_size = x.size(0)

        t_now = 1-t/self.num_timesteps
        t_next = max(1 - (t + 1 + td) / self.num_timesteps, 0)
        gamma_t = gamma(t_now)
        gamma_t = torch.FloatTensor([gamma_t]).to(x.device)
        gamma_t_batch = gamma_t.repeat(batch_size, 1)
        x_predict = self.denoise_fn(
                torch.cat(condition_x+[x], dim=1), gamma_t_batch)
        #through experiments,sigmoid should be removed here
        #apply condition normalization
        x_predict = torch.sigmoid(x_predict)
        x_predict = x_predict*2 -1 #prepare for clip
        print("infer: t:",t,"t now:",t_now,"gamma t", gamma_t)
        gamma_t_next = gamma(t_next)
        gamma_t_next = torch.FloatTensor([gamma_t_next]).to(x.device)
        gamma_t_next_batch =  gamma_t_next.repeat(batch_size, 1)
        if clip_denoised:
            x_predict = x_predict.clamp_(-1., 1.)#scale -0.1,0.1 is suggested in the general diffusion paper
        eps = 1/torch.sqrt(1-gamma_t) *(x-x_predict*torch.sqrt(gamma_t))
        x_next = torch.sqrt(gamma_t_next)*x_predict+torch.sqrt(1-gamma_t_next)*eps
        return x_next


    @torch.no_grad()
    def p_sample_loop(self, x_in, continous=False,pred_data=None):
        device = x_in.device
        sample_inter = (1 | (self.num_timesteps//10))
        x = x_in
        input_density = x

        x = self.encode_fn(pred_data)
        x_norm = torch.sigmoid(x)

        #shape = x.shape
        shape = (x.size(0), self.channels, self.box_size, self.box_size,self.box_size)
        img = torch.randn(shape, device=device)
        ret_img = img#to align shape well
        for i in tqdm(range(0, self.num_timesteps), desc='sampling loop time step', total=self.num_timesteps):
            condition_x = [x_norm,input_density]

            img = self.p_sample(img, i, condition_x=condition_x,clip_denoised=self.clip_flag)
            if i % sample_inter == 0:
                cat_img = img
                ret_img = torch.cat([ret_img, cat_img], dim=0)
        if continous:
            return ret_img
        else:
            return ret_img[-x_in.size(0):]
    @torch.no_grad()
    def super_resolution(self, x_in, continous=False,pred_data=None):
        return self.p_sample_loop(x_in, continous,pred_data=pred_data)
