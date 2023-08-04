import functools
import logging
from model.UNet3D import UNet3D
from model.GeneralDiffusion import GeneralDiffusion
from model.model_utils import init_weights
def define_G(opt):
    model_opt = opt['model']
    if ('norm_groups' not in model_opt['unet']) or model_opt['unet']['norm_groups'] is None:
        model_opt['unet']['norm_groups']=32
    if opt['phase']=="train":
        schedule_opt=model_opt['beta_schedule']['train']
    else:
        schedule_opt = model_opt['beta_schedule']['val']
    #use general framework for diffusion
    encode_model = UNet3D(
        in_channel=1,
        out_channel=1,
        norm_groups=model_opt['unet']['norm_groups'],
        inner_channel=model_opt['unet']['inner_channel'],
        channel_mults=model_opt['unet']['channel_multiplier'],
        attn_res=model_opt['unet']['attn_res'],
        res_blocks=model_opt['unet']['res_blocks'],
        dropout=model_opt['unet']['dropout'],
        box_size=model_opt['diffusion']['box_size'],
        with_noise_level_emb=False
        )

    diff_model = UNet3D(
        in_channel=model_opt['unet']['in_channel'],
        out_channel=model_opt['unet']['out_channel'],
        norm_groups=model_opt['unet']['norm_groups'],
        inner_channel=model_opt['unet']['inner_channel'],
        channel_mults=model_opt['unet']['channel_multiplier'],
        attn_res=model_opt['unet']['attn_res'],
        res_blocks=model_opt['unet']['res_blocks'],
        dropout=model_opt['unet']['dropout'],
        box_size=model_opt['diffusion']['box_size']
        )
    netG = GeneralDiffusion(
        encode_model,
        diff_model,
        opt = opt,
        box_size=model_opt['diffusion']['box_size'],
        channels=model_opt['diffusion']['channels'],
        loss_type= model_opt['diffusion']['loss_type'],
        schedule_opt=schedule_opt
        )

    init_weights(netG, init_type='orthogonal')
    return netG
