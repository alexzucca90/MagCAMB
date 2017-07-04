# MagCAMB
Version 0.0.

Patch for CAMB (http://camb.info) that computes CMB anisotropies sourced by Primordial Magnetic Fields (PMF)

- integrated with CosmoMC (MagCosmoMC). See https://github.com/alexzucca90/MagCosmoMC. 
- based on CAMB 2015
- part of the code is based on the "Magnetic CAMB" version by Richard Shaw, provided by private communication.

INSTALL: 
- Copy-paste the files in your CAMB directory 
- make

RUN:
- modify params_mag.ini according to what you want to calculate 
- ./camb params_mag.ini

If you use this code for a scientific work, please cite our paper, arXiv:1611.00757, and the original CAMB paper, http://arxiv.org/abs/astro-ph/9911177.

Alex Zucca


