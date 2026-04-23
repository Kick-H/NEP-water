#!/bin/bash
set -e
set -u

unzip 00_dlpoly_198K.zip

python dp2xyz-raw-npy-mix.py 00_dlpoly_198K nepxyz-from-npy

python get_cg_from_AA.py nepxyz-from-npy/train.xyz 00_dlpoly_198K.xyz
