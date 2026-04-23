#!/bin/bash
set -eu

ln -s ../CGNEP-Datasets/00_dlpoly_198K.xyz ./train.xyz
nep > nep.out
