#!/bin/bash
# Thu Jul 11 10:10:30 CEST 2019
# Made by L-F-S
# At the University Of Trento, Italy

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, '/home/lorenzo.signorini/utils/')
import filename_discrepancies


dataset="ZellerG_2014" #sys.argv[1] #TODO switch test dataset
outdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/Nstep" #TODO add nstepdir
tabledir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/4preprocesseddataviz" #TODO add nstepdir
CRISPRdir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/1crisprsearch/out"
annodir="/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/2casanno/crisprcasanno"
