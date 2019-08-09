#!/bin/bash
# Wed Jul 3 17:53:10 CEST 2019
# Made by L-F-S
# At the University Of Trento, Italy
# cycles through all datasets and launches any python script giving it a dataset name

#usage:

# beware! metti il giusto python script dentro
# chmod +x thisfilename
# sh thisfilename
#TODO occhio a zeevid, fatyto un ciclo solo per alcuni che non venivano, noi n ci sono piu tutti!


for dataset in BengtssonPalmeJ_2015 Castro-NallarE_2015 LawrenceA_2015 CM_caritro IjazUZ_2017 LiJ_2014 OhJ_2014 RaymondF_2016 SchirmerM_2016  VatanenT_2016 VincentC_2016  WenC_2017 XieH_2016 YuJ_2015 ZeeviD_2015; do
echo $dataset
python tabellazza.py $dataset &
done
