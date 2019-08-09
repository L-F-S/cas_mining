#!/bin/bash
# Fri Jun 21 10:37:37 CEST 2019
# Made by L-F-S
# At the University Of Trento, Italy
# cycles through all datasets and launches any python script giving it a dataset name

#usage:

# beware! metti il giusto python script dentro
# chmod +x thisfilename
# sh thisfilename


for dataset in AsnicarF_2017 BackhedF_2015 BengtssonPalmeJ_2015 BritoIL_2016 CM_cf CM_madagascar CM_periimplantitis Castro-NallarE_2015 ChengpingW_2017 ChngKR_2016 CosteaPI_2017 LawrenceA_2015 FengQ_2015 CM_caritro GeversD_2014 HMP_2012 HanniganGD_2017 HeQ_2017 IjazUZ_2017 KarlssonFH_2013 KosticAD_2015 LeChatelierE_2013 LiJ_2014 LiJ_2017 LiSS_2016 LiuW_2016 LomanNJ_2013 LoombaR_2017 LouisS_2016 NielsenHB_2014 Obregon-TitoAJ_2015 OhJ_2014 OlmMR_2017 QinJ_2012 QinN_2014 RampelliS_2015 RaymondF_2016 SchirmerM_2016 SmitsSA_2017 VatanenT_2016 VincentC_2016 VogtmannE_2016 WenC_2017 XieH_2016 YuJ_2015 ZeeviD_2015_A ZeeviD_2015_B ZellerG_2014; do
echo $dataset
time python compare_minced_prokka.py $dataset &
done
