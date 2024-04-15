#!/bin/bash
WORK_DIRECTORY=$1

REALPATH_WK_DIR=`realpath $WORK_DIRECTORY`

full_path_pip=`realpath $0`
path_pip=`dirname $full_path_pip`

echo "WORK_DIRECTORY : $REALPATH_WK_DIR"

python3 $path_pip/analyse_mapping.py "${REALPATH_WK_DIR}/OUTSIDER/TrEMOLO_SV_TE/INS/SV_INS_CLUST.bln" \
    "${REALPATH_WK_DIR}/1-UTILS/TE_SIZE.tsv" \
    "${REALPATH_WK_DIR}/OUTSIDER/TrEMOLO_SV_TE/INS/SV_SIZE.tsv" \
    "${REALPATH_WK_DIR}/data.json"

rm -f $path_pip/../back-end/data/data_GEN.json
rm -f $path_pip/../back-end/data/data_GEN_*
rm -f $path_pip/../back-end/index/*
rm -f $path_pip/../back-end/DEPTH/*
rm -f $path_pip/../back-end/TE_INFOS/*

ln -sr "${REALPATH_WK_DIR}/data.json" $path_pip/../back-end/data/data_GEN.json

echo "DONE !"
