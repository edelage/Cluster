#!/bin/sh

## Parametres d'entree
input=$1
cluster=$2
heatmap=$3
outputPng=$4
outputSvg=$5

## Constantes
CLUSTERDIR= '/home/galaxy/galaxy-dist/local_tools/Cluster'
CIRCOSBIN= '/home/galaxy/galaxy-dist/tools/circos/circos-0.66/bin/circos'


## Parsing du fichier d'entree
python ${CLUSTERDIR}/cluster.py $cluster $input $heatmap

## Lancement de circos
perl $CIRCOSBIN -conf ${CLUSTERDIR}/conf/circos.conf  -outputdir $(dirname $outputPng) -outputfile $(basename $outputPng)  


## Renommage du fichier de sortie pour Galaxy
mv ${outputPng}.png $outputPng
mv ${outputPng}.svg $outputSvg


## Suppression des fichiers temporaires créés par le script python
##rm /home/galaxy/galaxy-dist/local_tools/Cluster/data/tmp/*

