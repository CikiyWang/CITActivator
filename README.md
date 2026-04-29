# CITActivator
Cap-independent translation activators

## circos_heatmap:
Translation Efficiency, Luciferase Activity, Splicing Efficiency, and CircRNA Amount from raw experimental data of 110 RBPs.
Plot using Circos R package.

## SPLSDA_classification:
R scripts for SPLSDA model training and prediction, including cross validation and feature selection in each fold.
dir: ./SPLSDA_classification
usage: Rscripts splsda_train_and_evaluation.R
input matrix: 

## RBP_eCLIP:
python scripts for peak coverage calculation in each single nucleotide resolution from eCLIP peak data.

## IRES_overlap:
python scripts for overlap currently reported IRES postion with RBP eCLIP peaks.
dir: ./IRES_overlap
usage: python peak_ires_overlap.py {RBP} region(150)
output: {RBP}_peak_ires_overlap.txt

