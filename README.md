# Renin cell identity

Here we include a step-by-step breakdown of the analyses to reproduce the data presented in "The reninness score: integrative analysis of multi-omic data to define renin cell identity".

The goal is to identify a unique epigenetic landscape that defines renin cell identity; and develop a computational tool to use that unique epigenetic landscape to identify renin-expressing cells, and quantify the renin program of unknown cell samples. 

## Overview of the experimental design
![workflow](img/concept.svg)

## Read preprocessing 
We used the PEPATAC pipeline to process the raw ATAC-seq reads, including alignment, peak-calling, and quality control.  The input files to run PEPATAC for this study are stored in the [metadata](metadata) sub-folder. For more information on how to use PEPATAC, see: http://pepatac.databio.org/

## Consensus peak set generation
We used the genomic interval machine learning (geniml) Python package to construct a consensus region set, or the “universe”, using maximum likelihood approach. For more information on how to use `genimal`, see: https://docs.bedbase.org/geniml/tutorials/create-consensus-peaks/

## Differential accessibility analysis and differential accessibile region annotation
All code used for differential accessibility analysis and differential accessibile region annotation are stored in the [src](src) sub-folder.

## Reninness score calculation and model performace evaluation
All code used for model training and reninness score calulation can be found [here](https://github.com/databio/reninness_score/tree/master)). All code used for model performace evaluation are stored in the [src](src) sub-folder.
