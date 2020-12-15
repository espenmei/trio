## Trio-GCTA

This repository contains simulated data and code for fitting the model described in ref.

The data consists of 1000 simulated parent-offspring trios.
It is arranged so that the first 1000 individuals are mothers, the next 1000 are fathers and the last 1000 are children.
Within each group, the data is sorted according to trios.
For example, individual 1, 1001 and 2001 is the mother, father and offspring of one trio.
The data is stored in typical plink format with *.fam*, *.bed* and *.bim* files.

The file *trio.R* uses **OpenMx** to fit the model in **R**.

The file *trio.jl* uses **VCModels.jl** to fit the model in **Julia**.
This code is currently experimental and not recommended for use.

