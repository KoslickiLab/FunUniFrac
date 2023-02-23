#!/bin/bash
conda env create -f environment.yml --force
conda deactivate
conda activate fununifrac 

####################
# failures.
####################
# - Pip failures
# numpy => Removed the environment if exists. Or retried setup.sh again.
