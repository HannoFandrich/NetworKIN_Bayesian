#!/bin/bash
#SBATCH --job-name=NK_morpho_seqs
#SBATCH --output=output.txt
#SBATCH --partition=BigMem

# Run your Python script
python3 NetworKIN.py -n netphorest/netphorest -d data 9606 cured_morpho_seqs_v2.fa MS_Gaussian_updated_09032023.tsv