#!/bin/bash
#SBATCH --job-name=NK_test
#SBATCH --output=output.txt


# Run your Python script
python3 NetworKIN_CODE_v3.0/NetworKIN.py -n netphorest/netphorest -d NetworKIN_CODE_v3.0/data 9606 test1.fas