#!/bin/bash
#SBATCH --job-name=NK_csv
#SBATCH --output=output.txt
#SBATCH --partition=BigMem

# Run your Python script
python3 csv_column_rewrite.py