#!/bin/bash
#SBATCH -o "8brb"
#SBATCH --job-name=8brb
#SBATCH --account=yang_lab_csb
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=10G
#SBATCH --time=12:00:00
source $HOME/setup-accre.sh

#TODO(CJ)
cd /data/yang_lab/jurichc/enzy_rcd_benchmarking/enzy_rcd_no_qm/8brb && python main.py
