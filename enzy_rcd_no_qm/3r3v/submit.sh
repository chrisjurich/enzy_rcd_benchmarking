#!/bin/bash
#SBATCH -o "3r3v_no_qm.txt"
#SBATCH --job-name=3r3v_no_qm
#SBATCH --account=yang_lab_csb
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=10G
#SBATCH --time=12:00:00
source $HOME/setup-accre.sh


cd /data/yang_lab/jurichc/enzy_rcd_benchmarking/enzy_rcd_no_qm/3r3v && python main.py
