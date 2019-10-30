#!/bin/sh

#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --job-name=scMatch_run
#SBATCH --time=00-24:00:00
#SBATCH --mem=32GB

module load pandas
module load scipy
module load scmatch

scMatch.py \
--refType mouse \
--testType mouse \
--refDS PATH_TO_FANTO5_DATABASE/FANTOM5 \
--dFormat 10x \
--testDS PATH_TO_CELLRANGER_OUTPUT/outs/filtered_feature_bc_matrix


