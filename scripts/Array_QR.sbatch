#!/bin/bash
#SBATCH --job-name=Ar_QR_01
#SBATCH --time=2:00:0        
#SBATCH --partition=short
#SBATCH --hint=memory_bound
#SBATCH --mem-bind=verbose
#SBATCH --mem=50G 
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=20
#SBATCH --array=1              #e.g. 1-10%5 
#SBATCH --mail-type=FAIL
#SBATCH --output=
#SBATCH --error=

# export chr, base_pheno. set --array --job-name

IFS="$(printf '\n\t')"
module load gcc/6.2.0 R/3.6.1
module list # confirm modules used
Rscript Rserve.R

# end directories with slash /
dir= # directory containing 
snp_dir= #directory containing comma-separated lists of SNPs (fragmented for embarassingly parallel processing)
out_dir= #output directory for quantile regression data

bed=`find $dir -name "ukb_cal_"$chr"_v2.bed"` 
bim=`find $dir -iname "ukb_snp_"$chr"_*.bim"`
fam=`find $dir -iname "ukb2288_cal_"$chr"_*.fam"`

pheno= #path to phenotype 

srun --export=dir=$dir,snp_dir=$snp_dir,out_dir=$out_dir,chr=$chr,bed=$bed,bim=$bim,fam=$fam,pheno=$pheno,base_pheno=$base_pheno /home/bs263/Software/runQR/Array_rplink_02.sh

