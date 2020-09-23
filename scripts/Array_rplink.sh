#!/bin/bash

date
echo "Node:${SLURMD_NODENAME} Core:${SLURM_PROCID} Array:${SLURM_ARRAY_JOB_ID} Task:${SLURM_ARRAY_TASK_ID}"
export PATH=$PATH':Software/plink_linux_x86_64_20200616'
cd $dir

sub=`expr ${SLURM_ARRAY_TASK_ID}`
echo "sub: "$sub
snps=$snp_dir'snp_'$sub'.txt'
echo "snpfile: "$snps
snp_list=`cat $snps`

function='QR_plugin.R'
out=$out_dir$base_pheno"_"$chr"_QR_n"$sub
covar_path='ukb.cov'
covar_name='sex,age,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10'

plink --threads 1 --memory 49000 --bed $bed --bim $bim --fam $fam --pheno $pheno --snps $snp_list --d : --R $function --out $out --covar $covar_path --covar-name $covar_name

date
echo Done: $out $function 
