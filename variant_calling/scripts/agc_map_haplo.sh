#!/bin/bash

# Wrapper around agc files to perform extraction per-sample and then run map_haplo.sh
agc_input=$1
sample_name=$2
ref=$3
out_prefix=$4

# hard coding ... boo
agc=/users/u233287/scratch/misc_software/agc-2.1_x64-linux/agc
threads=8

project=hprc

sample=$(echo ${sample_name} | cut -f1 -d\.)
haplotag=H$(echo ${sample_name} | cut -f2 -d\.)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

name=${project}_${sample}_${haplotag}
out_dir=${out_prefix}/${name}


fa_path='$TMPDIR'/${name}.assembly.fa
echo "${agc} getset -t 8 $agc_input ${sample_name} >  $fa_path"
echo bash ${DIR}/map_haplo.sh ${fa_path} $ref ${sample}.${haplotag} ${out_dir}
