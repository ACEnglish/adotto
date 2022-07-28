mkdir -p jobs
mkdir -p logs

ref=$1

grep -vw 'rel_path' assemblies/relative_path_metadata.txt | while read sample project haplotag rel_path
do
    full_path=$(realpath assemblies/${rel_path})
    name=${project}_${sample}_${haplotag}
    out_dir=initial_alignment/${name}
    echo "#!/bin/bash" > jobs/aln_${name}.sh
    echo "bash scripts/map_haplo.sh ${full_path} $ref ${sample}.${haplotag} ${out_dir}" >> jobs/aln_${name}.sh
done

