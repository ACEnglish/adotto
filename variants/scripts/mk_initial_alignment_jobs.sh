mkdir -p initial_alignments
mkdir -p jobs
mkdir -p logs

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

ref=$1

grep -vw 'rel_path' assemblies/relative_path_metadata.txt | while read sample project haplotag rel_path
do
    full_path=$(realpath assemblies/${rel_path})
    name=${project}_${sample}_${haplotag}
    out_dir=initial_alignments/${name}
    echo "#!/bin/bash" > jobs/aln_${name}.sh
    echo "bash $DIR/map_haplo.sh ${full_path} $ref ${sample}.${haplotag} ${out_dir}" >> jobs/aln_${name}.sh
done

