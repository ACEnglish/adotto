# 0 - assume input is
#       reference
#       Output directory from map_haplo.sh
# 1 - Will create a new reference using fix-ref first (since we're using whatever)
# 2 - Will then just run map_haplo again

ref=$(realpath $1)
in_dir=$(realpath $2)
vcf=$in_dir/aln.covanno.vcf.gz
out_dir=${in_dir}/realign/

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

mkdir -p $out_dir

bcftools view -i "FILTER == 'PASS'" $vcf \
    | bcftools +fill-from-fasta /dev/stdin -- -c REF -f ${ref} \
    | bgzip > ${out_dir}/variants.vcf.gz
tabix ${out_dir}/variants.vcf.gz

# Real pipeline should have -i here
bcftools consensus -H1 -f $ref ${out_dir}/variants.vcf.gz > ${out_dir}/haplotype.fa

bash $DIR/map_haplo.sh ${out_dir}/haplotype.fa $ref sample ${out_dir}
