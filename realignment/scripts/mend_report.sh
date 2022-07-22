in_vcf=$1
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo Mendelian Reports
bcftools +mendelian -m c -T $DIR/mend_fams.txt $in_vcf
bcftools +mendelian -m a -T $DIR/mend_fams.txt $in_vcf | python $DIR/mend_err_bases.py

