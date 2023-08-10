trbench="/Users/english/code/adotto/benchmark/GIABTR_benchmark.6.26/GIABTR.HG002.benchmark.regions.bed.gz"

function make_nums()
{
    cmp=$1
    echo "Total"
    bedtools intersect -u -f 1 -a $trbench -b $cmp | wc -l
    echo "full tier1"
    bedtools intersect -u -f 1 -a $trbench -b $cmp | grep -c Tier1
    echo "full tier2"
    bedtools intersect -u -f 1 -a $trbench -b $cmp | grep -c Tier2
    echo "Total partially contained"
    bedtools intersect -u -a $trbench -b $cmp | wc -l
    echo "partial tier1"
    bedtools intersect -u -a $trbench -b $cmp | grep -c Tier1
    echo "partial tier2"
    bedtools intersect -u -a $trbench -b $cmp | grep -c Tier2
}

echo "SV v0.6"
make_nums HG002_SVs_Tier1_v0.6_liftover38.bed
echo "v4.2.1"
make_nums HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
echo "CMRG SV"
make_nums HG002_GRCh38_CMRG_SV_v1.00.bed
echo "CMRG small"
make_nums HG002_GRCh38_CMRG_smallvar_v1.00.bed
