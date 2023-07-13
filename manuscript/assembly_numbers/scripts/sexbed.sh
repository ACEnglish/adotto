#!/bin/bash
set -e

# Default values
sex=""
COV=1 # Target coverage value

# Reference files are stored beside the script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
genome=$DIR/GRCh38_1kg_mainchrs.genome_bedtools.bed
par_reg=$DIR/GRCh38.chrX.PAR.bed

# Function to display help message
function display_help() {
  echo "Usage: $0 [--mat <bed_file>] [--pat <bed_file>] [--sex <M|F>] [--out <bed_file>]"
  echo ""
  echo "Use bedtools to intersect haplotype coverage files and find single coverage"
  echo "regions of a genome."
  echo "Bed files have format 'chrom<tab>start<tab>end<tab>coverage'"
  echo 'Intermediate files are saved to $TMPDIR/$(basename $out)*'
  echo ""
  echo "Options:"
  echo "  --mat   Maternal haplotype bed file"
  echo "  --pat   Paternal haplotype bed file"
  echo "  --sex   Specify the sex as (M)ale or (F)emale"
  echo "  --out   Specify the output file"
  echo "  -h    Display this help message"
}

# Helper method to get region count, total length
function get_stats() {
    awk '{tot_cnt += 1; tot_len += $3 - $2} END {print tot_cnt "\t" tot_len }' $1
}

# Need bedtools
if ! command -v bedtools &> /dev/null
then
  echo "bedtools could not be found"
  exit 1
fi

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
  --mat)
    mat="$2"
    shift
    shift
    ;;
  --pat)
    pat="$2"
    shift
    shift
    ;;
  --sex)
    sex="$2"
    shift
    shift
    ;;
  --out)
    out="$2"
    shift
    shift
    ;;
  -h)
    display_help
    exit 0
    ;;
  *)
    echo "Unknown parameter: $1"
    exit 1
    ;;
  esac
done

# Check if required parameters are provided
if [[ -z $mat ]] || [[ -z $pat ]] || [[ -z $sex ]] || [[ -z $out ]]; then
  echo "Missing required parameters!"
  exit 1
fi

# Check if the file specified by mat exists
if [[ ! -f $mat ]]; then
  echo "The file specified by --mat does not exist: $mat"
  exit 1
fi

# Check if the file specified by pat exists
if [[ ! -f $pat ]]; then
  echo "The file specified by --pat does not exist: $mat"
  exit 1
fi

# Check if sex parameter is valid
if [[ $sex != "M" ]] && [[ $sex != "F" ]]; then
  echo "Invalid value for --sex parameter: $sex"
  exit 1
fi
echo "# Parameters:"
echo "#   mat: $mat"
echo "#   pat: $pat"
echo "#   sex: $sex"
echo "#   out: $out"

echo -e "set\tcount\tlength"
MTMP=$TMPDIR/$(basename ${out})
echo -e "mat\t$(get_stats $mat)"
echo -e "pat\t$(get_stats $pat)"

# Get 1x covered regions per-haplotype
awk '$4 == 1' $mat | cut -f1-3 > ${MTMP}_mat_cov.bed
echo -e "mat_1x\t$(get_stats ${MTMP}_mat_cov.bed)"
awk '$4 == 1' $pat | cut -f1-3 > ${MTMP}_pat_cov.bed
echo -e "pat_1x\t$(get_stats ${MTMP}_pat_cov.bed)"

# Get diploid coverage
cat ${MTMP}_mat_cov.bed ${MTMP}_pat_cov.bed \
  | bedtools sort \
  | bedtools genomecov -i - -g $genome -bga > ${MTMP}_diploid_cov.bed
echo -e "diploid\t$(get_stats ${MTMP}_diploid_cov.bed)"

# Clean out alt contigs with '_
# Non-chrX/chrY requires 1x per haplotype
awk '$1 != "chrX" && $1 != "chrY" && $4 == 2' ${MTMP}_diploid_cov.bed \
  | cut -f1-3 > ${MTMP}_cov.bed
echo -e "auto\t$(get_stats ${MTMP}_cov.bed)"

if [[ $sex == "F" ]] # Female workflow
then
  # No chrY - and 2x on chrX
  awk '$1 == "chrX" && $4 == 2' ${MTMP}_diploid_cov.bed \
    | cut -f1-3 > ${MTMP}_chrX_cov.bed
  echo -e "chrX\t$(get_stats ${MTMP}_chrX_cov.bed)"
  cat ${MTMP}_chrX_cov.bed >> ${MTMP}_cov.bed

elif [[ $sex == "M" ]]  # Male workflow
then
  # chrY requires 1x
  awk '$1 == "chrY" && $4 == 1' ${MTMP}_diploid_cov.bed \
    | cut -f1-3 > ${MTMP}_chrY_cov.bed
  echo -e "chrY\t$(get_stats ${MTMP}_chrY_cov.bed)"

  # chrX requires 2x on PAR regions, 1x everywhere else
  #    extract chrX
  grep chrX ${MTMP}_diploid_cov.bed > ${MTMP}_diploid_cov.chrX.bed
  #    subtract PAR from chrX coverage and get only 1x covered
  bedtools subtract -a ${MTMP}_diploid_cov.chrX.bed -b $par_reg \
    | awk '$4 == 1' | cut -f1-3 > ${MTMP}_chrX_np.bed
  echo -e "chrX\t$(get_stats ${MTMP}_chrX_np.bed)"

  #    compliment PAR to make non-PAR
  bedtools complement -i $par_reg -g $genome > ${MTMP}_nonpar.bed
  #    subtract non-PAR from chrX coverage and get only 2x covered
  bedtools subtract -a ${MTMP}_diploid_cov.chrX.bed -b ${MTMP}_nonpar.bed \
    | awk '$4 == 2' | cut -f1-3 > ${MTMP}_chrX_par.bed
  echo -e "PARchrX\t$(get_stats ${MTMP}_chrX_par.bed)"

  cat ${MTMP}_chrY_cov.bed >> ${MTMP}_cov.bed
  cat ${MTMP}_chrX_par.bed >> ${MTMP}_cov.bed
  cat ${MTMP}_chrX_np.bed >> ${MTMP}_cov.bed
fi

bedtools sort -i ${MTMP}_cov.bed > $out
echo -e "final\t$(get_stats $out)"
