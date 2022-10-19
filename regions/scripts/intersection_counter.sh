echo v0.1
echo pbsv
bedtools intersect -u -a data/pbsv/merged.bed.gz -b delme_v0.1/data/tr_regions.bed.gz | wc -l
echo usc
bedtools intersect -u -a data/usc/merged.bed.gz -b delme_v0.1/data/tr_regions.bed.gz | wc -l
echo trgt
bedtools intersect -u -a data/trgt/merged.bed.gz -b delme_v0.1/data/tr_regions.bed.gz | wc -l

echo v0.2
echo pbsv
bedtools intersect -u -a data/pbsv/merged.bed.gz -b adotto_TRannotations_v0.2.bed.gz | wc -l
echo usc
bedtools intersect -u -a data/usc/merged.bed.gz -b adotto_TRannotations_v0.2.bed.gz | wc -l
echo trgt
bedtools intersect -u -a data/trgt/merged.bed.gz -b adotto_TRannotations_v0.2.bed.gz | wc -l

echo v0.3
echo pbsv
bedtools intersect -u -a data/pbsv/merged.bed.gz -b adotto_TRannotations_v0.3.bed.gz | wc -l
echo usc
bedtools intersect -u -a data/usc/merged.bed.gz -b adotto_TRannotations_v0.3.bed.gz | wc -l
echo trgt
bedtools intersect -u -a data/trgt/merged.bed.gz -b adotto_TRannotations_v0.3.bed.gz | wc -l
