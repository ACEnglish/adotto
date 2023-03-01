Found this paper:
https://www.nature.com/articles/s41598-021-82050-5

They have 
  6101 41598_2021_82050_MOESM2_ESM.txt

Interesting CGG repeats on hg19. After liftover to grch38 we have
  6093 aut_cgg.grch38.bed

After turning the start positions into ranges (repeat length * median length)
We intersect with adotto_TRregions v1 and see:

$ bedtools intersect -c -a aut_cgg.grch38.bed -b ../adotto_TRannotations_v1.0.bed.gz  | cut -f6 | sort | uniq -c
   2043 0
   3939 1
    100 2
      6 3
      2 4
      3 5

Only 66.4% are in our catalog - Rebuild the results with these.


