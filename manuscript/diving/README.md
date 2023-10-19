
For patho/codis,
Collect all allele deltas

I already have TRVIZ

Need a data-structure that's:

Locus   SAMPLE1_AD1  SAMPLE1_AD2 SAMPLE2_AD1    SAMPLE2_AD2

I can then make a boxplot with:

X - allele delta (-N to N)
Y - Each locus
Item_A - Boxplot of AD
Point_B - HG002 AD Maternal
Point_C - HG002 AD Paternal

Hue - Codis/Pathogenic/Phenotypic (or Tier?)
If Tier2 are <20% of sites (its not, I think) that can be hue.
Then we put up blocks along Y-right-axis of Codi/Path/Phen
