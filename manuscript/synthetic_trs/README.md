I want to make two experiments exploring the impact of Truvari's improvements for TRs.

Experiment 1: Unroll

For a number of tandem repeats, create two alternate allele that are:
a) within the TR's boundaries
b) of the same type (expansion/contraction)
c) of the same size.
d) the 'base' variant will be at the 5' most position of the TR
e) the 'comp' variant will be placed at some number of positions downstream
f) we'll collect:
- size of the variant
- distance between the two
- refcontext similarity
- unroll similarity

we can repeat this over the same TR regions any number of times with the assumption that the
randomization will make them different each iteration.

I don't have to mess around with trying to make VCFs. I can do this all in-memory with the Truvari API


Experiment 2: split-variants

Same 'random variant' maker used above (without the shifting)
But, for the second one, we want to make it 'split'.

For this, I think I actually do want VCFs? This is weird because I'm making the same haplotype, just with different
representations, and then i'm recreating that haplotype to show that the output VCF has one joint representation.

Okay, so yes I want VCFs.
And I'll show their counts and base/comparison similarity. Run it through phab. Show that they have identical
counts/100% similarity.


