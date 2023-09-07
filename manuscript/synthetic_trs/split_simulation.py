"""
Similar to seq similarity, we're going to make some random TRs.

But, for the rep2, we're going to split them up.
For example, +2 copies, we put them in two different positions.
(Do I we want to allow partials splits?)

Then we want to run them through mafft and see if variants are starting/ending
in the same position, having the same counts.

We can again do all this in-memory.

Ho - harmonized variants are equally comparable us unharmonized vriants
Ha - harmonization improves comprability
"""
