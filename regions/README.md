In progress

We need to define what regions of the genome contain tandem-repeats.

To accomplish this, we collected tandem repeat region bed files from various sources.

Formatting scripts to turn them all into bed files are located in
`scripts/formatters/`

Next, the bed files are merged

`stats` somewhere

Then I added a 75bp slop

Then I removed the variants >= .. 50kb

Then I removed variants within Nbp of reference gaps.

Then I ran TRF on the reference sequence of the remaining regions

Then translate that TRF output back to genomic coordinates and format to a bed/gz/tbi

Then I checked the input source beds against this set of regions to ensure that our new set
at least somewhat is representative of the input beds. For example, source ABC has a region
on chr:pos-end with motif GGG. Does this final produced bed have that same (or a similar) 
repeat description?

`stats` of the remaining regions
