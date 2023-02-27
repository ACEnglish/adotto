# use docker to setup R easily
# docker run --rm -p 8787:8787 --cpus=4 -e PASSWORD=yourpasswordhere rocker/rstudio

# Update container with required libraries
# docker exec -it <container_name> /bin/bash
# apt-get -qq update && apt-get install -yq libxml2 liblzma-dev libbz2-dev libz-dev

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("regioneR")

library(regioneR)

genome <- toGRanges("/data/GRCh38_1kg_mainchrs.fa.fai")
my.regions <- toGRanges("/data/grch38.simplerepeat_merge_max50.bed", genome=genome)
my.mask <- toGRanges("/data/grch38.exclude_regions.bed", genome=genome)

all.gene <- toGRanges("/data/grch38.proteinGenes_merged.bed", genome=genome)

ptO_genes <- overlapPermTest(A=my.regions, B=all.gene, ntimes=100, genome=genome, per.chromosome=TRUE,
                        mask=my.mask, verbose=TRUE, force.parallel=TRUE)
plot(ptO_genes)

all.gene <- toGRanges("/data/grch38.proteinGenes_merged_promoters.bed", genome=genome)
ptO_promoter <- overlapPermTest(A=my.regions, B=all.gene, ntimes=100, genome=genome, per.chromosome=TRUE,
                                mask=my.mask, verbose=TRUE, force.parallel=TRUE)
plot(ptO_promoter)

