mkdir -p data/

#gene_gtf_url=http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_35/gencode.v35.annotation.gtf.gz
#gene_gtf=data/gencode.annotation.gtf.gz
#wget $gene_gtf_url -O $gene_gtf

# protein coding transcripts
#gunzip -c $gene_gtf \
#    | grep 'transcript_type "protein_coding"' \
#    | awk '($3=="transcript") {printf("%s\t%s\t%s\n",$1,int($4)-1,$5);}' \
#    | sort -T . -t $'\t' -k1,1 -k2,2n \
#    | bedtools merge > data/grch38.protein_coding_transcript.bed

#gunzip -c $gene_gtf \
#    | grep 'transcript_type "protein_coding"' \
#    | awk '($3!="transcript") {printf("%s\t%s\t%s\n",$1,int($4)-1,$5);}' \
#    | sort -T . -t $'\t' -k1,1 -k2,2n \
#    | bedtools merge > data/grch38.protein_coding_nonintrons.bed

# bedtools subtract -a data/grch38.protein_coding_transcript.bed \
#                   -b data/grch38.protein_coding_nonintrons.bed \
#                   > data/grch38.protein_coding_introns.bed

# EPD promoters
python -c """
import urllib.request
import json
import pandas as pd

url = 'https://api.genome.ucsc.edu/getData/track?genome=hg38;track=epdNewPromoter'
with urllib.request.urlopen(url) as response:
    html = response.read()
    data = json.loads(html)['epdNewPromoter']
    for i in data:
        print(i['chrom'], i['chromStart'], i['chromEnd'], sep='\t')
""" > data/grch38.epd_promoters.bed

# CpG islands
python -c """
import urllib.request
import json
import itertools

url = 'https://api.genome.ucsc.edu/getData/track?genome=hg38;track=cpgIslandExt'
with urllib.request.urlopen(url) as response:
    html = response.read()
    data = itertools.chain(*json.loads(html.decode())['cpgIslandExt'].values())
    for i in data:
        print(i['chrom'], i['chromStart'], i['chromEnd'], sep='\t')
""" > grch38.cpg_islands.bed
