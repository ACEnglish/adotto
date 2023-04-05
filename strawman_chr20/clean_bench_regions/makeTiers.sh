paste */refine.regions.txt > all_refine.bed
cut -f13,26,39 all_refine.bed | sort | uniq -c | sort -nr > states.txt
python make_tiers.py
