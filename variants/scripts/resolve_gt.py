import sys
for line in sys.stdin:
    if line.startswith("##"):
        sys.stdout.write(line)
        continue
    data = line.strip().split('\t')
    if line.startswith("#CHR"):
        sys.stdout.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (from BPDP)">\n')
        sys.stdout.write('##FORMAT=<ID=FT,Number=1,Type=String,Description="Genotype filter">\n')
        data[9] = sys.argv[1]
        sys.stdout.write("\t".join(data[:10]) + '\n')
        continue
    
    # I need to check the coverage again...
    # But for both
    # Then I can do the AD/BPDP and set FT and stuff
    # so.. how.. there's no 
    data[8] = 'GT'#':AD:FT'
    gt, bpdp = data[9].split(':')
    comb_bpdp = bpdp
    bpdp_simple = [_ for _ in bpdp.split(',') if _ != '1']
    ad1 = '1' if len(bpdp_simple) == 0 else max(bpdp_simple)
    new_gt = ""
    if gt.startswith(("1/1", "1|1")):
        new_gt += "1"
    else:
        new_gt += "."
    new_gt += "|"

    gt, bpdp = data[10].split(':')
    comb_bpdp += ',' + bpdp
    bpdp_simple = [_ for _ in bpdp.split(',') if _ != '1']
    ad2 = '1' if len(bpdp_simple) == 0 else max(bpdp_simple)
    if gt.startswith(("1/1", "1|1")):
        new_gt += "1"
    else:
        new_gt += "."
    ad = ad1 + ',' + ad2
    if ad != "1,1":
        ft = 'FAIL'
    else:
        ft = 'PASS'
    #sys.stdout.write("\t".join(data[:9]) + "\t" + new_gt + ":" + comb_bpdp + ':' + ad + ':' + ft + "\n")
    sys.stdout.write("\t".join(data[:9]) + "\t" + new_gt + "\n")
