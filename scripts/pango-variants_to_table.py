import sys
import vcf
import vcf.utils
import pandas
import os.path

# sites
output = sys.argv[1]
sites_f = sys.argv[2]
vcf_files = sys.argv[3:]

sites = set()
reader = vcf.Reader(open(sites_f, 'rb'))
for i, r in enumerate(reader):
	assert len(r.ALT)==1
	sites.add((r.CHROM, r.POS, str(r.REF), str(r.ALT[0])))

# SAMPLE, CHROM, POS, REF, ALT, ADREF, ADALT, DP
muts = []
for f in vcf_files:
	reader = vcf.Reader(open(f, 'rb'))
	for i, r in enumerate(reader):
		for j, mp in enumerate(r.ALT):
			for k, s in enumerate(reader.samples):
				ref, alt = vcf.utils.trim_common_suffix(str(r.REF), str(r.ALT[j]))
				mut = (r.CHROM, r.POS, ref, alt)
				if mut in sites:
					muts.append([s, r.CHROM, r.POS, ref, alt, r.samples[k]['AD'][0], r.samples[k]['AD'][j+1], r.samples[k]['DP']])

df = pandas.DataFrame(muts, columns=['sample','chrom','pos','ref','alt','adref','adalt','dp'])
df.to_csv(output, sep='\t', index=False)
