import argparse as ap
import difflib
from scipy import stats
import gzip
import io

#get samples from index


#polarize to ancestral or major allele freq 
def get_polarized_genotypes(line, index1, index2, index_out):
	
###############
	selected_samples_out = [line[i] for i in index_out]
	samples_genotypes_out = []
	for sample in selected_samples_out:
		sample = sample.split(':')[0]
		if sample == '.':
			samples_genotypes_out.append('.')
			samples_genotypes_out.append('.')
		else:
			sample = sample.split('/')
			if len(sample) == 1:
				samples_genotypes_out.append(sample[0])
				samples_genotypes_out.append(sample[0])
			else:
				samples_genotypes_out.append(sample[0])
				samples_genotypes_out.append(sample[1])

###############
	selected_samples1 = [line[i] for i in index1]
	samples_genotypes1 = []
	for sample in selected_samples1:
		sample = sample.split(':')[0]
		if sample == '.':
			samples_genotypes1.append('.')
			samples_genotypes1.append('.')
		else:
			sample = sample.split('/')
			if len(sample) == 1:
				samples_genotypes1.append(sample[0])
				samples_genotypes1.append(sample[0])
			else:
				samples_genotypes1.append(sample[0])
				samples_genotypes1.append(sample[1])

###############
	selected_samples2 = [line[i] for i in index2]
	samples_genotypes2 = []
	for sample in selected_samples2:
		sample = sample.split(':')[0]
		if sample == '.':
			samples_genotypes2.append('.')
			samples_genotypes2.append('.')
		else:
			sample = sample.split('/')
			if len(sample) == 1:
				samples_genotypes2.append(sample[0])
				samples_genotypes2.append(sample[0])
			else:
				samples_genotypes2.append(sample[0])
				samples_genotypes2.append(sample[1])

	out_geno = [int(i) for i in samples_genotypes_out if i != '.']
	p1_geno = [int(i) for i in samples_genotypes1 if i != '.']
	p2_geno = [int(i) for i in samples_genotypes2 if i != '.']
	p12_geno = p1_geno+p2_geno

	total_alleles = p1_geno+p2_geno+out_geno

	###Version 27/03/2021 - More precise when there is no info for the outgroup
	
	if len(set(out_geno)) == 0 and len(set(p1_geno+p2_geno)) == 0: #all missing
		ref_allele = 9
		flag = 'allMiss'

	elif len(set(p1_geno+p2_geno)) == 0: #ingroup species missing
		ref_allele = int(stats.mode(out_geno)[0])
		flag = 'inMiss'

	elif len(set(out_geno)) == 0 and len(set(p1_geno+p2_geno)) > 0: #outgroup is useless
		if len(set(p1_geno)) == 0:
			ref_allele = int(stats.mode(p2_geno)[0]) #not really relevant as the target pop1 will be made of missing data
			flag = 'in2Fold'
		elif len(set(p2_geno)) == 0: #also p2 is useless
			ref_allele = int(stats.mode(p1_geno)[0])
			flag = 'in1Fold'
		elif len(set(p1_geno)) == 1 and len(set(p2_geno)) == 1: 
			ref_allele = p1_geno[0] #pop1 allele as reference = conservative. See also below
			flag = 'InFixOutMiss'
		elif len(set(p1_geno)) == 2 and len(set(p2_geno)) == 1:
			ref_allele = p2_geno[0]
			flag = 'unfoldOutMiss'
		elif len(set(p1_geno)) == 1 and len(set(p2_geno)) == 2:
			ref_allele = p1_geno[0]
			flag = 'unfoldOutMiss'
		elif len(set(p1_geno)) == 2 and len(set(p2_geno)) == 2:
			ref_allele = int(stats.mode(p1_geno+p2_geno)[0]) ### not random! Will always be 0 in case of 50/50
			flag = 'inFold'

	elif len(set(out_geno)) == 1 and len(set(p1_geno+p2_geno)) == 1: 
		ref_allele = p12_geno[0] #pop1 allele as reference = conservative.
		flag = 'allFix'
	elif len(set(out_geno)) == 1 and len(set(p1_geno+p2_geno)) == 2:
		ref_allele = out_geno[0]
		flag = 'unfolded'

	elif len(set(out_geno)) == 2 and len(set(p1_geno+p2_geno)) == 1:
		ref_allele = p12_geno[0]
		flag = 'InFixAnc'
	elif len(set(out_geno)) == 2 and len(set(p1_geno+p2_geno)) == 2:
		ref_allele = int(stats.mode(p1_geno+p2_geno+out_geno)[0]) ### not random! Will always be 0 in case of 50/50
		flag = 'allFold'

	else:
		ref_allele = 9
		print('WARNING: Non biallelic locus in the vcf. Replaced with missing value for all samples')
		print(p1_geno, p2_geno, out_geno)
		flag = 'triall'

	p1_geno_polarized = []
	if ref_allele == 0:
		for i in samples_genotypes1:
			if i == '.':
				p1_geno_polarized.append(9)
			elif int(i) == 0:
				p1_geno_polarized.append(0)
			elif int(i) == 1:
				p1_geno_polarized.append(1)
			elif int(i) == 2:
				p1_geno_polarized.append(2)###just in case
			elif int(i) == 3:
				p1_geno_polarized.append(3)###just in case
	elif ref_allele == 1:
		for i in samples_genotypes1:
			if i == '.':
				p1_geno_polarized.append(9)
			elif int(i) == 0:
				p1_geno_polarized.append(1)
			elif int(i) == 1:
				p1_geno_polarized.append(0)
			elif int(i) == 2:
				p1_geno_polarized.append(2)###just in case
			elif int(i) == 3:
				p1_geno_polarized.append(3)###just in case
	else:
		for i in samples_genotypes1:
			p1_geno_polarized.append(9)

	p2_geno_polarized = []
	if ref_allele == 0:
		for i in samples_genotypes2:
			if i == '.':
				p2_geno_polarized.append(9)
			elif int(i) == 0:
				p2_geno_polarized.append(0)
			elif int(i) == 1:
				p2_geno_polarized.append(1)
			elif int(i) == 2:
				p2_geno_polarized.append(2)###just in case
			elif int(i) == 3:
				p2_geno_polarized.append(3)###just in case
	elif ref_allele == 1:
		for i in samples_genotypes2:
			if i == '.':
				p2_geno_polarized.append(9)
			elif int(i) == 0:
				p2_geno_polarized.append(1)
			elif int(i) == 1:
				p2_geno_polarized.append(0)
			elif int(i) == 2:
				p2_geno_polarized.append(2)###just in case
			elif int(i) == 3:
				p2_geno_polarized.append(3)###just in case
	else:
		for i in samples_genotypes2:
			p2_geno_polarized.append(9)

	out_geno_polarized = []
	if ref_allele == 0:
		for i in samples_genotypes_out:
			if i == '.':
				out_geno_polarized.append(9)
			elif int(i) == 0:
				out_geno_polarized.append(0)
			elif int(i) == 1:
				out_geno_polarized.append(1)
			elif int(i) == 2:
				out_geno_polarized.append(2)###just in case
			elif int(i) == 3:
				out_geno_polarized.append(3)###just in case
	elif ref_allele == 1:
		for i in samples_genotypes_out:
			if i == '.':
				out_geno_polarized.append(9)
			elif int(i) == 0:
				out_geno_polarized.append(1)
			elif int(i) == 1:
				out_geno_polarized.append(0)
			elif int(i) == 2:
				out_geno_polarized.append(2)###just in case
			elif int(i) == 3:
				out_geno_polarized.append(3)###just in case
	else:
		for i in samples_genotypes_out:
			out_geno_polarized.append(9)

	return p1_geno_polarized, ref_allele, total_alleles, p2_geno_polarized, out_geno_polarized, flag


def replace_all(text):
	table = {'0/0':'1|1','0/1':'1|0','1/0':'0|1','1/1':'0|0'}
	for i, j in table.items():
		text = text.replace(i, j)
	return text


parser = ap.ArgumentParser()
parser.add_argument('-a', '--ann', help='An annotated gzipped vcf file', required=True, type=str)
parser.add_argument('-p1', '--popin1', help='Provide txt file with one individual per line', required=True, type=str)
parser.add_argument('-p2', '--popin2', help='Provide txt file with one individual per line', required=True, type=str)
parser.add_argument('-p3', '--popout', help='Provide txt file with one individual per line', required=True, type=str)
args = parser.parse_args()

#define input file
vcf = args.ann
scaf = args.ann.replace('.ann.vcf.gz','')
output = scaf+'.'+args.popin1

out_annFreq = open(output+'.annFreq','w')

#fields in the output
out_annFreq.write('scaffold\tposition\tflagPol\tflagQual\tderived1\ttotal1\tderived2\ttotal2\tderived_out\ttotal_out\tavgCov\tref\talt\tvartype\teffect\n')

#get the header
with gzip.open(vcf, 'rb') as input_file:
	with io.TextIOWrapper(input_file, encoding='utf-8') as handle:
		for line in handle:
			line = line.rstrip()
			if line.startswith('##'):
				continue
			elif line.startswith('#CHROM'):
				header = line.rstrip().split('\t')
			else:
				break

#get index of ingroup and outgroup samples
samples1 = []
for line in open(args.popin1, 'r').readlines():
    if line != '\n' :
        samples1.append(line.strip())
index_samples1 = []
for sample in samples1:
	index_sample = header.index(sample)
	index_samples1.append(index_sample)

samples2 = []
for line in open(args.popin2, 'r').readlines():
    if line != '\n' :
        samples2.append(line.strip())
index_samples2 = []
for sample in samples2:
	index_sample = header.index(sample)
	index_samples2.append(index_sample)

samples_out = []
for line in open(args.popout, 'r').readlines():
    if line != '\n' :
        samples_out.append(line.strip())

index_samples_out = []
for sample in samples_out:
	index_sample = header.index(sample)
	index_samples_out.append(index_sample)


oldqual = 0
#go through the vcf
with gzip.open(vcf, 'rb') as input_file:
	with io.TextIOWrapper(input_file, encoding='utf-8') as handle:
		for line in handle:
			if line.startswith('#') or line == '':
				continue
			else:

				line = line.rstrip()
				line = line.replace('|','/') # to correct for phased+unphased
				line = line.split('\t')
				scaffold = line[0]
				pos = int(line[1])
				ref = line[3]
				alt = line[4]
				qual = line[5]
			
				if qual == oldqual:
					flagQ = 'haplo'
				else:
					flagQ = 'snp'

				oldqual = qual

				polarized = get_polarized_genotypes(line,index_samples1,index_samples2,index_samples_out)
				derived1 = polarized[0].count(1)
				derived2 = polarized[3].count(1)
				derived3 = polarized[4].count(1)
				anc_allele = polarized[1]
				tot_alleles = polarized[2]
				p1_total_allele = len([i for i in polarized[0] if i != 9])
				p2_total_allele = len([i for i in polarized[3] if i != 9])
				p3_total_allele = len([i for i in polarized[4] if i != 9])
				fl = polarized[5]

				if anc_allele == 1:
					newref = alt
					newalt = ref

				elif anc_allele == 0:
					newref = ref
					newalt = alt

				if len(tot_alleles) < len(index_samples1+index_samples2+index_samples_out)/2:
					continue
				else:
					avgCov = round(int(line[7].split(';')[7].replace('DP=', '').split(',')[0])/len(tot_alleles), 1)
					snpeff = line[7].split(';')[42].split('/')
					annotation = snpeff[1]
					effect = snpeff[2]
					if 'missense' in annotation:
						vartype = 'missense'
					elif 'synonymous' in annotation:
						vartype = 'synonymous'
					elif 'intergenic' in annotation:
						vartype = 'intergenic'
					elif 'intron' in annotation:
						vartype = 'intron'
					else:
						vartype = 'else'

					out_annFreq.write(scaffold+'\t'+str(pos)+'\t'+fl+'\t'+flagQ+'\t'+str(derived1)+'\t'+str(p1_total_allele)+'\t'+str(derived2)+'\t'+str(p2_total_allele)+'\t'+str(derived3)+'\t'+str(p3_total_allele)+'\t'+str(avgCov)+'\t'+newref+'\t'+newalt+'\t'+vartype+'\t'+effect+'\n')

out_annFreq.close()


