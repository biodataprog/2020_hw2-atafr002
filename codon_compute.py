#!/usr/bin/env python3

import os, gzip, itertools

def isheader(line):
    return line[0] == '>'

def aspairs(f):
    seq_id = ''
    sequence = ''
    for header,group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence

url1="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2/cds/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
url2="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/cds/Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"
file1="Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
file2="Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"

if not os.path.exists(file1):
    os.system("curl -O %s"%(url1))

if not os.path.exists(file2):
    os.system("curl -O %s"%(url2))

with gzip.open(file1,"rt") as fh:
	seqs = aspairs(fh)

	total_genes1 = 0
	gene_length1 = 0
	G1 = 0
	C1 = 0
	codon1_dict = {'TTT':0 , 'TTC': 0, 'TTA' : 0, 'TTG' : 0, 'CTT': 0, 'CTC': 0, 'CTA': 0, 'CTG': 0,'ATT': 0,
	'ATC': 0,'ATA': 0,'ATG': 0,'GTT': 0,'GTC': 0,'GTA': 0,'GTG': 0,'TCT': 0,'TCC': 0,'TCA': 0,'TCG': 0,'CCT': 0,'CCC': 0,
	'CCA': 0,'CCG': 0,'ACT': 0,'ACC': 0,'ACA': 0,'ACG': 0,'GCT': 0,'GCC': 0,'GCA': 0,'GCG': 0,'TAT': 0,'TAC': 0,'TAA': 0,
	'TAG': 0,'CAT': 0,'CAC': 0,'CAA': 0,'CAG': 0,'AAT': 0,'AAC': 0,'AAA': 0,'AAG': 0,'GAT': 0,'GAC': 0,'GAA': 0,'GAG': 0,
	'TGT': 0,'TGC': 0,'TGA': 0,'TGG': 0,'CGT': 0,'CGC': 0,'CGA': 0,'CGG': 0,'AGT': 0,'AGC': 0,'AGA': 0,'AGG': 0,'GGT': 0,
	'GGC': 0,'GGA': 0,'GGG': 0}
	for seq in seqs:
		seqname  = seq[0]
		seqstring= seq[1]
		codon = []
#Counts the total number of genes:
		total_genes1 += 1
#Counts the total gene length:
		gene_length1 += len(seqstring)

#Counts the number of occurance of G and C:
		for i in seqstring:
			if i == "G":
				G1 += 1
			if i == "C":
				C1 += 1
			else:
				continue

#Updates the dictionary for the codons:
		for i in range(0,len(seqstring),3):
			codon1_dict[str(seqstring[i:i+3])] += 1


with gzip.open(file2, "rt") as f:
	seqs = aspairs(f)

	total_genes2 = 0
	gene_length2 = 0
	G2 = 0
	C2 = 0
	codon2_dict = {'TTT':0 , 'TTC': 0, 'TTA' : 0, 'TTG' : 0, 'CTT': 0, 'CTC': 0, 'CTA': 0, 'CTG': 0,'ATT': 0,
        'ATC': 0,'ATA': 0,'ATG': 0,'GTT': 0,'GTC': 0,'GTA': 0,'GTG': 0,'TCT': 0,'TCC': 0,'TCA': 0,'TCG': 0,'CCT': 0,'CCC': 0,
        'CCA': 0,'CCG': 0,'ACT': 0,'ACC': 0,'ACA': 0,'ACG': 0,'GCT': 0,'GCC': 0,'GCA': 0,'GCG': 0,'TAT': 0,'TAC': 0,'TAA': 0,
        'TAG': 0,'CAT': 0,'CAC': 0,'CAA': 0,'CAG': 0,'AAT': 0,'AAC': 0,'AAA': 0,'AAG': 0,'GAT': 0,'GAC': 0,'GAA': 0,'GAG': 0,
        'TGT': 0,'TGC': 0,'TGA': 0,'TGG': 0,'CGT': 0,'CGC': 0,'CGA': 0,'CGG': 0,'AGT': 0,'AGC': 0,'AGA': 0,'AGG': 0,'GGT': 0,
        'GGC': 0,'GGA': 0,'GGG': 0}

	for seq in seqs:
		seqname = seq[0]
		seqstring= seq[1]
		total_genes2 += 1
		gene_length2 += len(seqstring)

		for i in seqstring:
			if i == "G":
				G2 += 1
			if i == "C":
				C2 += 1
			else:
				continue

		for i in range(0,len(seqstring),3):
			codon2_dict[str(seqstring[i:i+3])] += 1


print("\n Question3_Part1")
print("Total number of genes for Salmonella enterica: %d" %(total_genes1))
print("Total number of genes for Mycobacterium tuberculosis: %d" %(total_genes2))   


print("\n Question3_Part2")
print("Total length of genes for Salmonella enterica: %d" %(gene_length1))
print("Total length of genes for Mycobacterium tuberculosis: %d" %(gene_length2))

print("\n Question3_Part3")
print("The G+C percentage for the whole dataset: %f" %((G1 + C1 + G2 + C2)/(gene_length1 + gene_length2)*100))

print("\n Question3_Part4")
print("Total number of codons in Salmonella enterica's genome: %d" %((gene_length1)/3))
print("Total number of codons in Mycobacterium tuberculosis's genome: %d" %((gene_length2)/3))

print("\n Question3_Part5")

for k,v in codon1_dict.items():
	print(k + "\t|\t" + str(codon1_dict[k]) + "\t|\t" + str(codon2_dict[k]) + "\n --------------------------------------")

