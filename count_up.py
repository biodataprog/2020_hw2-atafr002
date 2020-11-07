#!/usr/bin/env python3


gff="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"
fasta="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz"

import os,gzip,itertools,csv,re

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
            sequence = b''.join(line.strip() for line in group)
            yield seq_id, sequence



if not os.path.exists(gff):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz")

if not os.path.exists(fasta):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz")

ch_length = 0
with gzip.open(fasta,"r") as f:
    pairs = aspairs(f)
    seqs  = dict(pairs)
    for k,v in seqs.items():
        ch_length += len(v)
    print("Chromosome length is: %d " %(ch_length)) 


gene_count = 0
sum_length = 0
CDS_length = 0

with gzip.open(gff,"rt") as fh:
    gff = csv.reader(fh,delimiter="\t")
    for row in gff:
        if row[0].startswith("#"):
            continue
        if row[2] == "gene":
            gene_count = gene_count + 1
            length = int(row[4]) - int(row[3])
            sum_length = sum_length + length
        if row[2] == "CDS":
            length = int(row[4]) - int(row[3])
            CDS_length = CDS_length + length
    print("Total number of genes: %s " %(str(gene_count)))
    print("Total number of nucleotides in the CDS: %s" %(str(sum_length)))


print("The percentage of genome which is coding: %f" %((CDS_length/ch_length)*100))
