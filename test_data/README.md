# Test data for Dotplotic

This directory contains sample input files for testing **Dotplotic** program

## Files

- `Blast.tsv`
  
A BLAST output file in tabular format (format '6 std') generated through the following process: 
1) Only chromosomes extracted from three yeast genomes, *Saccharomyces cerevisiae* (GCF_000146045.2) as query, and *S. eubayanus* (GCF_001298625.1) and *S. pastorianus* (GCA_011022315.1) as subject.
2) A nucleotide BLAST search was performed, and alignment results shorter than 1kbp were filtered out.
   

```bash
seqkit seq -m 100000 GCA_011022315.1_ASM1102231v1_genomic.fna > CBS1483.fa
seqkit seq -m 100000 GCF_000146045.2_R64_genomic.fna > TwoYeastGenomes.fasta
seqkit seq -m 100000 GCF_001298625.1_SEUB3.0_genomic.fna >> TwoYeastGenomes.fasta
blastn -task blastn -query CBS1483.fa -subject TwoYeastGenomes.fasta -outfmt "6 std qlen slen" | awk '1000<$4 {print}' > Blast.tsv
```

Note: Only Blast.tsv file exists here.

- `CDS.gff` and `tRNA.gff`
  
These files are extracted from the GFF file of the GenBank assembly R64 (`GCF_000146045.2_R64_genomic.gff`). 

```bash 
cat GCF_000146045.2_R64_genomic.gff | grep -v "^#" | awk '$3=="CDS" {print}' > CDS.gff
cat GCF_000146045.2_R64_genomic.gff | grep -v "^#" | awk '$3=="tRNA" {print}' > tRNA.gff
```

- `query.bed` and `subject.bed`
  
The start and end positions of `NC_001136.10` is `200000` and `0` in the file `subject.bed`, respectively, so all BLAST alignments are horizontally inverted in the plot.

