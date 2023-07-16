# UnidirectionalGeneBlockFinder
Dinoflagellate genomes commonly exhibit gene clusters with genes oriented in the same coding direction, resulting in the formation of unidirectional gene blocks. Recent studies have revealed the intriguing association of this unique gene organization with topologically associated domains (TADs) and transcription. Notably, TAD boundaries often coincide with the region between two converging unidirectional gene blocks and frequently display a distinct GC% dip. These findings highlight the significance of unidirectional gene blocks for investigating genome structure and gene regulation. Here, I developed a Python script capable of identifying unidirectional gene blocks and verifying the presence of potential GC% dips between them. See details in https://doi.org/10.3390/microorganisms10081662.

# Required library
Biopython \
matplotlib

# Usage
```
# extract gene records from gene annotation in gff3 format
awk '$3=="gene"' YOUR_GENE_ANNOTATION.gff3 | sort -k 1,1 -k4,4g > test.gene.gff3
python3 find_dino_tads.py test.gene.gff3 test.genome.fa 10  # requiring at least 10 genes in all unidirectional gene blocks.
```
