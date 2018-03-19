# Motif Mark

Motif Mark assignment for BI 625 (Winter 2018).

This script creates a visualization of one to seven sequence motifs across a sequence or multiple sequences provided in UCSC format. Requires a FASTA file of the sequences and a list of motifs to be visualized.

Since a random color generator is used, if the output color scheme is unfavorable, run script again.

**Input**:

- FASTA file containing UCSC-formatted intron-exon-intron sequences to be searched for motifs
- Text file of sequence motifs to be searched for, formatted with one motif per line

**Output**:

- `sequences_twofer.fasta`: A FASTA file of your provided sequences that has two lines per entry (i.e. header and sequence)
- `motif_mark.svg`: An SVG visualization of where your motifs are found across your provided sequences