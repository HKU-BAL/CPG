# CPG

Step1. Assembly and filtering of novel contigs
•	Align reads to reference 
bwa index -p ref GRCh38_primary.fa
bwa mem ref read1.fq read2.fq > alignment.sam

•	Extract unaligned reads and corresponding mates 
samtools fastq -f 12 alignment.sam -1 R1_Unalignedmate.fq  -2 R2_Unalignedmate.fq 
samtools fastq -f 68 -F 8 alignment.sam > R1_alignedmate.fq
samtools fastq -f 132 -F 8 alignment.sam > R2_alignedmate.fq
samtools view -f 8 -F 4 alignment.sam > alignedmate_GRCh38.sam

•	Assemble unaligned reads into contigs 
Megahit -r R1_Unalignedmate.fq , R2_Unalignedmate.fq, R1_alignedmate.fq, R2_alignedmate.fq  -o sample_1 

•	Remove contaminations and contigs aligned to reference 
makeblastdb -in contaminations.fa -dbtype nucl -out contamination
blastn -db contamination -query contig.fa -outfmt 6 -max_target_seqs 1  -max_hsps 1  -out  contig_contamination.tsv
makeblastdb -in GRCh38_alt.fa -dbtype nucl -out ref_alt_Id
blastn -db ref_alt_Id -query contig.fa -outfmt 6 -max_target_seqs 1  -max_hsps 1  -out  contig_ref.tsv

Step2. Positioning of contigs in GRCh38

