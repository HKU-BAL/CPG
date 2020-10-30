# CPG
 
**Step1. Assembly and filtering of novel contigs** <br>
•	Align reads to reference <br>
bwa index -p ref GRCh38_primary.fa <br>
bwa mem ref read1.fq read2.fq > alignment.sam
 
•	Extract unaligned reads and corresponding mates 
samtools fastq -f 12 alignment.sam -1 R1_Unalignedmate.fq  -2 R2_Unalignedmate.fq  <br> 
samtools fastq -f 68 -F 8 alignment.sam > R1_alignedmate.fq <br>
samtools fastq -f 132 -F 8 alignment.sam > R2_alignedmate.fq <br>
samtools view -f 8 -F 4 alignment.sam > alignedmate_GRCh38.sam <br>
 
•	Assemble unaligned reads into contigs 
Megahit -r R1_Unalignedmate.fq , R2_Unalignedmate.fq, R1_alignedmate.fq, R2_alignedmate.fq  -o sample_1 
 
•	Remove contaminations and contigs aligned to reference 
makeblastdb -in contaminations.fa -dbtype nucl -out contamination <br>
blastn -db contamination -query contig.fa -outfmt 6 -max_target_seqs 1  -max_hsps 1  -out  contig_contamination.tsv <br>
makeblastdb -in GRCh38_alt.fa -dbtype nucl -out ref_alt_Id <br>
blastn -db ref_alt_Id -query contig.fa -outfmt 6 -max_target_seqs 1  -max_hsps 1  -out  contig_ref.tsv <br>
 
**Step2. Positioning of contigs in GRCh38**
•	Align reads to contigs 
bowtie2-build filteredcontig.fa contig_Id<br>
bowtie2 -x contig_Id -U R1_alignedmate.fq, R2_alignedmate.fq  -S readtocontig.sam<br>
 
•	Determine the placement region by reads and mates
samtools view -h -F 2304 readtocontig.sam  | samtools sort -n -O bam | bedtools bamtobed -i stdin | awk '{OFS="\t"} {print $4,$1,$6,$2,$3}' | sed -e 's/\/[1-2]//g' |sort > readtocontig.txt
samtools view -H alignedmate_GRCh38.sam | cat - <(awk 'NR==FNR{ a[$1]; next }$1 in a{ print $0 ; delete a[$1]; next }' readtocontig.txt <( samtools view alignedmate_GRCh38.sam )) | samtools sort -n -O bam | bedtools bamtobed -i stdin | awk '{OFS="\t"}{print $4,$1,$6,$2,$3}' | sed -e 's/\/[1-2]//g' | sort > pass_mates.txt
join -j 1 readtocontig.txt pass_mates.txt > mates_region.txt
 
