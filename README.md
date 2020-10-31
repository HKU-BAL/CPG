# CPG

# Introduction <br> 
The current human reference genome, GRCh38, is built from several individuals, most of whom are of Caucasian and African ancestry, which limits its usefulness for genomic analyses of distinct populations. For example, some population-specific variants cannot be detected when comparing genomes with the reference. Therefore, in recent years, the importance of capturing and representing sequencing data from diverse populations has been emphasized.  <br>
An ideal way to address the limitation is to create a pan-genome, a representation of both the core and variably distributed genomes of a species. Due to the computation complexity of assembling many deeply sequenced human genomes de novo and combining them into a genome and just a few gaps in GRCh38, we focused on finding large insertions, which made pan-genome assemble feasible.<br>


# Workflow for construction of a Chinese Pan-genome <br>
![Workflow](http://www.bio8.cs.hku.hk/RNA/CPG_workflow.png)<br> 
Firstly, we aligned the sequencing reads of 486 Han Chinese to the GRCh38.p13 reference genome individually and gained reads that could not be mapped to the reference genome. The unaligned reads were assembled into contigs (continuous sequences). Any contigs identified as contaminants or mapped to GRCh38 were eliminated. Based on the alignment positions of contigs’ reads and mates to GRCh38, we classified these contigs into two types: placed and unplaced ones. The exact insertion breakpoints of the placed contigs were determined and the placed contigs were then separated into three parts: right-end-placed (REP), left-end-placed (LEP) and both-end-placed (BEP). Secondly, we compared the placed contigs against one another and clustered together the similar contigs that placed close to each other. Unplaced contigs aligned closely to the placed ones with high identities were also included in the placed clusters. The remaining unplaced contigs were clustered by the cd-hit-est method. <br>

# Construction <br> 
## Data format <br> 
The name format of raw reads is Prefix+"_"+read_ID+"_"+sample_ID, like K14_17534_T1203 <br> 

## Step1. Assembly and filtering of novel contigs <br>
### 1.	Align reads to reference <br>
``` 
bwa index -p ref GRCh38_primary.fa 
bwa mem ref read1.fq read2.fq > alignment.sam
```  
### 2. Extract unaligned reads and corresponding mates <br>
``` 
samtools fastq -f 12 alignment.sam -1 R1_Unalignedmate.fq  -2 R2_Unalignedmate.fq  
samtools fastq -f 68 -F 8 alignment.sam > R1_alignedmate.fq  
samtools fastq -f 132 -F 8 alignment.sam > R2_alignedmate.fq  
samtools view -f 8 -F 4 alignment.sam > alignedmate_GRCh38.sam 
``` 
### 3. Assemble unaligned reads into contigs <br>
``` 
Megahit -r R1_Unalignedmate.fq , R2_Unalignedmate.fq, R1_alignedmate.fq, R2_alignedmate.fq  -o sample_1 
```  
### 4. 	Remove contaminations and contigs aligned to reference <br>
``` 
makeblastdb -in contaminations.fa -dbtype nucl -out contamination 
blastn -db contamination -query contig.fa -outfmt 6 -max_target_seqs 1  -max_hsps 1  -out  contig_contamination.tsv  
makeblastdb -in GRCh38_alt.fa -dbtype nucl -out ref_alt_Id 
blastn -db ref_alt_Id -query contig.fa -outfmt 6 -max_target_seqs 1  -max_hsps 1  -out  contig_ref.tsv 
```  
## Step2. Positioning of contigs in GRCh38  <br>
### 1.	Align reads to contigs <br>
``` 
bowtie2-build filteredcontig.fa contig_Id 
bowtie2 -x contig_Id -U R1_alignedmate.fq, R2_alignedmate.fq  -S readtocontig.sam 
``` 
### 2. Determine the placement region by reads and mates<br>
``` 
samtools view -h -F 2304 readtocontig.sam  | samtools sort -n -O bam | bedtools bamtobed -i stdin | awk '{OFS="\t"} {print $4,$1,$6,$2,$3}' | sed -e "s/\/[1-2]//g" |sort > readtocontig.txt 
samtools view -H alignedmate_GRCh38.sam | cat - <(awk 'NR==FNR{ a[$1]; next }$1 in a{ print $0 ; delete a[$1]; next }' readtocontig.txt <( samtools view alignedmate_GRCh38.sam )) | samtools sort -n -O bam | bedtools bamtobed -i stdin | awk '{OFS="\t"}{print $4,$1,$6,$2,$3}' | sed -e "s/\/[1-2]//g" | sort > pass_mates.txt 
join -j 1 readtocontig.txt pass_mates.txt > mates_region.txt 
``` 
### 3.	Examine links to contig ends only, and filter based on unambiguity criteria<br> 
&ensp;&ensp;Place_region.py<br> 

### 4. Extracted contig ends and GRCh38 regions with samtools faidx<br> 
``` 
samtools faidx GRCh38_no_alt.fa Place_region > GRCh38_Region.fa 
``` 
### 5. Align contigs to the region determined by the linking mates <br>
``` 
nucmer  --maxmatch -l 15 -b 1 -c 15 -p alignment_contig GRCh38Regions.fa end_contig.fa<br>
delta-filter -q -r -o 0 -g aliged_info.delta > filtered_info.delta <br> 
``` 
### 6. Obtain BEP/LEP/REP contigs and the corresponding placedment positions  <br> 
&ensp;&ensp;Contig_type.py <br> 
&ensp;&ensp;Please remove contigs that the both end aligned to reference from the LEP/REP file. The remaining contigs are unplaced. <br>

## Step3. Cluster placed contigs <br>
### 1.	Cluster placed contigs <br>
1.1.  Get the bed file (placed_contigs.sorted.bed)<br>
``` 
For BEP contigs, 
awk '{OFS="\t"} {split(FILENAME,b,"."); if($7=="reverse") print $2,$3-1,$5,$1"_"b[1],"-";  else print $2,$3-1,$5,$1"_"b[1],"+"}' BEP_folder/* |bedtools sort -i > BEP_contigs.bed 
For LEP/REP contigs, 
awk '{OFS="\t"} {split(FILENAME,b,"."); if($4=="reverse") print $2,$7-1,$8,$1"_"b[1],"-";  else print $2,$7-1,$8,$1"_"b[1],"+"}' LEP/REP_folder/* |bedtools sort -i > LEP/REP_contigs.bed
``` 
1.2.  Merge contigs in same type.<br>
``` 
bedtools merge -d 20 -c 4 -o distinct -i  placed_contigs.sorted.bed > merge_contigs.bed 
``` 
### 2. Choose the longest one as the representatives and get the corresponding clusters <br>
&ensp;&ensp;Rep_obtain.py <br>

### 3. Remove contigs with no alignments to representatives <br>
``` 
nucmer -p align_info  rep.fa cluster.fa<br>
``` 
### 4. Add other types of contigs to sequences in current clusters <br>
4.1.  Align contigs to sequences in the clusters.<br>
``` 
makeblastdb -in remaining_cluster.fa -dbtype nucl -out remainingcontigs_Id 
blastn -db remainingcontigs_Id -query othertype_contig.fa -outfmt "6  qseqid sseqid pident qlen slen length qstart qend sstart send mismatch g
apopen gaps evalue bitscore" -max_target_seqs 1  -max_hsps 1  -out  othertype_contig.tsv 
``` 
4.2.  Obtain contigs that can be added to the clusters.<br> 
 
&ensp;&ensp;Two types of contigs:
``` 
awk '{OFS="\t"}{if($3>99 && ($6-$13)/$4>=0.99 && ($6-$13) /$5>=0.8 ) print $2,$1}' othertype_contig.tsv > Ensure_contigs.txt
awk '{OFS="\t"}{if($3>99 && ($6-$13)/$5<0.8 && ($6-$13)/$4>=0.99 ) print $2,$1}' thertype_contig.tsv > candidate_contigs.txt
``` 
&ensp;&ensp;Get contigs that satisy two contiditions from the list of candidate contigs. <br> 
&ensp;&ensp;Pass_contigs.py <br> 
&ensp;&ensp; Output name: pass_contigs.txt

4.3.  Add other types of contigs into the current cluster (contigs from Ensure_contigs.txt and pass_contigs,txt)<br>
&ensp;&ensp;&ensp;&ensp;Move_contigs.py <br>

### 5. Merge left-end placed and right-end placed contigs into a longer insertion<br>
5.1. If an LEP contig and an REP contig were within 100 bp in the same orientation, please align the two contigs with each other. <br> 
``` 
nucmer -f  -p align_info left_placed.fa  right_placed.fa 
delta-filter -q  -r -g -m -1 align_info > filterdalign_info.delta 
show-coords -H -T -l -c -o filterdalign_info.delta > filterdalign_info.coords 
``` 

5.2. Classfiy the alignment result into four types:<br>
``` 
Identity:
awk '{OFS="\t"}{if ($NF=="[IDENTITY]") print $0}' filterdalign_info.coords | sort |uniq > Identity.txt 
Contained (the default value of identity_cutoff is 97): 
awk '{OFS="\t"}{if ($7>=identity_cutoff && ($NF=="[CONTAINED]" || $NF=="[CONTAINS]")) print $0}' filterdalign_info.coords |sort |uniq  > Contained.txt 
Overlap (the default value of identity_cutoff is 90 and the default value of minimun_cov_cutoff is 5 ): 
awk '{OFS="\t"}{if ($7>=identity_cutoff && $11>= minimun_cov_cutoff && $NF=="[END]") print $0}' filterdalign_info.coords |sort|uniq  > Overlap.txt 
Partially map (the default value of coverage_cutoff is 50): 
awk '{OFS="\t"}{if (($10>=coverage_cutoff || $11>=coverage_cutoff) && $NF!="[IDENTITY]" && $NF!="[CONTAINS]" && $NF!="[CONTAINED]") print $0}' filterdalign_info.coords|sort|uniq  > Part.txt 
``` 
**Noted:** <br>
For the fourth situation, please further check wehther there is at least one contig shared by the two clusters.<br>
``` 
nucmer -p Lrep_Rcluster  REP_cluster.fa LEP_rep.fa   
nucmer -p Rrep_Lcluster LEP_cluster.fa  REP_rep.fa 
delta-filter  -r -q -g LEP_rep_REP_cluster.delta > LEP_rep_REP_cluster_filter.delta 
delta-filter  -r -q -g REP_rep_LEP_cluster.delta > REP_rep_LEP_cluster_filter.delta 
show-coords -H -T -l -c -o LEP_rep_REP_cluster_filter.delta > LEP_rep_REP_cluster_filter.coords 
show-coords -H -T -l -c -o REP_rep_LEP_cluster_filter.delta > REP_rep_LEP_cluster_filter.coords  
```         
c. Merge pass LEP and REP contigs into one contigs.<br>
``` 
popins merge -c LEP_REP.fa <br>
``` 
d. Move reads with sev<br>

### 6. Remove the redundancy of placed contigs
``` 
makeblastdb -in all_placed.fa -dbtype nucl -out all_placed_Id<br>
blastn -db all_placed_Id -query all_placed.fa -outfmt "6  qseqid sseqid  pident slen qlen length qstart qend sstart send mismatch gapopen gaps evalue bitscore" -max_target_seqs 1  -max_hsps 1  -out  all_placed_aligned.tsv<br>
``` 
If contig 


### 7. Cluster the unplaced contigs<br>
``` 
cd-hit-est -i remain_unplaced.fa -o unplaced_cluster  -c 0.9 -n 8 <br>
``` 

## Step4. Analysis of CPG <br>
### 1. Component of CPG  <br>
The workflow provides two classification method.<br>
1. Placed / Unplaced <br>
2. Common / Individual-sepcific <br>

### 2. Call variants  
``` 
bwa index -p new_ref_Id  new_ref.fa<br>
bwa mem new_ref_Id read1.fq read2.fq > alignment.sam<br>
java -jar picard.jar MarkDuplicates I=alignment.sam O=alignment.markdup.sam M=alignment.markdup.txt<br>
java  -jar picard.jar BuildBamIndex I=alignment.markdup.sam<br>
gatk HaplotypeCallerSpark -R GRCh38_decoy.fa -I alignment.markdup.sam -O vcffile<br>
``` 
### 3. Align unaligned reads of 486 individuals to common sequences 
``` 
bwa index -p common_seq_id  common_seq.fa<br>
bwa mem common_seq_id unaligned_reads.fa > alignment.sam <br>
samtools view -h  -F 2304  alignment.sam  | htsbox samview -pS - > Filter_aligned.paf  <br>
``` 
### 4. Annotate placed contigs
``` 
vep -i contig_insertion_points.vcf -o contig_annotation --dir Cache_path --cache --offline --fasta GRCh38_primary.fa --species homo_sapiens --everything --plugin StructuralVariantOverlap,file=gnomad_v2_sv.sites.vcf.gz<br>
``` 
### 5. Compare with other genomes
``` 
bwa index -p other_genome_Id  other_genome.fa<br>
bwa mem other_genome_Id CPG.fa > alignment.sam<br>
``` 

