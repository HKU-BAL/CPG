# CPG

# Introduction <br> 
The current human reference genome, GRCh38, is built from several individuals, most of whom are of Caucasian and African ancestry, which limits the genomic analyses of distinct populations. For example, some population-specific variants cannot be detected when comparing genomes with the reference.  <br>
An ideal way to address the limitation is to create a pan-genome, a representation of both the core and variably distributed genomes of a species. Pan-genomic analysis better captures unexplored or missed variants to improve the decoding of the genetic basis of human diseases. Considering the computation complexity of assembling deep sequenced human genomes _de novo_ and combining them, a feasible strategy for improved references is to focus on large insertion and construct population-specific pan-genomes.<br>


# Workflow for construction of a Chinese Pan-genome <br>
![Workflow](http://www.bio8.cs.hku.hk/RNA/CPG_workflow.png)<br> 
Firstly, we aligned the sequencing reads of 486 Han Chinese to the GRCh38.p13 reference genome individually and gained reads that could not be mapped to the reference genome. The unaligned reads were assembled into contigs (continuous sequences). Any contigs identified as contaminants or mapped to GRCh38 were eliminated. Based on the alignment positions of contigs’ reads and mates to GRCh38, we classified these contigs into two types: placed and unplaced ones. The exact insertion breakpoints of the placed contigs were determined and the placed contigs were then separated into three parts: right-end-placed (REP), left-end-placed (LEP) and both-end-placed (BEP). Secondly, we compared the placed contigs against one another and clustered together the similar contigs that placed close to each other. Unplaced contigs aligned closely to the placed ones with high identities were also included in the placed clusters. Finally, the remaining unplaced contigs were clustered by cd-hit-est. <br>

# Construction <br> 
## Data format <br> 
The name format of raw reads is `Prefix+"_"+read_ID`, like K14_17534. <br> 

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
### 1.	Align reads to filtered contigs <br>
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
Based on the alignment information, get the unambiguous placement for each contig. <br>

### 3. Extract contig ends and GRCh38 regions<br> 
``` 
# For files in unambiguous_placed_regions_folder/LEP or unambiguous_placed_regions_folder/REP folder, 
awk '{print $2":"$3"-"$4}' unambiguous_placed_regions_folder/Placed/contig_ID.txt > contig_ID_LEP/REP_region.txt
samtools faidx GRCh38_no_alt.fa contig_ID_LEP_region.txt > GRCh38_Region.fa 
samtools faidx contig_ID.fa region > LEP/REP_contig.fa

``` 
### 4. Align contigs to the region determined by the linking mates <br>
``` 
nucmer  --maxmatch -l 15 -b 1 -c 15 -p contig_ID GRCh38_Regions.fa REP_contig.fa/LEP_contig.fa  
delta-filter -q -r -o 0 -g contig_ID.delta > filtered_info.delta 
``` 
### 5. Determine BEP/LEP/REP contigs and the corresponding placedment positions  <br> 
``` 
python contig_type.py  --ref_name_id GRCH38.fa.fai --alignment_info PATH_filtered_info.delta  --LEP_contigs LEP_folder --REP_contigs REP_folder --BEP_contigs BEP_folder --BEP_contigs_all all_BEP_folder
```
Please remove contigs in BEP_contigs_all folder from the LEP/REP folder. The remaining contigs are unplaced. <br>

## Step3. Cluster placed contigs <br>
### 1.	Cluster placed contigs <br>
1.1.  Get the placement locations of contigs <br>
``` 
# For BEP contigs, 
awk '{OFS="\t"} {split(FILENAME,b,"."); if($7=="reverse") print $2,$3-1,$5,$1"_"b[1],"-";  else print $2,$3-1,$5,$1"_"b[1],"+"}' BEP_folder/* |bedtools sort -i > BEP_contigs.bed 

# For LEP/REP contigs, 
awk '{OFS="\t"} {split(FILENAME,b,"."); if($4=="reverse") print $2,$7-1,$8,$1"_"b[1],"-";  else print $2,$7-1,$8,$1"_"b[1],"+"}' LEP/REP_folder/* |bedtools sort -i > LEP/REP_contigs.bed
``` 
1.2.  Group contigs based on placement positions<br>
``` 
bedtools merge -d 20 -c 4 -o distinct -i  placed_contigs.sorted.bed > merge_contigs.bed 
``` 
### 2. Choose the longest one as the representatives and get the corresponding clusters <br>
``` 
python rep_obtain.py --seq_path LEP/REP/BEP_seq_path --path_merge_bed merge_contigs.bed --path_rep save_rep_folder --path_contig save_cluster_folder
``` 
### 3. Remove contigs with no alignments to representatives <br>
``` 
nucmer -p align_info  rep.fa cluster.fa<br>
``` 
Only save contigs that hit to representative (save folder: `remain_cluster_folder`) <br>

### 4. Move other types of contigs to sequences to the current clusters <br>
4.1.  Align contigs to sequences in the clusters <br>
``` 
makeblastdb -in remaining_cluster.fa -dbtype nucl -out remainingcontigs_Id 
blastn -db remainingcontigs_Id -query othertype_contig.fa -outfmt "6  qseqid sseqid pident qlen slen length qstart qend sstart send mismatch g
apopen gaps evalue bitscore" -max_target_seqs 1  -max_hsps 1  -out  othertype_contig.tsv 
``` 
4.2.  Obtain contigs that can be added to the clusters <br> 
``` 
# Two types of contigs
awk '{OFS="\t"}{if($3>99 && ($6-$13)/$4>=0.99 && ($6-$13) /$5>=0.8 ) print $2,$1}' othertype_contig.tsv > Ensure_contigs.txt
awk '{OFS="\t"}{if($3>99 && ($6-$13)/$5<0.8 && ($6-$13)/$4>=0.99 ) print $2,$1}' othertype_contig.tsv > candidate_contigs.txt
``` 
Remove some candidate contigs based on the placement of their linking mates. The remaining contigs are saved in `pass_contigs.txt`, file format: `Rep_ID  pass_contig_ID` <br> 

4.3.  Add other types of contigs to the current cluster (contigs from Ensure_contigs.txt and pass_contigs,txt)<br>
```
python move_contigs.py --ensure_contigs -Ensure_contigs.txt -pass_contigs pass_contigs.txt --cluster_folder remain_cluster_folder --contig_path othertype_contig_path
```
### 5. Merge left-end placed and right-end placed contigs into a longer insertion<br>
5.1. For two contigs within 100 bp in the same orientation, align the two contigs <br> 
``` 
nucmer -f  -p align_info left_placed.fa  right_placed.fa 
delta-filter -q  -r -g -m -1 align_info > filterdalign_info.delta 
show-coords -H -T -l -c -o filterdalign_info.delta > filterdalign_info.coords  
``` 
Classify the alignment result into four situtaions:<br>
situation1. The two representatives are identical `identity.coords` <br> 
situation2. One representative is contained, the identity cutoff is over 90%  `contained.coords` <br>
situation3. The ends of two representatives overlap in the correct arrangement and orientation `overlap.coords` <br>
situation4. One representatives covering at least 50% of the other one `part.coords`<br>

5.2. Update the alignment result<br>
``` 
python reorg_align_info.py  --LEP_bed LEP_contigs.bed --REP_bed REP_contig.bed --Identity_path identity.coords --Contained_path contained.coords  --Overlap_path overlap.coords --Part_align_path part.coords --save_folder updated_alignment_folder
``` 

5.3 Remove false-positive potentially merging clusters.
for the fourth alignment result, further check wehther there is at least one contig shared by the two clusters. 
``` 
nucmer -p Lrep_Rcluster  REP_cluster.fa LEP_rep.fa   
nucmer -p Rrep_Lcluster LEP_cluster.fa  REP_rep.fa 
delta-filter  -r -q -g LEP_rep_REP_cluster.delta > LEP_rep_REP_cluster_filter.delta 
delta-filter  -r -q -g REP_rep_LEP_cluster.delta > REP_rep_LEP_cluster_filter.delta 
show-coords -H -T -l -c -o LEP_rep_REP_cluster_filter.delta > LEP_rep_REP_cluster_filter.coords 
show-coords -H -T -l -c -o REP_rep_LEP_cluster_filter.delta > REP_rep_LEP_cluster_filter.coords  
```         
Only REP_rep_LEP_cluster_filter/LEP_rep_REP_cluster_filter.coords reports `CONTAINED/IDENTITY` can the alignment result be saved 
Results are saved in `updated_alignment_folder/final_part.txt`, file format: `LEP_rep REP_rep New_rep.{r/l}`. <br>

5.3. Merge overlaping a LEP and REP representative into one longer contig <br>
``` 
popins merge -c LEP_REP.fa <br>
``` 
Save files in `merged/` 

5.4.  Update the the placed representatives and corresponding clusters <br>
``` 
python update_ref.py --LEP_folder LEP_folder/ --REP_folder  REP_folder/  --contigs_fai  contigs_fai_path  --LEP_cluster_folder LEP_cluster_folder --REP_cluster_folder  --LEP_rep_folder LEP_rep/ --REP_rep_folder REP_rep/ --align_folder updated_alignment_folder/  --LEP_cluster_update_folder LEP_cluster_update/ --LEP_rep_update_folder LEP_rep_update --REP_cluster_update_folder REP_cluster_update/ --REP_rep_update_folder REP_rep_update/ --BEP_rep_folder BEP_rep/ --BEP_cluster_folder BEP_cluster/ --merged_contig_folder merged/
``` 

### 6. Remove the redundancy of placed contigs
6.1 Align placed contigs against each other <br>
``` 
makeblastdb -in all_placed.fa -dbtype nucl -out all_placed_Id 
blastn -db all_placed_Id -query all_placed.fa -outfmt "6  qseqid sseqid  pident slen qlen length qstart qend sstart send mismatch gapopen gaps evalue bitscore" -max_target_seqs 1  -max_hsps 1  -out  all_placed_aligned.tsv
``` 
6.2 Obtain contigs that can be merged together  <br>
``` 
python deduplcate_placed.py  --alignment_path  all_placed_aligned.tsv  --BEP_bed  BEP_bed_path --LEP_bed  LEP_bed_path --REP_bed  REP_bed_path --pass_alignment  placed_aligned.update.tsv
``` 
6.3 Generate the final placed representatives and clusters <br> 
``` 
python final_ref.py  --placed_align_path final_aligned.txt --contigs_fai contigs.fa.fai --BEP_cluster_folder  BEP_cluster_updater/  --LEP_cluster_folder  LEP_cluster_update/  --REP_cluster_folder  REP_cluster_update/ --BEP_rep_folder  BEP_rep_update/ --LEP_rep_folder LEP_rep_update/  --REP_rep_folder REP_rep_update/  --BEP_cluster_update_folder  BEP_final_cluster/ --LEP_cluster_update_folder  LEP_final_cluster/  --REP_cluster_update_folder REP_final_cluster/ --BEP_rep_update_folder  BEP_final_rep/ --LEP_rep_update_folder LEP_final_rep/  --REP_rep_update_folder  REP_final_rep/
``` 

### 7. Cluster the unplaced contigs <br>
``` 
cd-hit-est -i remain_unplaced.fa -o unplaced_cluster  -c 0.9 -n 8 <br>
``` 

## Step4. Analysis of CPG <br>
### 1. Call variants  
``` 
bwa index -p new_ref_Id  new_ref.fa<br>
bwa mem new_ref_Id read1.fq read2.fq > alignment.sam
java -jar picard.jar MarkDuplicates I=alignment.sam O=alignment.markdup.sam M=alignment.markdup.txt
java  -jar picard.jar BuildBamIndex I=alignment.markdup.sam
gatk HaplotypeCallerSpark -R GRCh38_decoy.fa -I alignment.markdup.sam -O vcffile
``` 
### 2. Align unaligned reads of 486 individuals to common sequences 
``` 
bwa index -p common_seq_id  common_seq.fa
bwa mem common_seq_id unaligned_reads.fa > alignment.sam 
samtools view -h  -F 2304  alignment.sam  | htsbox samview -pS - > Filter_aligned.paf  
``` 
### 3. Annotate placed contigs
``` 
vep -i contig_insertion_points.vcf -o contig_annotation --dir Cache_path --cache --offline --fasta GRCh38_primary.fa --species homo_sapiens --everything --plugin StructuralVariantOverlap,file=gnomad_v2_sv.sites.vcf.gz
``` 
### 4. Compare with other genomes
``` 
bwa index -p other_genome_Id  other_genome.fa
bwa mem other_genome_Id CPG.fa > alignment.sam
``` 

