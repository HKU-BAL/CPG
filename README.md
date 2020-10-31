# CPG
 
# Step1. Assembly and filtering of novel contigs <br>
•	Align reads to reference <br>
bwa index -p ref GRCh38_primary.fa <br>
bwa mem ref read1.fq read2.fq > alignment.sam
 
•	Extract unaligned reads and corresponding mates <br>
samtools fastq -f 12 alignment.sam -1 R1_Unalignedmate.fq  -2 R2_Unalignedmate.fq  <br> 
samtools fastq -f 68 -F 8 alignment.sam > R1_alignedmate.fq <br>
samtools fastq -f 132 -F 8 alignment.sam > R2_alignedmate.fq <br>
samtools view -f 8 -F 4 alignment.sam > alignedmate_GRCh38.sam <br>
 
•	Assemble unaligned reads into contigs <br>
Megahit -r R1_Unalignedmate.fq , R2_Unalignedmate.fq, R1_alignedmate.fq, R2_alignedmate.fq  -o sample_1 <br>
 
•	Remove contaminations and contigs aligned to reference <br>
makeblastdb -in contaminations.fa -dbtype nucl -out contamination <br>
blastn -db contamination -query contig.fa -outfmt 6 -max_target_seqs 1  -max_hsps 1  -out  contig_contamination.tsv <br>
makeblastdb -in GRCh38_alt.fa -dbtype nucl -out ref_alt_Id <br>
blastn -db ref_alt_Id -query contig.fa -outfmt 6 -max_target_seqs 1  -max_hsps 1  -out  contig_ref.tsv <br>
 
# Step2. Positioning of contigs in GRCh38  <br>
•	Align reads to contigs <br>
bowtie2-build filteredcontig.fa contig_Id<br>
bowtie2 -x contig_Id -U R1_alignedmate.fq, R2_alignedmate.fq  -S readtocontig.sam<br>
 
•	Determine the placement region by reads and mates<br>
samtools view -h -F 2304 readtocontig.sam  | samtools sort -n -O bam | bedtools bamtobed -i stdin | awk '{OFS="\t"} {print $4,$1,$6,$2,$3}' | sed -e "s/\/[1-2]//g" |sort > readtocontig.txt<br>
samtools view -H alignedmate_GRCh38.sam | cat - <(awk 'NR==FNR{ a[$1]; next }$1 in a{ print $0 ; delete a[$1]; next }' readtocontig.txt <( samtools view alignedmate_GRCh38.sam )) | samtools sort -n -O bam | bedtools bamtobed -i stdin | awk '{OFS="\t"}{print $4,$1,$6,$2,$3}' | sed -e "s/\/[1-2]//g" | sort > pass_mates.txt<br>
join -j 1 readtocontig.txt pass_mates.txt > mates_region.txt<br>

•	Examine links to contig ends only, and filter based on described unambiguity criteria<br> 
Place_region.py<br> 

•	Extracted contig ends and GRCh38 regions with samtools faidx<br> 
samtools faidx GRCh38_no_alt.fa Place_region > GRCh38_Region.fa<br> 

•	Align contigs to the region determined by the linking mates <br>
nucmer  --maxmatch -l 15 -b 1 -c 15 -p alignment_contig GRCh38Regions.fa end_contig.fa<br>
delta-filter -q -r -o 0 -g aliged_info.delta > filtered_info.delta <br> 

•	Obtain BEP/LEP/REP/ contigs and the corresponding placedment positions  <br> 
Contig_type.py <br> 
Please remove contigs that the both end aligned to reference from the LEP/REP file. <br>

# Step3. Cluster placed contigs <br>
•	Cluster placed contigs <br>
Bedtools sort -i  placed_contigs.bed >  placed_contigs.sorted.bed <br>
bedtools merge -d 20 -c 4 -o distinct -i  placed_contigs.sorted.bed > merge_contigs.bed <br>

•	Choose the longest one as the representatives and get the corresponding clusters <br>
Rep_obtain.py <br>

•	Remove contigs with no alignments to representatives <br>
nucmer -p align_info  rep.fa cluster.fa<br>

•	Align other types of contigs to sequences in current clusters <br>
makeblastdb -in remaining_cluster.fa -dbtype nucl -out remainingcontigs_Id<br>
blastn -db remainingcontigs_Id -query othertype_contig.fa -outfmt 6 -max_target_seqs 1  -max_hsps 1  -out  othertype_contig.tsv<br>
a. Contigs can be added to the current cluster if they are fully contained with >99% identity and covered >80% of the aligned cluster. <br>
b. If the coverage of the current cluster is under 80%, we record the ID of the current cluster and the contig ID. (in candidate_contigs.txt) <br>
c. Obatain the contigs that satisy several contiditions. The pass contigs need to be added to the current cluster.<br>
Pass_contigs.py<br>
 
•	Merge left-end placed and right-end placed contigs into a longer insertion<br>
a. If an LEP contig and an REP contig were within 100 bp in the same orientation, please align the two contigs with each other. <br> 
nucmer -f  -p align_info left_placed.fa  right_placed.fa<br>
delta-filter -q  -r -g -m -1 align_info > filterdalign_info.delta<br>
show-coords -H -T -l -c -o filterdalign_info.delta > filterdalign_info.coords<br>
b. CLassfiy the alignment result into four types:
Identity: awk '{OFS="\t"}{if ($NF=="[IDENTITY]") print $0}' filterdalign_info.coords | sort |uniq > Identity.txt<br>
Contained identity>= cutoff : awk '{OFS="\t"}{if ($7>=cutoff && ($NF=="[CONTAINED]" || $NF=="[CONTAINS]")) print $0}' filterdalign_info.coords |sort |uniq  > Contained.txt<br>
Overlap: awk '{OFS="\t"}{if ($7>=90 && $11>=5&& $NF=="[END]") print $0}' filterdalign_info.coords |sort|uniq  > Overlap.txt<br>
Partially map: awk '{OFS="\t"}{if (($10>=coverage_cutoff || $11>=coverage_cutoff) && $NF!="[IDENTITY]" && $NF!="[CONTAINS]" && $NF!="[CONTAINED]") print $0}' filterdalign_info.coords|sort|uniq  > Part.txt<br>

Note: <br>
For the fourth situation, please further check wehther there is at least one contig shared by the two clusters.<br>
        nucmer -p Lrep_Rcluster  REP_cluster.fa LEP_rep.fa   <br>
        nucmer -p Rrep_Lcluster LEP_cluster.fa  REP_rep.fa<br>
        delta-filter  -r -q -g LEP_rep_REP_cluster.delta > LEP_rep_REP_cluster_filter.delta<br>
        delta-filter  -r -q -g REP_rep_LEP_cluster.delta > REP_rep_LEP_cluster_filter.delta<br>
        show-coords -H -T -l -c -o LEP_rep_REP_cluster_filter.delta > LEP_rep_REP_cluster_filter.coords<br>
        show-coords -H -T -l -c -o REP_rep_LEP_cluster_filter.delta > REP_rep_LEP_cluster_filter.coords  <br>      
        
c. Merge pass LEP and REP contigs.
popins merge -c LEP_REP.fa <br>

•	Remove the redundancy of placed contigs<br>
makeblastdb -in all_placed.fa -dbtype nucl -out all_placed_Id<br>
blastn -db all_placed_Id -query all_placed.fa -outfmt 6 -max_target_seqs 1  -max_hsps 1  -out  all_placed_aligned.tsv<br>

•	Cluster the unplaced contigs<br>
cd-hit-est -i remain_unplaced.fa -o unplaced_cluster  -c 0.9 -n 8 <br>

# Additional programs used to analyze CPG  <br>
•	Annotate placed contigs <br>
vep -i contig_insertion_points.vcf -o contig_annotation --dir Cache_path --cache --offline --fasta GRCh38_primary.fa --species homo_sapiens --everything --plugin StructuralVariantOverlap,file=gnomad_v2_sv.sites.vcf.gz <br>

•	Compare with other genomes <br>
bwa index -p other_genome_Id  other_genome.fa <br>
bwa mem other_genome_Id CPG.fa > alignment.sam <br>

•	Call variants  <br>
bwa index -p new_ref_Id  new_ref.fa <br>
bwa mem new_ref_Id read1.fq read2.fq > alignment.sam <br>
java -jar picard.jar MarkDuplicates I=alignment.sam O=alignment.markdup.sam M=alignment.markdup.txt <br>
java  -jar picard.jar BuildBamIndex I=alignment.markdup.sam <br>
gatk HaplotypeCallerSpark -R GRCh38_decoy.fa -I alignment.markdup.sam -O vcffile <br>

