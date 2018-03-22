# DETAILS ABOUT THE LIB PREP PROTOCOL!

# Quality and adapter trimming (observe the custom adapter sequence):
for file in *fastq.gz; do echo $file && trim_galore --adapter "ATCTCGTATGCCG" $file; done

# Trim UMIs (first 8 nt of each read), append them to Fastq headers:
python3 UMI_to_Fastq_header.py . 8

# Align to TAIR10 using STAR:
for file in ./*UMI.fq.gz; do echo $file && STAR --genomeDir /index/tair10/star --readFilesIn $file --runThreadN 4 --outFileNamePrefix ${file/.fq.gz/_} --outSAMmultNmax 1 --alignEndsType EndToEnd --readFilesCommand zcat; done

# Sort SAM files and convert to BAM:
for file in *sam; do echo $file && samtools view -hu $file | samtools sort - -o ${file/.sam/_sorted.bam} && rm $file; done

# Filter out rRNA, tRNA and sn/snoRNA reads:
for file in *sorted.bam; do echo $file && bedtools intersect -v -abam $file -b Araport11_rRNA_tRNA_snRNA_snoRNA.bed > ${file/.bam/_filt.bam}; done

# Filter out multimapper reads:
for file in *filt.bam; do echo $file && samtools view -h -q 10 $file -o ${file/.bam/_mapq10.bam}; done

# Deduplicate on UMIs:
python3 Deduplicate_BAM_files_on_UMI.py .

# Sort again by coordinates:
for file in *dedup.bam; do echo $file && samtools sort $file -o ${file/.bam/_sorted.bam}; done

# Make stranded Bedgraph files (only the first base of each read is considered):
for str in "+" "-"; do
  echo $str
  [ "$str" = "+" ] && n="fw" || n="rev"
  for file in *dedup_sorted.bam; do
    echo $file && bedtools genomecov -ibam $file -bg -5 -strand $str | sort -k 1,1 -k 2,2n > ${file/.bam/}_${n}.bg
  done
done

# Remove singletons (positions covered by a single tag):
for file in *bg.gz; do echo $file && zcat $file | awk '{if ($4>1) print $0}' > ${file/.bg.gz/_cov2.bg}; done

# Ensure that all positions have 1 bp width:
python3 Expand_bedGraph_to_single_base_resolution.py .

# Convert Bedgraph files to Bigwig (requires kentUtils):
for file in *cov2.bg; do echo $file && bedGraphToBigWig $file TAIR10.chrom.sizes ${file/bg/bw}; done