# This script removes PCR duplicate reads from BAM file based on UMI barcodes;
# Requires pysam module installed;
# Input: directory containing coordinate-sorted single-end BAM files (UMI sequences should have been appended to read names by the UMI_to_Fastq_header.py script prior to alignment);
# Output BAM is not guaranteed to be sorted in the same order as the input BAM, thus it is recommended to sort the output by Samtools;

progress_chunk = 1000000
buffer_limit = 10000

import os, sys, random
try:
    import pysam
except:
    print('Check if pysam module is installed!')
    exit

def writeReads(reads):
    groups = {}
    n = 0
    for read in reads:
        umi, alignment = read[3], read[4]
        if umi not in groups:
            groups[umi] = [alignment]
        else:
            groups[umi].append(alignment)
    for umi in groups:
        alignments = groups[umi]
        output_file.write(random.choice(alignments))
        n += 1
    return n

def splitReverseBuffer(input_buffer):
    reads_to_write, output_buffer = [], []
    last_read = input_buffer[-1]
    for read in input_buffer:
        chrom, first_base = read[0], read[2]
        if chrom != last_read[0] or first_base < last_read[1]: # if 5' base (rightmost_coord) is less than leftmost coordinate of the last read
            reads_to_write.append(read)
        else:
            output_buffer.append(read)
    return reads_to_write, output_buffer

def TempBuffer(sorted_reads):
    temp_buffer = []
    n = 0
    for sorted_read in sorted_reads:
        if not temp_buffer or (temp_buffer[-1][0] == sorted_read[0] and temp_buffer[-1][2] == sorted_read[2]):
            temp_buffer.append(sorted_read)
        else:
            n += writeReads(temp_buffer)
            temp_buffer = [sorted_read]
    n += writeReads(temp_buffer)
    del temp_buffer
    return n

input_dir = sys.argv[1]
filenames = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f)) and (f.endswith('bam'))]
for filename in filenames:
    print(filename)
    input_file = pysam.AlignmentFile(filename, 'rb')
    alignments = input_file.fetch(until_eof = True)
    output_file = pysam.AlignmentFile(filename[:-4] + '_dedup.bam', 'wb', template = input_file)
    progress, written = 0, 0
    reads_fw, reads_rev = [], []
    for alignment in input_file:
        progress += 1
        if progress % progress_chunk == 0:
            print('.', end = '', flush = True)
        umi = alignment.query_name.split(':')[-1]
        chrom = alignment.reference_name
        leftmost_coord = alignment.reference_start
        rightmost_coord = alignment.reference_end
        read = (chrom, leftmost_coord, rightmost_coord, umi, alignment)
        if alignment.is_reverse == False:
            if not reads_fw or (reads_fw[-1][0] == read[0] and reads_fw[-1][1] == read[1]):
                reads_fw.append(read)
            else:
                written += writeReads(reads_fw)
                reads_fw = [read]
        elif alignment.is_reverse == True:
            reads_rev.append(read)
            if len(reads_rev) == buffer_limit:
                old_reads, reads_rev = splitReverseBuffer(reads_rev)
                old_reads = sorted(old_reads, key = lambda x : x[2]) # sort by 5' base (rightmost_coord)
                written += TempBuffer(old_reads)
    written += writeReads(reads_fw)
    reads_rev = sorted(reads_rev, key = lambda x : x[2])
    written += TempBuffer(reads_rev)
    input_file.close()
    output_file.close()
    del reads_fw, reads_rev
    if progress > progress_chunk:
        print('\n')
    print(progress, 'input reads,', written, 'output_reads,', round((progress - written)/progress*100, 1), '% duplicates;')
print('Done!')