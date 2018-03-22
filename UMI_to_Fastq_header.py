# This script removes the random barcode (first N bases) from each read sequence and appends it to the read name (the Fastq header)
# The first argument: working directory to look for gzipped Fastq files (with fastq.gz or fq.gz extension);
# The second argument: the length of random barcode to trim;

import os, sys, gzip
from itertools import islice

input_dir, barcode_length = sys.argv[1], int(sys.argv[2])

filenames = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f)) and (f.endswith('fastq.gz') or f.endswith('fq.gz'))]
for filename in filenames:
    print(filename)
    suffix = filename[filename.find('.'):] # part of input_filename after the first '.'
    output = filename[:-len(suffix)] + '_UMI' + suffix
    with gzip.open(filename, 'rt') as input_file, gzip.open(output, 'wt') as output_file:
        count = 0
        while True:
            lines = list(islice(input_file, 4)) # Fastq files expected as input
            if not lines:
                break
            count += 1
            read_name, read_seq, read_qual = lines[0].rstrip('\n'), lines[1].rstrip('\n'), lines[3].rstrip('\n')
            barcode = read_seq[:barcode_length]
            read_name = read_name.replace(' ', '_') + ':' + barcode
            read_seq, read_qual = read_seq[barcode_length:], read_qual[barcode_length:]
            output_file.write(read_name + '\n' + read_seq + '\n+\n' + read_qual + '\n')
        print('\t', count, 'reads processed;')
print('Done!')