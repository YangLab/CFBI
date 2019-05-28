#!/bin/bash
if [ ! $# = 3 ]; then
    echo "Usage: `basename $0` target.fa[target sequences] sample_R1.fq[FASTQ reads] sample_R1[output SAM/BAM]";
    echo "";
    echo "   BWA-MEM mapping with DNA-seq reads (FASTQ). For paired-end sequencing, only R1 reads were used.";
    echo "   eg: sh 02_bwa_mem_mapping.sh target.fa sample_R1.fq sample_R1";
    echo "";
    echo "Prerequisites: BWA-MEM Version: 0.7.17, SAMtools Version: 1.9";
    echo "";
    exit 0;
fi

`bwa mem $1 $2 -t 12 > $3.sam`;
`samtools view -bh -F 4 $3.sam |samtools sort --threads 12 -o $3.bam`;
`samtools index $3.bam`;
`rm $3.sam`;
