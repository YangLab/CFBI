#!/bin/bash
if [ ! $# = 1 ]; then
    echo "Usage: `basename $0` target.fa[target sequences]";
    echo "";
    echo "   BWA index target sequences.";
    echo "   eg: sh 01_bwa_mem_index.sh target.fa";
    echo "";
    echo "Prerequisites: BWA-MEM Version: 0.7.17, SAMtools Version: 1.9";
    echo "";
    exit 0;
fi

`bwa index -a is $1`;
`samtools faidx $1`;
