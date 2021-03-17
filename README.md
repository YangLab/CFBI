# CFBI
**CFBI** (**C**alculating **F**requencies of **B**ase-substitutions and **I**ndels) is the custome Perl and Shell scripts for calculating frequencies of base subsutitution and indels.
Version: 1.0.0

# Features

Calaulate base substitution and indel frequencies of target gene for CRISPR genome editing.

# BAM format file

BAM file was originally mapped form BWA-MEM.

# Prerequisites

**Software / Package**

* [Perl](https://www.perl.org/) v5.26.2
* [BioPerl](https://bioperl.org/) v1.7.2
* [BWA](http://bio-bwa.sourceforge.net/) Version: 0.7.17
* [SAMtools](http://samtools.sourceforge.net/) Version: 1.9
* [GCC](https://gcc.gnu.org/)
* [GNU coreutils](http://www.gnu.org/licenses/gpl.html)

# Usage: 
-----------------------------------
To calculate base substitution and indel frequencies, BAM format file was generated firstly.

All Perl and Shell scripts were marked ***'bold italic'***.
* 1. Please export the CFBI directory and BioPerl to your **$PATH**.
```bash or zsh
export PATH="~/CFBI-master:$PATH";
export PATH="~/Bio:$PATH";

```

* 2. BWA index target sequences. 
***'01_bwa_mem_index.sh'***.
```bash
sh 01_bwa_mem_index.sh target.fa
```

* 3. BWA-MEM mapping with DNA-seq reads (FASTQ). For paired-end sequencing, only R1 reads were used.
***'02_bwa_mem_mapping.sh'***.
```bash
sh 02_bwa_mem_mapping.sh target.fa sample_R1.fq sample_R1
```

* 4. Calculate base substitution with mapped reads (BAM).
***'03_base_substitution.sh'***.
```bash
sh 03_base_substitution.sh sample_R1.bam target.fa sample_R1
```
Output file [**sample_R1.xls**](https://github.com/YangLab/CFBI/blob/master/example_output/sample_R1.xls) is an example result of base substitution.

* 5. Calculate indel frequencies with target gene InDel location.
***'04_indel_frequencies.sh'***.
```bash
sh 04_indel_frequencies.sh sample_R1.bam EXM1 250 300 ascii.txt
```
Output file [**EMX1_indel_frequencies.txt**](https://github.com/YangLab/CFBI/blob/master/example_output/EMX1_indel_frequencies.txt) is an example for the number of indel frequencies of target EMX1.

-----------------------------------

# Input files
1. Target sequences with FASTQ format. (target.fa)
2. DNA-seq R1 reads with FASTQ format. (sample_R1.fq)
3. Name of BAM file. (sample_R1)
4. Target gene symbol. (EMX1)
5. Start location for indel frequencies calculation of target gene. (250, cutting site -25 bp for EMX1)
5. End location for indel frequencies calculation of target gene. (300, cutting site +25 bp for EMX1)
6. The sequence quality of reads by ASCII code. [ascii.txt](https://github.com/YangLab/CFBI/blob/master/ascii.txt)


# Output files
See details in [sample_R1.xls](https://github.com/YangLab/CFBI/blob/master/example_output/sample_R1.xls), the example output file for base substitution.

| Field       	          | Description                                  |
| :---------------------: | :------------------------------------------: |
| Gene symbol  	          | Name of gene                                 |
| Mut location    	      | Location of mutant site                      |
| Ref base       	        | Sequence of base   	                         |
| Flanking base           | ± 1 bp flanking sequence of base	           |
| # of mapped reads       | Number of mapped reads for target gene       |
| # of total mutant reads | Number of total mutant reads for target gene |
| % of total mutant reads | Ratio of total mutant reads for target gene	 |
| # of mutant A           | Number of reads for mutant A	               |
| % of mutant A           | Ratio of reads for mutant A	                 |
| # of mutant C           | Number of reads for mutant C	               |
| % of mutant C           | Ratio of reads for mutant C	                 |
| # of mutant G           | Number of reads for mutant G	               |
| % of mutant G           | Ratio of reads for mutant G	                 |
| # of mutant N           | Number of reads for mutant N	               |
| % of mutant N           | Ratio of reads for mutant N	                 |
| # of mutant T           | Number of reads for mutant T	               |
| % of mutant T           | Ratio of reads for mutant T	                 |

See details in [EMX1_indel_frequencies.txt](https://github.com/YangLab/CFBI/blob/master/example_output/EMX1_indel_frequencies.txt), the example output file for indel frequencies.

| Field       	          | Description                                  |
| :---------------------: | :------------------------------------------: |
| Gene symbol  	          | Name of gene                                 |
| Type of reads    	      | Type of all reads or indel reads             |
| # of reads       	      | Number of all reads or indel reads   	       |


# Citation

The related paper about CFBI is submitted.

# License

Copyright (C) 2021 YangLab. Licensed GPLv3 for open source use or contact YangLab (yanglab@picb.ac.cn) for commercial use.
