#!/bin/bash
if [ ! $# = 5 ]; then
    echo "Usage: `basename $0` sample_SE_R1.bam[input BAM] Gene_ID[Gene symbol] target_Start[target start location] target_End[target end location] ascii.txt[ASCII file]";
    echo "";
    echo "   Calculate indel frequencies with target gene InDel location.";
    echo "   eg: sh 04_indel_frequencies.sh sample_R1.bam EXM1 250 300 ascii.txt";
    echo "";
    echo "Prerequisites: BWA-MEM Version: 0.7.17, SAMtools Version: 1.9";
    echo "";
    exit 0;
fi

`samtools view $1 $2:$3-$4 > 01_all.sam`;
`cut -f 1 01_all.sam |sort -u > $2_all_reads.txt`;
`awk '{if($6~"D" || $6~"I") print $0}' 01_all.sam > 02_indel.sam`;
`cut -f 1 02_indel.sam |sort -u > 02_indel.rid`;
`perl -alne '$"="\t";$n=0;for($F[5]=~/(\d+\D)/g){$s=$_;$s=~/(\d+)/;$m=$1;if($s=~/I/){for(1..$m){$seq=substr($F[9],$n,1);$qua=substr($F[10],$n,1);$n++;print "@F[0,2,3,5]\t$n\t$seq\t$qua\t$m"}}elsif($s=~/D/){for(1..$m){$num=$n+1;print "@F[0,2,3,5]\t$num\tseq\tqua\t$m"}};$n+=$m if $s=~/S|M/}' 02_indel.sam |awk '{if($7=="qua"){print $0"\t"$3+$5"\t"$3+$5+$8}else{print $0"\t"$3+$5"\t"$3+$5+1}}' > 03_indel.loc`;
`awk '{if((('$4'-$9)>=1)&&(($10-'$3')>=1))print $0}' 03_indel.loc > 04_indel_region.loc`;
`cut -f 1 04_indel_region.loc |sort |uniq -c |awk '{if($1<=0)print $2}' > 05_remove.rid`;
`select_v_ID.pl 04_indel_region.loc 05_remove.rid 1 > 06_retained.loc`;
`awk '{if($5>=6 && $5<=145)print $0}' 04_indel_region.loc > 07_overhang.loc`;
`cut -f 1 07_overhang.loc |sort -u > 07_overhang.rid`;
`join_ID.pl 07_overhang.loc $5 7 2 |awk '{if($13>=0)print $1}' |sort -u > 08_I_Q0.rid`;
`join_ID.pl 07_overhang.loc $5 7 2 |awk '{if($13<0)print $1}' |sort -u > 08_delete_I_Q0.rid`;
`cat 08_I_Q0.rid 08_delete_I_Q0.rid |sort |uniq -d > 08_retained.rid`;
`cat 08_delete_I_Q0.rid 08_retained.rid |sort |uniq -u > 08_remove.rid`;
`cat 07_overhang.rid 08_remove.rid |sort |uniq -u > $2_indel_reads.txt`;
`wc -l $2_all_reads.txt |awk '{print "'$2'\tall_reads\t"$1}' > $2_indel_frequencies.txt`;
`wc -l $2_indel_reads.txt |awk '{print "'$2'\tindel_reads\t"$1}' >> $2_indel_frequencies.txt`;
`rm 01_all.sam 02_indel.sam 02_indel.rid 03_indel.loc 04_indel_region.loc 05_remove.rid 06_retained.loc 07_overhang.loc 07_overhang.rid 08_I_Q0.rid 08_delete_I_Q0.rid 08_retained.rid 08_remove.rid $2_all_reads.txt $2_indel_reads.txt`;
