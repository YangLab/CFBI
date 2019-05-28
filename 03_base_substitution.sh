#!/bin/bash
if [ ! $# = 3 ]; then
    echo "Usage: `basename $0` sample_SE_R1.bam[input BAM] target.fa[target sequences] sample_SE_R1[output]";
    echo "";
    echo "   Calculate base substitution with mapped reads (BAM).";
    echo "   eg: sh 03_base_substitution.sh sample_R1.bam target.fa sample_R1";
    echo "";
    echo "Prerequisites: BWA-MEM Version: 0.7.17, SAMtools Version: 1.9";
    echo "";
    exit 0;
fi

`pileupBam_for_target.pl -i $1 -s $2 -o 0 -minBQ 0 -HPB 0 -v 0 -eSignal 0 --cRatio 0 -depth 8000 > $3.vcf`;
`awk '{print $1":"$2}' $3.vcf > $3.sites`;
`awk '{print $1":"$2"\t"$0}' $3.vcf |cut -f 1,4-19 > $3.link`;
`awk -F"\t" '{if(($2=="A")||($2=="a")){printf "%s\t%s\t%d\t%d\t%.8f\t%d\t%.8f\n", $1,$2,$3,$6,($6/$3),($3-$6),(($3-$6)/$3)};if(($2=="C")||($2=="c")){printf "%s\t%s\t%d\t%d\t%.8f\t%d\t%.8f\n", $1,$2,$3,$8,($8/$3),($3-$8),(($3-$8)/$3)};if(($2=="G")||($2=="g")){printf "%s\t%s\t%d\t%d\t%.8f\t%d\t%.8f\n", $1,$2,$3,$10,($10/$3),($3-$10),(($3-$10)/$3)};if(($2=="T")||($2=="t")){printf "%s\t%s\t%d\t%d\t%.8f\t%d\t%.8f\n", $1,$2,$3,$14,($14/$3),($3-$14),(($3-$14)/$3)}}' $3.link > $3.out`
`awk '{if($7>=0)print $1}' $3.out > $3.filter.sites`;
`select_ID.pl $3.link $3.filter.sites 1 |awk '{if($3>=0) print $0}' |cut -f 1-3,6-17 > $3.unique.mut`;
`cut -f 1 $3.unique.mut |awk -F":" '{print $1"\t"$2}' |paste - $3.unique.mut |cut -f 1-15 > $3.unique.out1`;
`join_ID.pl $3.unique.out1 $3.out 3 1 |cut -f 1-15,19-22 > $3.unique.out2`;
`awk -F":" '{print $1"\t"$2-1"\t"$2+1}' $3.filter.sites > $3.unique.agp`;
`get_agp_sequence.pl $2 $3.unique.agp |xargs |sed 's/>/\n/g' |awk -F":" '{print $1"\t"$2}' |awk '{print $1"\t"$2"\t"$3}' |cut -f 2 |awk -F"-"  '{print $1+1}' |sed '1d' > $3.unique.loc`;
`get_agp_sequence.pl $2 $3.unique.agp |xargs |sed 's/>/\n/g' |awk -F":" '{print $1"\t"$2}' |awk '{print $1"\t"$2"\t"$3}' |cut -f 1,3 |sed '1d' > $3.unique.base`;
`paste $3.unique.base $3.unique.loc |awk '{print $1":"$3"\t"$2}' > $3.unique.seq`;
`join_ID.pl $3.unique.out2 $3.unique.seq 3 1 |awk '{OFS="\t"}{print $1,$2,$4,$21,$5,$18,$19,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' |awk '{if($7>=0)print $0}' > $3.xls`;
`rm $3.vcf $3.sites $3.link $3.out $3.filter.sites $3.unique.mut $3.unique.out1 $3.unique.out2 $3.unique.agp $3.unique.loc $3.unique.base $3.unique.seq`;
