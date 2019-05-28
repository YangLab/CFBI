#!/usr/bin/env perl
#
# Name:        pileupBam_for_target.pl
# License:     GPL
#
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Bio::DB::Sam;

my $version='1';
my $command=join"\t",@ARGV;
warn "Command:\t$0\t$command\n"; 
local $|=1;
sub usage {
    print <<"END_USAGE";
Usage: perl $0
        --input         FILENAME for bam file
        --seq           Genome sequence
        --overhang      delete the bases on the ends(default >=1: all)
        --reads         more than n reads(default >=1)
        --HPB           more than n HPB(default >5)
        --variants      more than n variants(default >=2)
        --eSignal       effective signal(default >0.95)
        --minBQ         minimum base quality to filter reads ( default >=0,25=ord(:)-33,s//N/)
        --cRatio        change ratio(default >0.05)
        --depth         max per-BAM depth to avoid excessive memory usage (default 8000)
        --all           No cutoff(overhang=1,HPB=0,eSignal=0,cRatio=0,var=0,minBQ=0), 
                        suppress all other parameters.
        eg: 
            pileupBam_for_target.pl -i 01_sample_SE_R1.bam -s target.fa -o 0 -minBQ 0 -HPB 0 -v 0 -eSignal 0 --cRatio 0 -depth 8000 > 02_sample_SE_R1.vcf
        --help
END_USAGE
    exit;
}
####################### get options to overlap default values;
my ($input,$seq,$region,$list,$num,$over,$reads,$HPB,$variants,$eSignal,$minBQ,$cRatio,$depth,$all,$help,$debug);
$over=0; $HPB=0; $variants=0; $eSignal=0;$minBQ=0; $cRatio=0;$reads=0;$depth=8000;
GetOptions (
        'input=s'		=>	\$input,
        'seq=s'         =>  \$seq,
        'overhang=i'    =>  \$over,
        'HPB=f'         =>  \$HPB,
        'reads=i'       =>  \$reads,
        'variants=i'    =>  \$variants,
        'eSignal=f'     =>  \$eSignal,
        'minBQ=i'       =>  \$minBQ,
        'cRatio=f'      =>  \$cRatio,
        'depth=i'       =>  \$depth,
        'help'			=>	\$help,
        ) or usage();
usage() if $help or !$input or !$seq;
if($all){
    $over=0; $HPB=0; $variants=0; $eSignal=0; $cRatio=0;$reads=0;
}
my $paras="";
$paras.="-r $region " if $region;
$paras.="-l $list " if $list;

my %idx=(A=>7,C=>9,G=>11,N=>13,T=>15);
{###### close
    my $count=0;
    open my $in, "samtools idxstats $input |" or die;
    while(<$in>){
        my @F=split /\t/;
        $count+=$F[2];
    }
    #### get rlength
    open $in,"samtools view $input |" or die;
    my $rlength=0;
    while(<$in>){
        my @F=split /\t/;
        if($F[5]=~/^[\dM]*$/){
            $rlength=scalar(split //,$F[9]);
            last;
        }
    }
    my $ratio=1000000000/$count/($rlength-2*$over);
    open my $pileup,"samtools mpileup -Q 0 -B -d $depth -O -f $seq $paras $input|" or die;
    while(<$pileup>){
        my @F=split /\t/;
        $F[4]=~s/(\^.)+|(\$)+//g;
        my $temp=$F[4];
        while($temp=~/(\+|-)(\d+)/g){
            my $ss=$1;
            $ss="\\$ss" if $ss eq "+";
            $F[4]=~s/$ss$2 @{["." x $2]}//x;}
        my @a=split //,$F[4];
        my @positions=map{chomp;$_} split /,/,$F[6]; 
        my @scores=split //,$F[5];
        my @p_sites=grep {$a[$_]!~/[><\*]/ and $positions[$_]>=$over and $positions[$_]<=$rlength-$over+1 } 0..$#positions;   ##### >4 and <96, not only junction
        my $ct_total=@p_sites;
        if($minBQ){
            @p_sites=grep {ord($scores[$_])-33>=$minBQ} @p_sites;
        }
        my $ct=@p_sites; 
        next if $ct<1;
        my @a_new=@a[@p_sites];
        my @positions_new=@positions[@p_sites];
        my @scores_new=@scores[@p_sites];
        my (%ct,$vct,@locs,@ids); $vct=0;
        map{$ct{$_}=0} qw(A C G N T);
        map{$ct{$a_new[$_]=~/[,.]/?uc($F[2]):($vct++,(push @locs,$_),uc $a_new[$_])}++} 0..$#a_new;
        ##### have variants (total variants > 0)
        @ids=sort {$ct{$b}<=>$ct{$a}} grep {$_ ne uc($F[2])} qw(A C G N T);      #### get max variant
        print join("\t",@F[0..2],$ct,join("",@a_new),join("",@a_new[@locs])."^///^".join("",@scores_new[@locs])."^///^".join(",",@positions_new[@locs]),(map {$_?$_:0} map { $ct{$_},$ct{$_}/$ct} qw(A C G N T)),$ct_total,$ct_total*$ratio)."\n";  
    }
}

sub getFiles{
    my $input=$_[0];
    my (@files,%files);
    for (split /,/,$input){
        if(/\|/){
            my $bas=basename($_);
            my $dir=dirname($_);
            for my $f(split /\|/,$bas){
                map {push @files,$_ if ! $files{$_}++} glob "$dir/$f";
            }
        }else{ map {push @files,$_ if ! $files{$_}++} glob $_; }
    }
    return \@files;
}
