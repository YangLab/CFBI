#!/usr/bin/perl -w
use strict;
if( @ARGV != 2 ) {   
    print "Usage: perl $0 file.fasta all-agp.location \n";
    print "Destination: find all sequence from file.fasta in all-agp.location(column 2 ~ column 3) \n";
    exit 0;
}
use Bio::Seq;
use Bio::SeqIO;
use Bio::PrimarySeq;

my $fh1=shift @ARGV;
my $fh2=shift @ARGV;
my $name=shift @ARGV;

my $in=Bio::SeqIO->new(-file=>"$fh1",'-format'=>'fasta');
my $seq=$in->next_seq();
my $string=$seq->seq();
my $SEQ;
my $disp=$seq->id();
while ($seq)
{
  $disp=$seq->id();
  $string=$seq->seq();
  $SEQ->{$disp}=$string;
  $seq=$in->next_seq();
}
open FH2,"<$fh2";
my $j=1;
while (<FH2>)
{
  chomp($_);
  my @rec=split(/\s+/,$_);
  
  print ">".$rec[0].":".$rec[1]."-".$rec[2]."\n";
  print substr($SEQ->{$rec[0]},($rec[1]-1),($rec[2]-$rec[1]+1))."\n";

}  
close FH2;
