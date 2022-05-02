#! /usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use Cwd;


my $usage=<<USAGE;

Describe:find DNAe site in a sequence... 

 Author: Yun Wang, wangyun\@genomic.org.cn
Version: 1.0, Date: 2012-11-11

Usage:perl $0 <infile> <outfile> [-option]
  
      -DNAe     <str>  : A inducible DNA elements,such as loxp and rox; #Only loxp is available.
 
      -help            : show this message

Example:
  perl $0 <infile> <outfile> -DNAe loxp

USAGE

my ($DNAe,$help);
GetOptions(
   "DNAe:s"    => \$DNAe,  
   "help"    => \$help,
   );

$DNAe ||= "loxp";
   
die $usage if (@ARGV != 2 ||$help);

open IN, "$ARGV[0]" ||die $!;
open OUT, ">$ARGV[1]" ||die $!;

my $chrid = <IN>;
chomp $chrid;
$chrid =~ s/^>(.+)/$1/; 
my $seq;
while(<IN>){
  chomp $_;
  $seq .= $_;
}

$seq = uc ($seq);
my $len = length($seq);

my %Deseq;
$Deseq{"loxp"} = "ATAACTTCGTATAATGTACATTATACGAAGTTAT";

my $offset = 0;

my $result = index($seq, $Deseq{$DNAe}, $offset);
my $star =  $result + 1;
my $END  =  $star + 33;
print OUT "FTR\tloxp\t$chrid\t$star\t$END\n";

while ($result != -1) {
$offset = $result + 1;
$result = index($seq, $Deseq{$DNAe}, $offset);
$star =  $result + 1;
$END  =  $star + 33;
print OUT "FTR\tloxp\t$chrid\t$star\t$END\n" if $star > 0;
}
close IN;
close OUT;


