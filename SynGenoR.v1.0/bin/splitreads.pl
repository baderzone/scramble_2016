#! /usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use Cwd;

my $usage=<<USAGE;

Describe:split umapped loxp reads...

  Author: Yun Wang, wangyun\@genomic.org.cn
 Version: 1.0, Date: 2012-11-05

Usage:perl $0 <loxpfq> <loxpsp>  [-option]

      -spid     <str>  : Sample ID/prefix of output file  
      -minlen   [int]  : minimal lengnth of split reads 
      -method   [str]  : left or remove the loxpsym site
      -DNAe     <str>  : A inducible DNA elements,such as loxp and rox; #Only loxp is available.
      -outdir          : Directory of output [default: ./]   
      -help            : show this message

USAGE

my ($samplename,$len,$method,$DNAe,$outdir,$help);
GetOptions(
   "minlen:s" => \$len, #output the reads which is more than L;  
   "method:s"  => \$method, #left or remove the loxpsym site
   "spid:s" => \$samplename,
   "DNAe:s"    => \$DNAe, #Only loxp is available. 
   "outdir:s" => \$outdir,
   "help" => \$help,
);
die $usage if (@ARGV != 1 ||$help);

print "#$samplename\n";
print "#split ummapped loxp reads\n";

my $loxpumfq = $ARGV[0];

open IN,$loxpumfq ||die $!;

$outdir ||= getcwd;
$len ||= 15;
$method ||= "left";
my %Deseq;
$Deseq{"loxp"} = "ATAACTTCGTATAATGTACATTATACGAAGTTAT";

open OUT, ">$outdir/$samplename.splitted.fq" ||die $!;

$/ = "\n@";

my $num =1;
while (<IN>){
  chomp;
  my $chunkinfo = $_;
  if ($num == 1){
    $chunkinfo =~ s/^@//;
  }
  
	my @fqinfo = split /\n/,$chunkinfo;
  if (@fqinfo < 4) {
      my $addinfo = <IN>;
      chomp $addinfo;
      $addinfo = "\@"."$addinfo";
      $chunkinfo .= "\n$addinfo";
  }
  splitfq($chunkinfo);
  $num++;
}
close IN;
close OUT;
print "Reads number:$num\n";

sub splitfq{
	my $chunk  = $_[0];
	my @chunkline = split /\n/,$chunk;
  my $site = index($chunkline[1],$Deseq{"loxp"});
  #my $rdslen = length ($chunkline[1]);
  my ($read_1,$read_2,$qual_1,$qual_2);
  if($method =~ /remove/i){
    $read_1 = substr ($chunkline[1],0,$site);
    $read_2 = substr ($chunkline[1],$site+34);
    $qual_1 = substr ($chunkline[3],0,$site);
    $qual_2 = substr ($chunkline[3],$site+34);
  }elsif($method =~/left/i){	
    $read_1 = substr ($chunkline[1],0,$site + 34);
    $read_2 = substr ($chunkline[1],$site);
    $qual_1 = substr ($chunkline[3],0,$site + 34);
    $qual_2 = substr ($chunkline[3],$site); 
  }
  my $readid_1 = $chunkline[0]."_1";
  my $readid_2 = $chunkline[0]."_2";
  
  if (length($read_1)>=$len && length($read_2)>=$len){
    print OUT "\@$readid_1\n$read_1\n+\n$qual_1\n";
    print OUT "\@$readid_2\n$read_2\n+\n$qual_2\n";
  }
}
 
 close OUT;
