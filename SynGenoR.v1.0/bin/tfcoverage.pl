#! /usr/bin/perl -w

use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use List::Util qw/max min/;
use Cwd;

my $usage=<<USAGE;

Describe: tansforming coverage.....

  Author: Yun Wang, wangyun\@genomic.org.cn
 Version: 1.0, Date: 2012-11-11

Usage:perl $0 <single cvg> [-option]

      -spid     <str>  : Sample ID/prefix of output file
      -wind     [int]  : length of a window [default: 10]
      -chrid    [str]  : chromosome id
      -outdir          : Directory of output [default: ./]   
      -help            : show this message

USAGE

my ($sample,$wind,$chrid,$outdir,$help);
 
GetOptions(
   "spid:s" => \$sample,
   "wind:s" => \$wind,
   "chrid:s" => \$chrid,
   "outdir:s" => \$outdir,
   "help" => \$help,
);

die $usage if (@ARGV != 1 ||$help);

my $cvg = $ARGV[0];
$outdir ||= getcwd;
$wind ||= 10;

open IN,$cvg || die $!;
open OUT,">$outdir/$sample.depth" ||die $!;

#=====Finding max value====##
$/= ">";

<IN>;
my @line;

while(<IN>){
	chomp;	
	@line = split /\n/;
	my $chrID = shift @line;
	last if $chrID eq "$chrid";
	}

my @sgdepth;
foreach (@line){
	my @lnsgdep = split /\s+/; 
	push @sgdepth, @lnsgdep;
}
my $numbase = @sgdepth; 

my $maxdep = max (@sgdepth);
print "$numbase\t$maxdep\n";

#my $site = 1;

#foreach (@sgdepth){
#	my $depR = $_/$maxdep;
#	print OUT "sc9\t$site\t$site\t";
#	printf OUT "%5.4f",$depR;
#	print OUT  "\n";
#	$site ++;
#}


my $num = @sgdepth;
my $depth = 0;
my $cycle = int($num/$wind);
foreach (1..$cycle){
   my $n = $_;
   foreach (1..$wind){
     $depth +=	$sgdepth[($n-1)*$wind+$_-1];
   }
   my $depR = $depth/($wind*$maxdep);
   
   my $site1 = ($n-1)*$wind+1;
   my $site2 = $site1+$wind-1;
   print OUT "$chrid\t$site1\t$site2\t";
   printf OUT "%5.4f",$depR;
   print OUT "\n";
   $depth = 0;
}
my $lstnum = $num -$cycle*$wind;
if($lstnum){

foreach (1..$lstnum){
  $depth +=	$sgdepth[$cycle*$wind+$_-1];
}
my $depR = $depth/($lstnum*$maxdep);
my $site = $cycle*$wind+1;  
print OUT "$chrid\t$site\t$num\t";
printf OUT "%5.4f",$depR;
print OUT "\n"; 
}

close IN;
close OUT;

