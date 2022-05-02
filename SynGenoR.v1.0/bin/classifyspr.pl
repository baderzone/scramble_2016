#! /usr/bin/perl -w
use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use List::Util qw/max min/;
use Cwd;


my $usage=<<USAGE;

Describe: classifying spr reads...

  Author: Yun Wang, wangyun\@genomic.org.cn
 Version: 1.0, Date: 2012-11-11

Usage:perl $0  <mpi> <sort_spr>  [-option]
          
      -spid       <str>  : Sample ID/prefix of output file  
      -insertsize <str>  : insert size
      -outdir            : Directory of output [default: ./]   
      -help              : show this message

USAGE

my ($sample,$insertsize,$outdir,$help);
GetOptions(
  
   "spid:s" => \$sample,
   "insertsize:s" => \$insertsize,
   "outdir:s" => \$outdir,
   "help" => \$help,
);

$insertsize ||= 500;
die $usage if (@ARGV != 2 ||$help);
 
 my ($mpi,$sortspr) =  @ARGV;
 
 open MPI, $mpi||die $!;
 open SPR, $sortspr||die $!;
 
$outdir ||= getcwd;
open OUT, ">$outdir/$sample.varid" ||die $!;

print OUT "#\"0\" represent splited reads id\n";
print OUT "#\"1\" represent codirectional mapping pair-end reads id\n";
print OUT "#\"2\" represent the pair-end reads id with mapping distance >= $insertsize\n"; 

while(<MPI>){
	chomp;
	my $id = $_;
	if ($id =~ /(.+)\/\d+$/){$id = $1;}
#	print "$1\n";
	
	print OUT "$id\t0\n";
 }
 
 while(<SPR>){
 	chomp;
 	my ($id,$site1,$ori1,$site2,$ori2) = (split /\s+/,$_)[0,2,3,6,7];
	my $min = min($site1,$site2);
	my $max = max($site1,$site2);
	if ($id =~ /(.+)\/\d+$/){$id = $1;}
  if($ori1 eq $ori2){
	   print OUT "$id\t1\n";

	}elsif(($max-$min)>= $insertsize){
		 print OUT "$id\t2\n";
	}else{next;}
}
close MPI;
close SPR;
close OUT;
