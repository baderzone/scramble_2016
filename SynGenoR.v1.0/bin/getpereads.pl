#! /usr/bin/perl -w

use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use PerlIO::gzip;
use Cwd;

my $usage=<<USAGE;

Describe: getting pair-end reads from SOAP alignment results...

  Author: Yun Wang, wangyun\@genomic.org.cn
 Version: 1.0, Date: 2012-11-11

Usage:perl $0 <soap_result>   [-option]
          
      -spid       <str>  : Sample ID/prefix of output file
      -insertsize [int]  : insert size [default: 500]
      -refseq     <str>  : reference of SOAP alignment[*.fa]
      -dejoint           : remove joint maping if synthetic genome is circular
      -distance   [int]  : remove long distance mapping reads if multiple hits [default: 700]   
      -outfile    <str>  : output file
      -outdir            : Directory of output [default: ./]   
      -help              : show this message

example:
     1. not remove joint maping
     perl $0 <soap_result> -spid <str> -distance [int] -outdir
     2. remove joint maping
     perl $0 <soap_result> -spid <str> -dejoint -insertsize [int] -refseq <str> -distance [int] -outdir
     
USAGE

my ($sample,$insertsize,$refseq,$dejoint,$dist,$outdir,$outfile,$help);
GetOptions(  
   "spid:s" => \$sample,
   "insertsize:s" => \$insertsize,
   "refseq:s" => \$refseq,
   "dejoint" => \$dejoint,
   "distance:s" => \$dist,
   "outdir:s" => \$outdir,
   "outfile:s" => \$outfile,
   "help" => \$help,
);

die $usage if (@ARGV != 1 ||$help);
die $usage if ($dejoint && !$refseq);
	
my $reflen = 0;
my $tarchr;
if ($dejoint){
	my @line = (split /\n+/,`cat $refseq`);
	$tarchr = shift @line;
	$tarchr =~ s/>(.+)/$1/; 
	my $seqence .= join '',@line;
	$reflen = length($seqence);
 }

#print "reflen:$reflen..\n";

my $spfile = $ARGV[0];

if ($spfile =~ /.gz$/){
	open IN,"gzip -dc $spfile |" || die $!;
}else{
  open IN,$spfile||die $!;
}

$insertsize ||= 500;
$dist ||= 700;
$outdir ||= getcwd;
$outfile ||="$outdir/$sample.spr";
open OUT,">$outfile" ||die $!;

my $fline;
while(1){
	$fline = <IN>;
  if($fline =~ /^#/){
  	next;
   }else{last;}
}

chomp $fline;
my ($readid,$hitn1,$ori,$chr1,$site1) = (split /\s+/,$fline)[0,3,6,7,8];
my $prit1 = "$readid\t$chr1\t$site1\t$ori";
my $read1;
if ($readid =~ /(.+)\/\d+$/){$read1 = $1;}
my $Tnum = 1; #Total number of reads
my $Tdejnt = 0;#Total number of reads mapping on the joint
my $penum = 0; #Total number of pair reads
my $uqpenum =0; #Total number of unipue mapping pair reads 
my $scalarsize = $insertsize*1.4;

while (<IN>){
	next if($_ =~ /^#/);
	chomp;
	$Tnum++;
	my $hitn;
	my $site;
	my $chr2;
	($readid,$hitn,$ori,$chr2,$site) = (split /\s+/,$_)[0,3,6,7,8];
	my $prit2 = "$readid\t$chr2\t$site\t$ori";
	
	if ($readid =~ /(.+)\/\d+$/){$readid = $1;}
	$penum ++ if ($readid eq $read1);

	#remove joint maping if synthetic genome is circular
	if($dejoint){
		if($readid eq $read1 &&$chr2 eq $chr1 && $chr2 eq $tarchr){
		  
		  if(($site1 < $scalarsize && $site > $reflen - $scalarsize)||($site < $scalarsize && $site1 > $reflen - $scalarsize)){
		  	$Tdejnt ++;
		  	
		  	$prit1 = $prit2;
        $read1 = $readid;
        $hitn1 = $hitn;
        $site1 = $site;
        $chr1 = $chr2;
		  	next;
		  }
		}else{
			  $prit1 = $prit2;
        $read1 = $readid;
        $hitn1 = $hitn;
        $site1 = $site;
        $chr1 = $chr2;
			  next;}
	 }
	
	#remove long distance mapping reads if multiple hits 
	if ($readid eq $read1 && $hitn==1&&$hitn1==1&&$chr2 eq $chr1){
    print OUT "$prit1\t$prit2\n";
    $uqpenum++;
  #}elsif($readid eq $read1 &&$chr2 eq $chr1&& abs($site1-$site)<500){ # 500 was used for chrIXR vision1.0..
   }elsif($readid eq $read1 &&$chr2 eq $chr1&& abs($site1-$site) < $dist){ # 700 was used for loxp reads check..
   	print OUT "$prit1\t$prit2\n";
  }
  
   $prit1 = $prit2;
   $read1 = $readid;
   $hitn1 = $hitn;
   $site1 = $site;
   $chr1 = $chr2;
 }
 
 print "Total number of reads: $Tnum\n";
 print "Total number of pair reads: $penum\n";
 print "Total number of unipue mapping pair reads:$uqpenum\n";
 print "Total number of reads mapping on the joint:$Tdejnt\n" if ($dejoint);
 close IN;
 close OUT;
   
