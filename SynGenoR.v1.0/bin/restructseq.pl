#! /usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use Cwd;

my $usage=<<USAGE;
Describe:Only used to refactor SCRaMbLEd genome sequence.
  
  Author: Wang Yun, wangyun\@genomics.org.cn
 Version: 1.0  Date: 2012-11-05
 
   Usage:perl $0 <ref_fa> <loxpregion> <encod> [option]
         
         -spid     <str>  sample id
         -reseqid  <str>  reconstructed sequence id
         -outdir          output directory
         -help            show this message
         
    Note:This script only get a sequence of the first path in <encod> file.

USAGE
  
my ($sample,$reseqid,$outdir,$help);
GetOptions(
   "spid:s" => \$sample,
   "reseqid:s" => \$reseqid,
   "outdir:s" => \$outdir,
   "help" => \$help,
);


die $usage if (@ARGV != 3 ||$help);

my ($ref_fa,$loxpregion,$encod) = @ARGV;
 
open IN1, $ref_fa ||die $!;
open IN2, $loxpregion ||die $!;
open IN3, $encod || die $!;

$outdir ||= getcwd;
open OUT,">$outdir/$sample.reseq" || die $!;

#print "loading data...\n";

my %lxpregseq;
my $refseq=();
my $n=0;
my $chrid = <IN1>;

while (<IN1>){
	chomp;
	$refseq .= $_;
}
my $len= length($refseq);
print "Reference length: $len\n";

while(<IN2>){
	chomp;
	my ($star,$END,$indx)=(split /\t/)[3,4,8];
	my $reglen = $END - $star +1;
	 
	$lxpregseq{$indx} = substr($refseq,$star-1,$reglen);
  $n++;
 }
print "Total $n loxp region!\n";	

my $restructseq = <IN3>;

unless(-s "$outdir/$sample.encod"){
	 open OUT2,">$outdir/$sample.encod" || die $!;
	 print OUT2 "$restructseq";
}

chomp $restructseq;
my $encode = (split /\s+/,$restructseq)[1];
my @codeord = split /,/,$encode;

my $reseq;	
foreach (@codeord){
	if ($_ == 0){
		next;
	}elsif($_>0){
	   $reseq .= $lxpregseq{$_};
	}elsif($_<0){
	   my $idx = abs($_);
	   my $revcom = 	reverse($lxpregseq{$idx});
	   $revcom =~ tr/ATGCatgc/TACGtacg/;
	 	 $reseq .= $revcom;
  	}
 }

chomp $chrid; 

if($reseqid){ 
	print OUT ">$reseqid\n";
}else{
	$reseqid = "$chrid"."_scb";
  print OUT "$reseqid\n";
}

my $lenscb = length($reseq);
print "Refactoring sequence length:$lenscb\n";
print OUT "$reseq\n";

close IN1;
close IN2;
close IN3;
close OUT;

#===end===	
