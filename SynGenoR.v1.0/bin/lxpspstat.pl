#! /usr/bin/perl -w
use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use List::Util qw/max min/;
use Cwd;

my $usage=<<USAGE;

Describe: Only used to stat the mapping reads with loxp...

  Author: Yun Wang, wangyun\@genomic.org.cn
 Version: 1.0, Date: 2012-11-11

Usage:perl $0 <soap> <lxp.ftr> -prefix <str> -o <output_path> 
          
      -prefix     <str>  : Sample ID/prefix of output file
      -DNAe       <str>  : A inducible DNA elements,such as loxp and rox; #Only loxp is available.  
      -readlen    <str>  : length of each reads
      -splxpid           : extract loxp reads id 
      -outdir            : Directory of output [default: ./]   
      -help              : show this message
     
USAGE

my ($prefix,$DNAe,$readlen,$splxpid,$outdir,$help);
GetOptions(
   "prefix:s" => \$prefix,
   "DNAe:s"  => \$DNAe,
   "readlen:s" => \$readlen,
   "splxpid" => \$splxpid,
   "outdir:s"  => \$outdir,
   "help"  => \$help,
);

die $usage if(@ARGV != 2 ||$help);

$prefix ||= "lxpspstat";
$outdir ||= getcwd;
$DNAe ||= "loxp";
$readlen ||= 100;

my %Deseq;
$Deseq{"loxp"} = "ATAACTTCGTATAATGTACATTATACGAAGTTAT";

open SP, "$ARGV[0]" || die $!;
open FTR, "$ARGV[1]" || die $!;
open OUT, ">$outdir/$prefix.lxps" || die $!;
(open OUT2, ">$outdir/$prefix.loxp.spid" || die $!) if $splxpid; 
##print headline##

print OUT "# READ ME\n# The output format is as flllowing: \n";
print OUT "#chrID\tstart\tend\ttotal\t+mapping\t-mapping\n";

my %feature;

while (<FTR>){
	chomp;
	my ($chrid,$star,$ends) = (split /\s+/,$_)[2,3,4];
	my $mapstar = $ends - $readlen;
	my $mapend = $star;
	push @{$feature{$star}},($chrid,$star,$ends,$mapstar,$mapend);
}

my %lxpspstat = ();
my %positive = ();
my %negtive = ();
foreach (keys %feature){
	$lxpspstat{$_} = 0;
	$positive{$_} = 0;
	$negtive{$_} = 0;
}

while (<SP>) {
	chomp;
	next if ($_=~ /^#/);
	my $ck = 0; 
	my ($readid,$seq,$ori,$chrID,$site) = (split /\s+/,$_)[0,1,6,7,8];
	if($seq =~/$Deseq{$DNAe}/){
		 foreach my $key(keys %feature){
			  if($site >= ${$feature{$key}}[3] && $site <= ${$feature{$key}}[4]){
			    $lxpspstat{${$feature{$key}}[1]}++;
	
			    $positive{${$feature{$key}}[1]} ++ if $ori eq "+";
			    $negtive{${$feature{$key}}[1]} ++ if $ori eq "-";
			    $ck = 1;
			    last;
			   }
		 }
		 if($ck == 0){
		 	print "A read with loxp ,which mapped on $site of $chrID, is NOT counted\n";
		}
  print OUT2 "$readid\n" if ($splxpid); 
 }
}

my $Tlxpsp  = 0;

foreach (keys %lxpspstat){
 	   $Tlxpsp +=  $lxpspstat{$_};
 	   
 	}

print "Total number of $Tlxpsp reads with loxp map on reference\n";

#foreach ( keys %lxpspstat){
foreach ( sort { $a <=> $b } keys %lxpspstat){
	print OUT "${$feature{$_}}[0]\t${$feature{$_}}[1]\t${$feature{$_}}[2]\t$lxpspstat{$_}\t$positive{$_}\t$negtive{$_}\n";
 }

close OUT;

		   
