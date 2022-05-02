#! /usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use Cwd;

my $usage=<<USAGE;

Describe:pick up umapped loxp reads.

  Author: Yun Wang, wangyun\@genomic.org.cn
 Version: 1.0, Date: 2012-11-05

Usage:perl $0 <loxpfq> <loxpsp>  [-option]

      -spid     [int]  : Sample ID/prefix of output file    
      -outdir          : Directory of output [default: ./]   
      -help            : show this message

USAGE

my ($samplename,$outdir,$help);
GetOptions(
   "spid:s"   => \$samplename,
   "outdir:s" => \$outdir,
   "help"   => \$help,
);
die $usage if (@ARGV != 2 ||$help);
	
print "#$samplename\n#pick up umapped loxp reads\n";

my ($loxpfq,$loxpsp) = @ARGV;
open IN1, $loxpfq || die $!;
open IN2, $loxpsp || die $!;

$outdir ||= getcwd;
open OUT, ">$outdir/$samplename.umap.fq" ||die $!;

my %loxpread;
my $Tnum =1;
$/ = "\n@";

while (<IN1>){
  chomp;  
  my $chunkinfo = $_;
  if ($Tnum == 1){
    $chunkinfo =~ s/^@//;
     print "$chunkinfo\n";
  }
  
	my @fqinfo = split /\n/,$chunkinfo;
  if (@fqinfo < 4) {
      my $addinfo = <IN1>;
      chomp $addinfo;
      $addinfo = "\@".$addinfo;
      $chunkinfo .= "\n".$addinfo;
      #print "$chunkinfo\n";
  }
  
  my $readid = $fqinfo[0];
  $loxpread{$readid} = $chunkinfo;
  $Tnum++;
}

print "Total number of loxp reads:$Tnum\n";
close IN1;
$/ = "\n";
my $mapnum=0; 
while (<IN2>){
	chomp;
	my $readid = (split /\s+/,$_)[0];
	if (exists $loxpread{$readid}){
		delete $loxpread{$readid};
		$mapnum++;
   	}
  }
  
  print "Total of loxp mapped reads: $mapnum\n";

my $umapnum = 0;
my $ck = 1;
foreach (keys %loxpread){
	next unless defined ($loxpread{$_});
	print OUT "\@$loxpread{$_}" if($ck ==1);
	print OUT "\n\@$loxpread{$_}";
	$ck ++;
  $umapnum++;
}
print "Total number of umapped reads:$umapnum\n";
close IN2;
close OUT;
