#! /usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use Cwd;

my $usage=<<USAGE;

Describe:Statistic for basic reads... 

 Author: Yun Wang, wangyun\@genomic.org.cn
Version: 1.0, Date: 2012-11-05

Usage:perl $0 <cleanfq_1> <cleanfq_2> <soap> <soap2> [-option]

      -spid     <str>  : Sample ID/prefix of output file    
      -syntype  <str>  : w|L|R; w,whole; L,left arm; R,right arm  [default: R]   
      -chrsyn   <str>  : Synthetic chromosome id 
      -synlen   [int]  : Length of synthetic chromsome
      -chrwt    <str>  : Reference chromosome id
      -wtlen    [int]  : Length of synthetic reference chromsome
      -reflen   [int]  : Length of whole reference chromsome;
      -DNAe     <str>  : A inducible DNA elements,such as loxp and rox; #Only loxp is available.
      -outdir          : Directory of output [default: ./]   
      -help            : show this message

USAGE

my ($sample,$syntype,$chrsyn,$synlen,$chrwt,$wtlen,$reflen,$DNAe,$outdir,$help);
GetOptions(
   "spid:s"    => \$sample,
   "syntype:s" => \$syntype,
   "chrsyn:s"  => \$chrsyn,
   "synlen:s"  => \$synlen,
   "chrwt:s"   => \$chrwt,
   "wtlen:s"   => \$wtlen,
   "reflen:s"  => \$reflen,
   "DNAe:s"    => \$DNAe,   
   "outdir:s"  => \$outdir,
   "help"    => \$help,
);

die $usage if (@ARGV != 4 ||$help);
die $usage if (!$sample||!$chrsyn||!$synlen||!$chrwt||!$wtlen||!$reflen||!$DNAe);
$syntype ||= "R";
$outdir ||= getcwd;

my ($cleanfq_1,$cleanfq_2,$soap,$soap2) = @ARGV;
my %Deseq;
$Deseq{"loxp"} = "ATAACTTCGTATAATGTACATTATACGAAGTTAT";

if ($cleanfq_1 =~ /\.gz$/) {
  open IN1, "gunzip -c $cleanfq_1 |" || die $!;
  open IN2, "gunzip -c $cleanfq_2 |" || die $!;
}else{
  open IN1, $cleanfq_1 || die $!;
  open IN2, $cleanfq_2 || die $!;
  }
if($soap =~ /\.gz$/){
  open IN3, "gunzip -c $soap |"|| die $!;
}else{
  open IN3, $soap || die $!;
}

if ($soap2 =~ /\.gz$/){
  open IN4, "gunzip -c $soap2 |" || die $!;
}else{
  open IN4, $soap2 || die $!;
}

open OUT1, ">$outdir/$sample.loxp.fq" || die $!;
open OUT2, ">$outdir/$sample.stat" || die $!;
open OUT3, ">$outdir/$sample.loxp.sp" || die $!;

my $Tnum1 = 0;
my $TloxpR1 = 0;
my $Tnum2 = 0;
my $TloxpR2 = 0;
my $Tnum = 0;
my $TloxpR = 0;
my $TloxpN = 0; 

$/ = "\n@";
print "#$sample\n";
while (<IN1>) {
  chomp;
  $Tnum1 ++;
  my $fqinfo = $_;
  my @fqinfo = split /\n/, $fqinfo;
  if (@fqinfo < 4) {
      my $addinfo = <IN1>;
      chomp $addinfo;
      $addinfo = "\@".$addinfo;
      my @fqaddinfo = split /\n/,$addinfo;
      push @fqinfo, @fqaddinfo;
      $fqinfo .= "\n$addinfo";
  }
   
  my $readseq = $fqinfo[1]; 
  if ($readseq =~ /$Deseq{"loxp"}/gi){
  	$Tnum1 == 1 ? print OUT1 "$fqinfo\n": print OUT1 "\@$fqinfo\n"; 	
  	$TloxpR1 ++;
   }
}
close IN1;

while (<IN2>) {
  chomp;
  $Tnum2 ++;
   my $fqinfo = $_;
   my @fqinfo = split /\n/,$fqinfo;
   
  if (@fqinfo < 4) {
      my $addinfo = <IN2>;
      chomp $addinfo;
      $addinfo = "\@".$addinfo;
      my @fqaddinfo = split /\n/,$addinfo;
      push @fqinfo, @fqaddinfo;
      $fqinfo .= "\n$addinfo";
  }
   
  my $readseq = $fqinfo[1]; 
  if ($readseq =~ /$Deseq{"loxp"}/gi){
  	$Tnum2 == 1 ? print OUT1 "$fqinfo\n":print OUT1 "\@$fqinfo\n";
  	$TloxpR2 ++;
 }
}
close IN2;
close OUT1;

$Tnum = $Tnum1+$Tnum2;
$TloxpR = $TloxpR1+$TloxpR2;
$TloxpN = $Tnum - $TloxpR;
print "Total number of reads:$Tnum\nTotal number of reads with loxp:$TloxpR\n";
print "Total number of reads without loxp:$TloxpN\n";
print OUT2 "Total number of reads:$Tnum\nTotal number of reads with loxp:$TloxpR\n";
print OUT2 "Total number of reads without loxp:$TloxpN\n";

$/ = "\n";
my $spnum = 0;
my $sp2num = 0;

my $maploxp = 0;
my $mapnum = 0;
my $map9R = 0;
my $mapelse = 0;
my $umloxp = 0;

my ($synstar,$synend)=(0,0);
if ($syntype eq "R"){
	 $synstar =  $reflen - $wtlen;
}elsif($syntype eq "L"){
	 $synend = $wtlen;
	}
	  
print OUT3 "####From soapPE####\n";
while (<IN3>){
	chomp;
	my $read1 = $_; 
  $spnum ++;	
	my ($seq,$chr,$site) = (split /\s+/,$read1)[1,7,8];
  my $n = 0;
  
  if($syntype eq "R"){
   $map9R ++ if($chr eq $chrsyn || ($chr eq "$chrwt" && $site > $synstar));##chrIXL = 355300
 }elsif($syntype eq "L"){
 	 $map9R ++ if($chr eq $chrsyn || ($chr eq "$chrwt" && $site < $synend));
 }elsif($syntype eq "W"){
  	$map9R ++ if($chr eq $chrsyn || ($chr eq "$chrwt"));
 	}
 	
  if ($seq =~ /$Deseq{"loxp"}/gi){
  	print OUT3 "$read1\n";
  	$n++;
  	$maploxp ++;
  	}
  my $read2 = <IN3>;
  chomp $read2;
  $spnum ++;
  print OUT3 "$read2\n" if $n==1;
  
  ($seq,$chr,$site) = (split /\s+/,$read2)[1,7,8];
  
  if($syntype eq "R"){
   $map9R ++ if($chr eq $chrsyn || ($chr eq "$chrwt" && $site > $synstar));##chrIXL = 355300
 }elsif($syntype eq "L"){
 	 $map9R ++ if($chr eq $chrsyn || ($chr eq "$chrwt" && $site < $synend));
 }elsif($syntype eq "W"){
  	$map9R ++ if($chr eq $chrsyn || ($chr eq "$chrwt"));
 	}
  
  if ($seq =~ /$Deseq{"loxp"}/gi){
  	if ($n==0){
  		print OUT3 "$read1\n";
  	  print OUT3 "$read2\n";
  	  } 		
  	$maploxp ++;
  	}
  }
close IN3;
#print "PE loxp: $maploxp\n";
 
print OUT3 "####From soapSE####\n"; 
my $temp = <IN4>;
chomp $temp;
$sp2num ++;
my ($tempid,$seq,$chr,$site) = (split /\s+/,$temp)[0,1,7,8];

if($syntype eq "R"){
   $map9R ++ if($chr eq $chrsyn || ($chr eq "$chrwt" && $site > $synstar));##chrIXL = 355300
 }elsif($syntype eq "L"){
 	 $map9R ++ if($chr eq $chrsyn || ($chr eq "$chrwt" && $site < $synend));
 }elsif($syntype eq "W"){
  	$map9R ++ if($chr eq $chrsyn || ($chr eq "$chrwt"));
 	}
  
 if ($seq =~ /$Deseq{"loxp"}/gi){
    print OUT3 "$temp\n"; 		
  	$maploxp ++;
  	}

while (<IN4>){
  chomp;
  $sp2num ++;
  my($readid,$seq,$chr,$site) = (split /\s+/,$_)[0,1,7,8];
  $readid = (split /\//,$readid)[0];

  if($syntype eq "R"){
   $map9R ++ if($chr eq $chrsyn || ($chr eq "$chrwt" && $site > $synstar));##chrIXL = 355300
 }elsif($syntype eq "L"){
 	 $map9R ++ if($chr eq $chrsyn || ($chr eq "$chrwt" && $site < $synend));
 }elsif($syntype eq "W"){
  	$map9R ++ if($chr eq $chrsyn || ($chr eq "$chrwt"));
 	}
  
  if ($seq =~ /$Deseq{"loxp"}/gi){
   if($tempid =~ /$readid/){
     print OUT3 "$temp\n"; 		
     print OUT3 "$_\n";
     $temp = $_;
     $tempid = $readid;
   }else{
   	 my $line = <IN4>;
   	 $sp2num ++;
   	 chomp $line;
   	 my($readids,$seqs,$chrs,$sites) = (split /\s+/,$line)[0,1,7,8];
   	 $readids = (split /\//,$readids)[0];
    if($syntype eq "R"){
       $map9R ++ if($chr eq $chrsyn || ($chr eq "$chrwt" && $site > $synstar));##chrIXL = 355300
    }elsif($syntype eq "L"){
 	     $map9R ++ if($chr eq $chrsyn || ($chr eq "$chrwt" && $site < $synend));
    }elsif($syntype eq "W"){
  	   $map9R ++ if($chr eq $chrsyn || ($chr eq "$chrwt"));
 	  }
 	  
     $maploxp ++ if($seqs =~ /$Deseq{"loxp"}/gi);
     if($readid =~ /$readids/){
     	print OUT3 "$_\n"; 
      print OUT3 "$line\n";
     }else{
     	 if ($seqs =~ /$Deseq{"loxp"}/gi){print OUT3 "$line\n";}   	  
   	   print OUT3 "$_\n";
   	 }
   	 $temp = $line;
     $tempid = $readids; #FCC0CH9ACXX:7:2308:1873:179645#TAATTACC/2
   }
  	 $maploxp ++; 	 
  }else{
  	 $temp = $_;
     $tempid = $readid;}   
 }
$mapnum =  $spnum + $sp2num;
$mapelse = $mapnum - $map9R;
$umloxp = $TloxpR - $maploxp;
print OUT2 "Total number of reads mapped:$mapnum\n";
print OUT2 "Total number of reads mapping $chrwt $syntype:$map9R\n";
print OUT2 "Total number of reads mapping elsewhere:$mapelse\n";
print OUT2 "Total number of reads with loxp mapped:$maploxp\n";
print OUT2 "Total number of reads with recombination:$umloxp\n";
print "Total number of reads mapped:$mapnum\n";
print "Total number of reads mapping $chrwt $syntype:$map9R\n";
print "Total number of reads mapping elsewhere:$mapelse\n";
print "Total number of reads with loxp mapped:$maploxp\n";
print "Total number of reads with recombination:$umloxp\n";

close IN4;
close OUT2;
close OUT3;
