#! /usr/bin/perl -w
use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use List::Util qw/max min/;
use Cwd;

my $usage=<<USAGE;

Describe: sorting gvr or spr file...

  Author: Yun Wang, wangyun\@genomic.org.cn
 Version: 1.0, Date: 2012-11-05

Usage:perl $0  <infile>  [-option]
          
      -spid     <str>  : Sample ID/prefix of output file  
      -type     <str>  : run sorting pattern
      -outfile  <str>  : output file
      -outdir          : Directory of output [default: ./]   
      -help            : show this message

USAGE

my ($sample,$type,$outfile,$outdir,$help);
GetOptions( 
   "spid:s"   => \$sample, 
   "type:s"   => \$type,
   "outfile:s"=> \$outfile,
   "outdir:s" => \$outdir,
   "help:s"   => \$help,
);

die $usage if (@ARGV != 1 ||$help);
 
 my $sortfile = $ARGV[0];
 $outdir ||= getcwd;
 $outfile ||= "$outdir/$sample.sort.$type";
 
 open IN, $sortfile ||die $!;
 open OUT, ">$outfile"||die $!;
 
 if($type =~ /gvr/i){
 	  sortgvr();
 }elsif($type =~ /spr/i){
    sortspr();
}

 close IN;
 close OUT;

#######sub function########
 sub sortgvr
 ##===========##
 { 
   my %gvrsite;
   while(<IN>){
 	   chomp;
 	   my ($lesite,$risite) = (split/\s+/,$_)[4,7];
     my $minsite = min($lesite,$risite);
     $gvrsite{$minsite} .= "$_\n";
   }
 
   my @all_key = keys (%gvrsite);
   my @keysall = sort { $a <=> $b } @all_key;
   foreach my $key(@keysall){
     print OUT "$gvrsite{$key}"; 
  }
  system "rm  $outdir/$sample.gvr.temp";
}

 ##===========##
 sub sortspr
 ##===========##
 { 
 	 my %sprsite;
   while(<IN>){
   	 chomp;
   	 my ($lesite,$risite) = (split/\s+/,$_)[2,6];
     my $minsite = min($lesite,$risite);
     my $maxsite = max($lesite,$risite);
     $sprsite{$minsite} .= "$_\n" ;
   }
   my @all_key = keys (%sprsite);
   my @keysall = sort { $a <=> $b } @all_key;

   foreach my $key(@keysall){
        print OUT "$sprsite{$key}"; 
    }
}
#==END==#
