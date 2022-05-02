#! /usr/bin/perl -w
use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use Cwd;

my ($len,$help,$outdir);
GetOptions(
   "l:s" => \$len,
   "o:s" => \$outdir,
   "h:s" => \$help,
);
if (@ARGV != 1 ||$help){
	die"Describe:Only used for transforming loxplist to loxp region.
Usage:perl $0 <loxplist> -l [chr_len] [-o (output path)]";
}
print "loading...\n";
my $loxplist = $ARGV[0];
$outdir ||= getcwd;
open IN,$loxplist or die $!;

my @suffixlist = qw(.lxpftr.lst .txt .list);
my $basename = basename ($loxplist, @suffixlist);
open OUT,">$outdir/$basename.loxpreg";

my $regstar = 1; 
my $n = 1;
my $chrid;
while(<IN>){
	chomp;
	my ($star,$END);
	($chrid,$star,$END)= (split /\s+/)[2,3,4];
	my $regend = $star+16;
	print OUT "$chrid\tloxp\tloxpreg\t$regstar\t$regend\t.\t.\t.\t$n\n";
	$regstar = $regend+1;
	$n ++;
}

print OUT "$chrid\tloxp\tloxpreg\t$regstar\t$len\t.\t.\t.\t$n\n";
close IN;
close OUT;

	
