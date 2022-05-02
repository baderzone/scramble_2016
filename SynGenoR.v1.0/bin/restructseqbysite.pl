#! /usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use Cwd;

use Bio::SeqIO;

my $usage=<<USAGE;
Describe:Only used to refactor a sequence with structure varations.
  
  Author: Wang Yun, wangyun\@genomics.org.cn
 Version: 1.0  Date: 2012-11-05
 
   Usage:perl $0 <ref_fa> <region> <encod> [option]
         
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

my ($ref_fa,$region,$encod) = @ARGV;
 
open IN2, $region ||die $!;
open IN3, $encod || die $!;

$outdir ||= getcwd;
open OUT,">$outdir/$sample.reseq" || die $!;

#print "loading data...\n";

my %regseq;
my %refseq;
my $n=0;

my $seqio =  Bio::SeqIO->new(-file=>"$ref_fa",-format => 'Fasta')
     or die "couldn't create  Bio::SeqIO $!";
          

while (my $seqobj = $seqio->next_seq()){
   my $seq_id = $seqobj->id;
   my $seq_seq = $seqobj->seq;
   $seq_seq = uc($seq_seq);
   $refseq{$seq_id} = $seq_seq;
   my $seqlen = length ($seq_seq);
   
   print "$seq_id length: $seqlen\n";
}

open ANNOSEGM,">$outdir/annosegm.fasta" || die $!;

while(<IN2>){
  chomp;
  my ($seqid,$star,$END,$indx)=(split /\s+/)[1,2,3,4];
  my $reglen = $END - $star +1;
   
  $regseq{$indx} = substr($refseq{$seqid},$star-1,$reglen);
  print ANNOSEGM ">$indx\_$seqid\n$regseq{$indx}\n";
}

my $restructseq = <IN3>;
chomp $restructseq;
my $encode = (split /\s+/,$restructseq)[1];
my @codeord = split /,/,$encode;

my $reseq;	
foreach (@codeord){
  if ($_ == 0){
     next;
  }elsif($_>0){
     $reseq .= $regseq{$_};
  }elsif($_<0){
     my $idx = abs($_);
     my $revcom = reverse($regseq{$idx});
     $revcom =~ tr/ATGCatgc/TACGtacg/;
        $reseq .= $revcom;
  }
 }

if($reseqid){ 
  print OUT ">$reseqid\n";
}else{
  $reseqid = "synII"."_redup";
  print OUT "$reseqid\n";
}

my $lenscb = length($reseq);
print "Refactoring sequence length:$lenscb\n";
print OUT "$reseq\n";

close IN2;
close IN3;
close OUT;

#===end===	
