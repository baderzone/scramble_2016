#! /usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use lib "$Bin/../lib"; #for math::CDF model
use Math::CDF qw(:all);
use Bio::SeqIO;
use Cwd;

my $usage=<<USAGE;

Describe:Scan loxPsym... 

  Author: Yun Wang, wangyun\@genomic.org.cn
 Version: 1.0, Date: 2012-11-11

Usage:perl $0 [-option]
  
      -synreffa  <str>  : Synthetic reference <FA format >
  
      -aligns    <str>  : alignment file  <BAM/SAM or SOAP format >

      -outfile   <str>  : statistic output of file      
      
      -loxpdel   <str>  : loxp deletion file
      
      -logfile   <str>  : logs of for data trace 
      
      -DNAe      <str>  : A inducible DNA elements,such as loxp and rox; #Only loxp is available.
      
      -mrds      <str>  : minimal number of reads allowed for surpporting loxp.del 
      
      -pvalue    <float>: threshold of cumulative probabilities for a loxp.del #poisson model for count
      
      -help             : show this message

Example:
  perl $0 -synreffa yeast_chr02_3_22.FASTA -aligns synIIAY.clean.sort.bam -outfile synIIAY.loxp.align.stat -loxpdel synIIAY.loxp.align.stat -logfile synIIAY.loxp.align.sam -DNAe loxp -mrds 5  -pvalue 0.001

USAGE

my ($synreffa,$aligns,$outfile,$logfile,$loxpdel,$chrid,$DNAe,$mrds,$pvalue,$help);
GetOptions(
   "synreffa:s" =>  \$synreffa,
   "aligns:s"   =>  \$aligns,
   "outfile:s"  =>  \$outfile,
   "logfile:s"  =>  \$logfile,   
   "loxpdel:s"  =>  \$loxpdel,
   "chrid:s"    =>  \$chrid,
   "DNAe:s"     =>  \$DNAe,  
   "mrds:s"     =>  \$mrds,
   "pvalue:s"   =>  \$pvalue,
   "help"       =>  \$help,
   );
   
die $usage if (!$synreffa ||!$aligns||!$chrid||$help);

$DNAe ||= "loxp";
$outfile ||="loxp.align.stat";
$logfile ||="loxp.align.sam";
$loxpdel ||="loxp.del";
$mrds ||=5;
$pvalue ||=0.001;

##
my $inref  = Bio::SeqIO->new(-file   => "$synreffa",
                             -format => 'Fasta');
 
my $refseq;
while (my $inref_item = $inref->next_seq()){
   if ($inref_item -> id eq $chrid){
     $refseq = uc($inref_item -> seq);   
   }
}

die "Can't find the $chrid sequence in $synreffa \n" unless defined $refseq;	
my $reflen = length($refseq);

if($aligns =~ /.*\.sam/){
  open IN, $aligns ||die $!;
}elsif ($aligns =~ /.*\.bam/){
	use Env qw(SAMTOOLS);
	
   (-s $SAMTOOLS && -e $SAMTOOLS) ||  die "error:No executable samtools program for read BAM file...
please write a executable samtools to environment variable \$SAMTOOLS like this:
export SAMTOOLS=\"path_to_samtools\"
or input SAM format file
program quit.. 
";  
  open IN, "$SAMTOOLS view $aligns |" ||die $!;
}else{
  die print "unrecognized file can't to be read...\n$usage";	
}

open OUT, ">$outfile" ||die $!;
open LOGS, ">$logfile" ||die $!;
open LOST,">$loxpdel" ||die $!;

##

my %Deseq;
$Deseq{"loxp"} = "ATAACTTCGTATAATGTACATTATACGAAGTTAT";
 
my %loxpdep;
my %loxpcate;
my %loxpsite;
my %loxprds;

my $offset = 0;
my $result = 0;
while ($result != -1) {
   
  $result = index($refseq, $Deseq{$DNAe}, $offset);
  my $star =  $result + 1;
  my $ends  =  $star + 33;
  #print OUT "FTR\tloxp\t$chrid\t$star\t$ends\n" if $star > 0;

  ## hash loxpsite initialization
  if($star >0){
    foreach ($star..$ends){
       $loxpdep{$_} = 0; #loxp belongs and depth                      
       $loxpcate{$_} = $star;
    }
    push @{$loxpsite{$star}},($star..$ends);
    $loxprds {$star}{"loxp"} = 0;
    $loxprds {$star}{"noloxp"} = 0;
    $loxprds {$star}{"other"} = 0;
  }  
  $offset = $result + 1;
}

while(<IN>){
	#3I5RT:00855:02055       16      chr01   4044    42      201M1D169M      *       0       0       ACATCAAAT
  chomp;
  my $rdsitem = $_;
  my @rdsitems = split /\s+/,$rdsitem;
  
  if($rdsitems[2] eq $chrid){
    my $rdslen = length ($rdsitems[9]);
    my $n = int ($rdslen/30);
    my $chc = 0;
    
    my @loxpstar;
    my %temp;
    foreach my $i(0..$n){
      my $survey_site = $rdsitems[3] + $i *30;
      if(exists $loxpdep{$survey_site}){
        $chc = 1;
        push @loxpstar,$loxpcate{$survey_site} unless (exists $temp{$loxpcate{$survey_site}});             
      }  
    }
    my %match_base;
    my %match_base_no;
    foreach (@loxpstar){
       $match_base{$_} = 0;
       $match_base_no{$_} = 0;
    }
         
    if ($chc ==1){    	
      my @digital = split /[MID]/,$rdsitems[5];
      my @MDtype  = split /\d+/,$rdsitems[5];
      shift @MDtype;
      my $loxp_match = "";
      my $site = $rdsitems[3];
      
      foreach my $i(0..@digital-1){
          if ($MDtype[$i] eq "M"){
            foreach (0..$digital[$i] -1){
               if(exists $loxpdep{$site}){
                   $loxpdep{$site} ++;
                   $loxp_match .= "1"; 
                   $match_base{$loxpcate{$site}} ++;               
               }else{
                 $loxp_match .= "#";
               }
               $site ++;
            }
          }elsif($MDtype[$i] eq "D"){
             foreach (0..$digital[$i] -1){
               if(exists $loxpdep{$site}){
                   $loxp_match .= "0";
                   $match_base_no{$loxpcate{$site}} ++                  
               }else{
                 $loxp_match .= "-";
               }
               $site ++;          
             }         
         }   
      }
      $loxp_match = "X".$loxp_match."X";
      my ($upstr,$dnstr) =  (split /\d+/,$loxp_match)[0,-1];
      if (defined $upstr){
         shift @loxpstar if length ($upstr) -1 <10;
      }else{ shift @loxpstar;}
             
       
      if (defined $dnstr){
         pop @loxpstar if length($dnstr) -1 < 10; 
      }else{pop @loxpstar;} 
      	
     if(@loxpstar){
        foreach (@loxpstar){      
         if($match_base{$_} + $match_base_no{$_} >30){
           $match_base{$_} >= 24 || $match_base_no{$_} <= 5? $loxprds {$_}{"loxp"} ++:
           $match_base_no{$_} >= 15 ? $loxprds {$_}{"noloxp"} ++:
                                   $loxprds {$_}{"other"} ++
         }
        }
         print LOGS "$rdsitem\n";
      }
           
    }
  }else{
    next;
  }
}

foreach my $key(sort {$a<=>$b} keys %loxprds){
	 my $keyend = $key + 33;
	 my $mean = ($loxprds{$key}{loxp}+$loxprds{$key}{noloxp})/2;
	 my $p_value = 1 - ppois($loxprds{$key}{noloxp},$mean);  #poisson model for count
	 print OUT "FTR\tloxp\t$chrid\t$key\t$keyend\t$loxprds{$key}{loxp}\t$loxprds{$key}{noloxp}\t$loxprds{$key}{other}\t$p_value\n";

	 if($loxprds{$key}{noloxp} >= $mrds || $p_value < $pvalue){
	     my $refsite = $key -1;
	     my $refbase = substr ($refseq,$refsite-1,35);
	     my $seqbase = substr ($refseq,$refsite-1,1);
	     print LOST "$chrid\tloxpScan\tloxp.del\t$refsite\t$refbase\t$seqbase\t$loxprds{$key}{loxp}\t$loxprds{$key}{noloxp}\t$loxprds{$key}{other}\t$p_value\n";	     	     
	 }
}

close IN;
close OUT;
close LOGS;





