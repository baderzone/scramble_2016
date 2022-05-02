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

Describe:Scaning PCRtag... 

  Author: Yun Wang, wangyun\@genomic.org.cn
 Version: 1.0, Date: 2012-11-11

Usage:perl $0 [-option]
  
      -synreffa  <str>  : Synthetic reference <FA format>
  
      -taggff    <str>  : PCRtag annotation file <gff format>
       
      -aligns    <str>  : alignment file  <BAM or SAM format>

      -outfile   <str>  : statistic output of file   
      
      -wttag     <str>  : wild-type PCRtag file
      
      -logfile   <str>  : logs of for data trace <SAM format>
            
      -chrid     <str>  : synthetic chromsome id,which indicated the id is used to the alignment 
      
      -taglen    [int]  : PCRtag length [default 28]
       
      -mrds      <str>  : minimal number of reads allowed for surpporting wild-type PCRtag
      
      -pvalue    <float>: threshold of cumulative probabilities for a wild-type PCRtag #poisson model for count
      
      -help             : show this message

Example:
  perl $0 -synreffa yeast_chr02_3_22.FASTA -taggff yeast_chr02_3_22.TAGS.gff -aligns synIIAY.clean.sort.bam -outfile synIIAY.tags.align.stat -wttag synIIAY.wttag.stat -logfile synIIAY.loxp.align.sam -chrid yeast_chr02_3_22 -mrds 5  -pvalue 0.001

USAGE

my ($synreffa,$taggff,$aligns,$outfile,$logfile,$wttag,$chrid,$taglen,$mrds,$pvalue,$help);
GetOptions(
   "synreffa:s" =>  \$synreffa,
   "taggff:s"   =>  \$taggff,
   "aligns:s"   =>  \$aligns,
   "outfile:s"  =>  \$outfile,
   "logfile:s"  =>  \$logfile,
   "wttag:s"    =>  \$wttag,
   "chrid:s"    =>  \$chrid,
   "taglen:s"   =>  \$taglen,
   "mrds:s"     =>  \$mrds,
   "pvalue:s"   =>  \$pvalue,
   "help"       =>  \$help,
   );
   
die $usage if (!$synreffa ||!$taggff ||!$aligns ||!$chrid ||$help);

$outfile ||= "tags.align.stat";
$logfile ||= "tags.align.sam";
$wttag ||= "tags.wt";
$taglen ||= 28;
$mrds ||= 5;
$pvalue ||= 0.001;

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

open TAGS,"$taggff" || die $!;
open OUT,">$outfile" || die $!;
open LOGS,">$logfile" || die $!;
open WTTAG,">$wttag" || die $!;

##
my %tagseq;
my %tagbase_syn;
my %tagbase_wt;
my %tagdep;
my %tagcate;
my %tagsite;
my %tagrds;
my %temps;

while (<TAGS>) {
  chomp;
  my ($main_feat,$star,$ends,$anno) = (split /\s+/,$_)[2,3,4,8];  
  next if $main_feat eq "PCR_product";
  my @items = split /\;/,$anno;
  my $item = pop (@items);  
  my ($feat,$value) = (split /=/,$item);    
 
  if ($feat eq "wtseq"){
    $tagseq{$star}{"wt"} = $value; 
  }else{
    print "warning: reading wild-type tag seq error,$star..$ends\n";
  }
  $tagseq{$star}{"syn"} = substr ($refseq,$star-1,$ends-$star+1);
  
  my $i = $star;
  foreach (split //,$tagseq{$star}{"wt"}){
    $tagbase_wt{$i} = $_;
    $i ++;
  }
  
  $i = $star;
  foreach (split //,$tagseq{$star}{"syn"}){
    $tagbase_syn{$i} = $_;
    $i ++;
  }
    
  ## hash tagsite initialization

  foreach ($star..$ends){
    $tagdep{$_} = 0; #loxp belongs and depth                      
    push @{$tagcate{$_}},$star unless exists $temps{$_}{$star}; ##~~ doesn't exists for perl 5.8.8   
    #note: some tags have overlap...
    $temps{$_}{$star} = 1;
  }
    
  $tagrds {$star}{"syn"} = 0;
  $tagrds {$star}{"wt"} = 0;
  $tagrds {$star}{"other"} = 0;
}
undef %temps;

while(<IN>){
	#3I5RT:00855:02055       16      chr01   4044    42      201M1D169M      *       0       0       ACATCAAAT
  chomp;
  my $rdsitem = $_;
  my @rdsitems = split /\s+/,$rdsitem;
  
  if($rdsitems[2] eq $chrid){
    my $rdslen = length ($rdsitems[9]);
    my $n = int ($rdslen/$taglen);
    my $chc = 0;
    
    my @tagstar;
    my %temp;
    foreach my $i(0..$n){
      my $survey_site = $rdsitems[3] + $i * $taglen;
      if(exists $tagdep{$survey_site}){
        $chc = 1;
        foreach (@{$tagcate{$survey_site}}){
          if (exists $temp{$_}){
             next;
          }else{
             push @tagstar,$_;
             $temp{$_} = 1;
          }
        }
      }
    }
           
    if ($chc ==1){
      my @digital = split /[MID]/,$rdsitems[5];
      my @MDtype  = split /\d+/,$rdsitems[5];
      shift @MDtype;
      my %tag_syn_match;
      my %tag_wt_match;
      foreach (@tagstar){
         $tag_syn_match{$_} = 0;
         $tag_wt_match{$_} = 0;
      }
      
      my $tag_match = "";
      my $site = $rdsitems[3];
      my $rds_site = 1;
      
      foreach my $i(0..@digital-1){
          if ($MDtype[$i] eq "M"){
            foreach (0..$digital[$i] -1){
               if(exists $tagdep{$site}){
                   $tagdep{$site} ++;
                   $tag_match .= "1";
                   
                   my $seqbase = substr ($rdsitems[9],$rds_site - 1,1);
                   foreach(@{$tagcate{$site}}){
                     
                     $tag_syn_match{$_} ++ if $tagbase_syn{$site} eq $seqbase;
                     $tag_wt_match{$_} ++ if $tagbase_wt{$site} eq $seqbase;
                   }
               }else{
                  $tag_match .= "#";
               }
               $site ++;
               $rds_site ++;
            }
          }elsif($MDtype[$i] eq "D"){
             foreach (0..$digital[$i] -1){
               if(exists $tagdep{$site}){
                   $tag_match .= "0";             
               }else{
                   $tag_match .= "-";
               }
               $site ++;                         
             }         
         }elsif($MDtype[$i] eq "I"){
             foreach (0..$digital[$i] -1){             
               $rds_site ++;                         
             }    
        }   
      }
      
      $tag_match = "X".$tag_match."X";
      my ($upstr,$dnstr) =  (split /\d+/,$tag_match)[0,-1];

       if (defined $upstr){
          shift @tagstar if length ($upstr) -1 < 10;
       }else{ shift @tagstar;}
             
       
       if (defined $dnstr){
          pop @tagstar if length($dnstr) -1 < 10; 
       }else{pop @tagstar;}
      
       if(@tagstar){  
         foreach(@tagstar){
         	    
           if ( $tag_syn_match{$_} >= $tag_wt_match{$_} && $tag_syn_match{$_} >= 0.75*$taglen){
             $tagrds{$_}{"syn"} ++;
           }elsif( $tag_syn_match{$_} < $tag_wt_match{$_} && $tag_wt_match{$_} >= 0.60*$taglen){
             $tagrds{$_}{"wt"} ++;
           }else{
             $tagrds{$_}{"other"} ++;             	
           }
         }
       
         print LOGS "$rdsitem\n";
       }
        
    }
  }else{
    next;
  }
}

foreach my $key(sort {$a<=>$b} keys %tagrds){
	 my $keyend = $key + length ($tagseq{$key}{"syn"})- 1;
	 my $mean = ($tagrds{$key}{"syn"}+$tagrds{$key}{"wt"})/2;
	 my $p_value = 1 - ppois($tagrds{$key}{"wt"},$mean);  #poisson model for count
	 print OUT "$chrid\ttag\t$key\t$keyend\t$tagrds{$key}{syn}\t$tagrds{$key}{wt}\t$tagrds{$key}{other}\t$p_value\n";

	 if($tagrds{$key}{"wt"} >= $mrds || $p_value < $pvalue){
	     print WTTAG "$chrid\ttagScan\tPCRtag.wt\t$key\t$tagseq{$key}{syn}\t$tagseq{$key}{wt}\t$tagrds{$key}{syn}\t$tagrds{$key}{wt}\t$tagrds{$key}{other}\t$p_value\n";	     	     
	 }
}

close IN;
close OUT;
close LOGS;
close WTTAG;




