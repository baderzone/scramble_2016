package CNstat;
# This module is used to test copy number stat.
use strict;
use List::Util qw/max min/;


use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@EXPORT_OK = qw(cnstat); 

$VERSION     = 1.00;
@ISA         = qw(Exporter);


sub cnstat{
	my ($chrid,$cutlen,$edgeinfo,$coverage,$refcover,$loxpregion) = @_;
	
  open SCVG, "$coverage" ||die $!;
  open CVG, "$refcover" ||die $!;
  open LRG, "$loxpregion" ||die $!;
  
  my @cns;

  ##loadig loxp region###
  my %loxp_region;
  my $reflen=0;
  $/ = "\n";
  
  while (<LRG>){
     chomp;
	   my ($chrid,$star,$END,$regid) = (split /\s+/,$_)[0,3,4,8];
     push @{$loxp_region{$regid}},($chrid,$star,$END);   
    
     $reflen = $END  if $reflen <$END;
  }
  my $lastindex = max (keys %loxp_region); 
  close LRG;

  ##loading refcoverage##
  $/ =">";
  <CVG>;
  my %chrdepth;
  my %chrdepth_testee;

  while(<CVG>){
     chomp;
     next if ($_ =~ /^#/);
     my $depth = ();
     my @line = split /\n/;
     my $chrID = shift @line;
     next if ($chrID ne $chrid);
     foreach (@line){
       $depth .="$_";		
     }       

     @{$chrdepth{$chrid}} = split /\s+/,$depth;
     my $deplength =  @{$chrdepth{$chrid}};      
     last;
  }   
  $/ ="\n";
  close CVG;

  $/ =">";
  <SCVG>;
  while(<SCVG>){
     chomp;
     next if ($_ =~ /^#/);
     my $depth = ();
     my @line = split /\n/;
     my $chrID = shift @line;
     next if ($chrID ne $chrid);
     foreach (@line){
       $depth .="$_";		
     }       

     @{$chrdepth_testee{$chrid}} = split /\s+/,$depth;
     my $deplength =  @{$chrdepth{$chrid}};      
     last;
  }   
  $/ ="\n";
  close SCVG;
 
  my %testee_depth;
  my %copynum;
  my %ref_depth;
  my %regidxedg;

  foreach (keys %{$edgeinfo}){
	  chomp;
	  my @edge_info = split /\s+/,${$edgeinfo}{$_};
	  my @regidx = (split /,/,$edge_info[3]);
  
    my @site;
    foreach (@regidx){
  	
      $regidxedg{$_} = $edge_info[0];
      my $starsite = ${$loxp_region{$_}}[1];
      my $endsite = ${$loxp_region{$_}}[2];  
      
      my $average_depth_testee = averdep($chrid,$starsite,$endsite,\%chrdepth_testee);
      my $average_depth_ref = averdep($chrid,$starsite,$endsite,\%chrdepth);
      $testee_depth{$_} = $average_depth_testee;
      $ref_depth{$_} = $average_depth_ref;
      $copynum{$_} = $edge_info[1];
  
      my $idx = $_;    
      if($idx!=1&&$idx!= $lastindex &&($endsite-$starsite+1)>=2*$cutlen){
        my $n = int(($endsite-$starsite +1)/$cutlen);
        my $starsitenew = $starsite;
        my $endsitenew = $starsite+$cutlen - 1;
        my $key;
        foreach(1..$n-1){    	
          $key = "$idx".".$_";
          $regidxedg{$key} = $edge_info[0];
          $average_depth_testee = averdep($chrid,$starsitenew,$endsitenew,\%chrdepth_testee);
          $average_depth_ref = averdep($chrid,$starsitenew,$endsitenew,\%chrdepth);
          $testee_depth{$key} = $average_depth_testee;
          $ref_depth{$key} = $average_depth_ref;
          $copynum{$key} = $edge_info[1];
          $starsitenew  =  $endsitenew + 1;
          $endsitenew  = $starsitenew + $cutlen - 1;
        }
    
       $key = "$idx".".$n";
       $average_depth_testee = averdep($chrid,$starsitenew,$endsite,\%chrdepth_testee);
       $average_depth_ref = averdep($chrid,$starsitenew,$endsite,\%chrdepth);
       $testee_depth{$key} = $average_depth_testee;
       $ref_depth{$key} = $average_depth_ref;
       $copynum{$key} = $edge_info[1];        
       $regidxedg{$key} = $edge_info[0];
     }
   }
  }

  ##### reckon copy number(CN)#####
  my ($CNpot,$regkey) = reckonCN(\%testee_depth,\%ref_depth,\%copynum);

  foreach (@{$regkey}){ 
    push @cns, "$regidxedg{$_}\t${$CNpot}{$_}\n";
  }
  
 return @cns;
 
##=====================================########
##==========sub function===============########
##=====================================########
####caculate the average sequencing depth#####
  sub averdep{
    my ($chrID,$star,$END,$chrdepth)=@_;
    my $averdepth = 0;
    if(exists ${$chrdepth}{$chrID}){
  	  my $seqlen = $END -$star + 1;
      my $subtotaldep = 0;
      foreach($star-1..$END-1){
    		$subtotaldep += ${${$chrdepth}{$chrID}}[$_];
    	}
      $averdepth = int(($subtotaldep/$seqlen)*100)/100;    
    
    }else{
	    print "Chromosome ID error!\n";
    }
    return $averdepth;
  }

######sort keys########
  sub sortkey{
	  my ($hash,$opern) = @_;
    my @keyid = keys %{$hash};
    my %newkey;
    my @newkeyid;
    if ($opern ==1){
      foreach (@keyid){
      	my $findex = (split /_/,$_)[0];      	 
        push @{$newkey{$findex}},$_;
      }
    }elsif($opern==0){ 
  
      foreach (@keyid){
      	my $findex = (split /\//,$_)[0];      	 
        push @{$newkey{$findex}},$_;
      } 
    }
  
   foreach my $key(sort{$a<=>$b} keys %newkey){
     foreach (@{$newkey{$key}}){
       push  @newkeyid,$_;
     }
   }
  return @newkeyid;
  }

######caculate copy number####
  sub caculateCN{
    my ($testee_depths,$ref_depth,$testee_CN,$ref_CN,$regionkeys,$opern) = @_;
    my (%region_CN,%CN_CV);

    my %CNpot;  
    foreach my $i(@{$regionkeys}){
   	  my %copy_number;
   	  my $tempCNpot;
   	  foreach my $j(@{$regionkeys}){
   	  	  if ($i eq $j){
   	  	  $copy_number{$j} = "--"; next;}
          unless (${$testee_depths}{$j} == 0|| ${$ref_depth}{$i} ==0 || ${$ref_CN}{$j} == 0) {
              my $copy_num = (${$ref_depth}{$i} / ${$testee_depths}{$j}) * (${$testee_depths}{$i} / ${$ref_depth}{$i}) * (${$ref_CN}{$i} / ${$ref_CN}{$j}) * ${$testee_CN}{$j};
              $copy_num = sprintf("%.3f", $copy_num);
              $copy_number{$j} = $copy_num;
              $tempCNpot .= "\t$copy_num";

          }else{
              $tempCNpot .= "\t1";
          }
      }

      $CNpot{$i} =  $tempCNpot;
    }
   return (\%CNpot);       	
 }


#####reckon copy number####

  sub reckonCN{
	  my ($testee_depth,$ref_depth,$copynum) = @_;
    my %edge_CN;
   
    my @regionkeys = sortkey(\%{$testee_depth},1);	
 
    foreach (@regionkeys){
    	if (${$copynum}{$_}==0){
   		  delete ${$testee_depth}{$_};
   		  delete ${$ref_depth}{$_};
   		}
    }
    
    @regionkeys=();
    @regionkeys = sortkey(\%{$testee_depth},1);	
  
    my (%ref_CN,%testee_CN);
    foreach (@regionkeys){
   	  $ref_CN{$_} = 1;
   	  $testee_CN{$_} = ${$copynum}{$_};
   	}
  
    my $CNpot = caculateCN(\%{$testee_depth},\%{$ref_depth},\%testee_CN,\%ref_CN,\@regionkeys,1);  
    return (\%{$CNpot},\@regionkeys);
  }
 }
 