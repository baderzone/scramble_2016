#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use PerlIO::gzip;
use FindBin qw ($Bin);
use lib "$Bin/../lib";
use lib "$Bin/../lib/perllib";
use List::Util qw/max min/;
use Statistics::Distributions;
use Modelcheck qw /Ttestcheck Utestcheck Zscorecheck Mascheck/;
use CNstat qw /cnstat/;
use Cwd;

my $usage=<<USAGE;
Describe: This script is used to stat the edges/nodes and parse the input for PATH .
   
  Author: Yun Wang, wangyun\@genomic.org.cn
 Version: 1.0, Date: 2012-11-05

Usage :perl $0 
        -prefix     :   prefix of output file [default: parsepath] 
        -chrid      :   synthetic chromosome ID       
        -coverage   :   soapcoverage file of sample
        -refcover   :   soapcoverage file of reference sample(unscrambled sample)
        -abnorm     :   abnormal mapped split reads (*.sort.gvr)  
        -soap       :   pair_end soap result including normal mapped loxp reads
        -soap2      :   single_end soap result including normal mapped loxp reads
        -loxpregion :   loxpregion file
        -minread    :   minimal split reads for supporting a breakpoint [default: 2]
        -mdep       :   maximal average sequencing deth of copy number = 0 [default: 10] 
        -mcycle     :   maximal cycle for estimate copy number [default: 10] 
        -cutlen     :   cut each region into subregion for Mahalanobis Distance(MD) analysis [default: 500] #500 model only used
        -outdir     :   output directory [default: ./]
        -help       :   show this message
USAGE

my ($prefix,$chrID,$coverage,$refcover,$abnorm,$soap,$soap2,$loxpregion,$minread,$mdep,$mcycle,$cutlen,$outdir,$help);
GetOptions(
    "prefix:s"     => \$prefix,
    "chrid:s"      => \$chrID,
    "coverage:s"   => \$coverage,
    "refcover:s"   => \$refcover, 
    "abnorm:s"     => \$abnorm,
    "soap:s"       => \$soap,
    "soap2:s"      => \$soap2,
    "loxpregion:s" => \$loxpregion,
    "minread:s"    => \$minread,
    "mdep:s"       => \$mdep,
    "mcycle:s"     => \$mcycle,
    "cutlen:s"     => \$cutlen,
    "outdir:s"     => \$outdir,
    "help"       => \$help,
);

die $usage if (!$chrID||!$coverage||!$refcover||!$abnorm||!$soap||!$soap2||!$loxpregion||$help);

$prefix ||= "parsepath";
$minread ||= 2 unless ($minread==0);
$mdep ||= 10; 
$mcycle ||= 10;
$cutlen ||= 500;
$outdir ||= getcwd;

open CVG, "$coverage" ||die $!;
open ABN, "$abnorm" ||die $!;
open LRG, "$loxpregion" ||die $!;
open EDGE, ">$outdir/$prefix.edg" ||die $!;
open NODE, ">$outdir/$prefix.nod" ||die $!;


my %loxp_region;
my %loxpnod; #loxp node
my $reflen=0;
my $lastindex=0;
my %pointends;

while (<LRG>){
   chomp;
	 my ($chrid,$star,$END,$regid) = (split /\s+/,$_)[0,3,4,8];
   my $leftnod = $regid*2 -2;
   my $rightnod = $regid*2 -1;
   push @{$loxp_region{$regid}},($chrid,$star,$END);   
   push @{$loxpnod{$regid}},($leftnod,$rightnod);
   $reflen = $END  if $reflen <$END;
   $pointends{$star} = $leftnod;
   $pointends{$END}= $rightnod;
 }
$lastindex = max (keys %loxp_region); 
print "reflen:$reflen\tlastindex:$lastindex  \n";


#========================================
#===========Tansform the node============
#========================================
print "\nTansform the node...\n";
print "\nstep0:Collect breakpoints information...\n";
##find the breakpoints linkage##
my ($left,$right,%breakpoint);
while(<ABN>){
    chomp;
    my ($node_side1,$node_side2,$node);
    my @line = split(/\t/,$_);
    my $read_left_site_m = $line[4] + 17;
    my $read_left_site_l = $line[4] - 17;
    my $read_right_site_m = $line[7] + 17;
    my $read_right_site_l = $line[7] - 17;
           
    $node_side1 = getnod($read_left_site_m,$read_left_site_l,\%loxp_region,\%loxpnod);
    $node_side2 = getnod($read_right_site_m,$read_right_site_l,\%loxp_region,\%loxpnod);
    if($node_side1 && $node_side2 ){
      if($node_side1>$node_side2){$node = "$node_side2/$node_side1";}else{$node = "$node_side1/$node_side2";}
      $breakpoint{$node}++;    
    }else{ print "a illegal reads was ignored...\n";}
}



###check and rectify the abnormal node### 
my %points;

foreach (keys %breakpoint){
   if ($breakpoint{$_} < $minread){
     	delete $breakpoint{$_};
      next;}
   my ($nodeL,$nodeR) = split /\//,$_;
   if ($nodeR==$nodeL + 1 && ($nodeL%2)==1){
   	  delete $breakpoint{$_};
      next;
    }
   print "$_\t$breakpoint{$_}\n"; 
   
   #get the new edge####
   my $ck=0;
   foreach my $region_index(keys %loxpnod){
   	  my $newkeys;
   	  if($nodeL == ${$loxpnod{$region_index}}[0]){
   	  	$points{${$loxp_region{$region_index}}[1]}=$region_index;
   	  	$newkeys= ${$loxp_region{$region_index}}[1] - 1;
   	  	$points{$newkeys} = $region_index - 1;
   	  	$ck ++;
   	  }elsif($nodeL == ${$loxpnod{$region_index}}[1]){
   	  	$points{${$loxp_region{$region_index}}[2]} = $region_index;
   	  	$newkeys= ${$loxp_region{$region_index}}[2] + 1;
   	  	$points{$newkeys} = $region_index + 1;
   	  	$ck ++;
      }
      
      if($nodeR == ${$loxpnod{$region_index}}[0]){
   	  	$points{${$loxp_region{$region_index}}[1]}=$region_index;
   	  	$newkeys= ${$loxp_region{$region_index}}[1] - 1;
   	  	$points{$newkeys} = $region_index - 1;
   	  	$ck ++;
   	  }elsif($nodeR == ${$loxpnod{$region_index}}[1]){
   	  	$points{${$loxp_region{$region_index}}[2]}=$region_index;
   	  	$newkeys= ${$loxp_region{$region_index}}[2] + 1;
   	  	$points{$newkeys} = $region_index + 1;
   	  	$ck ++;
      }
   }
  if ($ck <2){print "Warning! you make a mistake! be carefull the bug3...";}
 }

my @newregion = (sort{$a<=>$b} keys %points);
#print join (",",@newregion);
print "\n";

my $breakpoint_num = @newregion;
print "Warning! you make a mistake! be carefull the bug4..."if (($breakpoint_num%2)==1);
$breakpoint_num = $breakpoint_num/2;

print "Detect $breakpoint_num breakpoints...\n";

print "\nstep1:Get node transformation...\n";

my %newnode;
my %newsididx;
my %newside;
my $newnodeid = 1; 
foreach(@newregion){
	$newsididx{$_} = $newnodeid;
	$newside{$newnodeid} = $_;
	$newnodeid ++;	
 }
 
foreach (keys %breakpoint){
   my ($node_left,$node_right) = split /\//,$_;
   my ($new_node_left,$new_node_right)=(0,0);
     foreach my $pointskey(keys %pointends){
   	    if($node_left == $pointends{$pointskey}){
   	      $new_node_left = $newsididx{$pointskey}; 
   	      last;
   	    }
   	}
   	 foreach my $pointskey(keys %pointends){
   	    if($node_right==$pointends{$pointskey}){
   	      $new_node_right = $newsididx{$pointskey}; 
   	      last;
   	    }
   	} 	
  my $newnodekey = "$new_node_left"."/$new_node_right";	
  $newnode{$newnodekey} = $breakpoint{$_};
  print "$newnodekey\t$newnode{$newnodekey}\n";  
 }

print "\nstep2:stat nomal mapping loxp reads on the breakpoint...\n";
####stat nomal mapping loxp reads on the breakpoint#######
my %nomal_node;
foreach my $n(1..$breakpoint_num){
	  my $nomal_node_left =2*$n-1;
	  my $nomal_node_right =2*$n;
	  my $nomal_node_key = "$nomal_node_left"."/$nomal_node_right";
	  $nomal_node{$nomal_node_key} = "$newside{$nomal_node_left}\t$newside{$nomal_node_right}";
	 
	}
my $nomal_node_stat = statnomalnod($chrID,\%nomal_node,$soap,$soap2);	

foreach (keys %{$nomal_node_stat}){
	
	 if(${$nomal_node_stat}{$_} >= $minread){
	 	  $newnode{$_} = ${$nomal_node_stat}{$_};
	 	}
	}

print "\nstep3:pharse the node information...\n";		
 	
my @sortnodekeys = sortkey(\%newnode,0);

foreach (@sortnodekeys){
	print NODE "$_\t$newnode{$_}\n";
	print "$_\t$newnode{$_}\n"; 	
}
	
	
#========================================
#===========Tansform the edge============
#========================================
print "\nTansform the edge...\n";
my $tempsite;
my %edge;
my $findex = regionindex ($newregion[0],\%points);
$edge{"0_1"} = "1\t$newregion[0]\t1\t$findex";
print "0_1\t1\t$newregion[0]\t1\t$findex\n";

foreach my $i(1..$breakpoint_num){
	  my $newedgeL = $i*2;
	  my $newedgeR = $i*2+1;
	  my $siteL = $i*2-1;
	  my $siteR = $i*2;
	  
	  my $newedge = "$newedgeL"."_"."$newedgeR";
	  my $regionindexL = regionindex ($newregion[$siteL],\%points);	  
	 
	  if ($i<$breakpoint_num){
	  	 my $regionindexR = regionindex ($newregion[$siteR],\%points);
	  	 $edge{$newedge} = "$newregion[$siteL]\t$newregion[$siteR]\t$regionindexL\t$regionindexR";     
       print "$newedge\t$newregion[$siteL]\t$newregion[$siteR]\t$regionindexL\t$regionindexR\n";
    }else{
   	   $edge{$newedge} = "$newregion[$siteL]\t$reflen\t$regionindexL\t$lastindex";     
       print "$newedge\t$newregion[$siteL]\t$reflen\t$regionindexL\t$lastindex\n";
    }
}

####caculate average sequencing depth###
my %seqdepth = readcvg($coverage,$chrID);	
my %testee_depth;
my %constant_depth;
foreach my $edgeid(keys %edge){
	my @edginfo = split /\t/,$edge{$edgeid};  
  my $average_depth = averdep($chrID,$edginfo[0],$edginfo[1],\%seqdepth);
  my $indexregion ="$edginfo[2]";
  $edginfo[2] = $edginfo[2]+1;
  foreach ($edginfo[2]..$edginfo[3]){
  	$indexregion .= ",$_";
  }
  $testee_depth{$edgeid} = "$average_depth\t$indexregion";
 }
 %constant_depth = %testee_depth;

%seqdepth = ();
%seqdepth = readcvg($refcover,$chrID);
my %ref_depth;
foreach my $edgeid(keys %edge){
  my @edginfo = split /\t/,$edge{$edgeid};  
  my $average_depth = averdep($chrID,$edginfo[0],$edginfo[1],\%seqdepth);
  $ref_depth{$edgeid} = $average_depth;
}      

##### reckon copy number(CN)#####
my ($edge_CN,$CN_pvalue) = reckonCN(\%testee_depth,\%ref_depth);
my @edge_key = sortkey(\%{$edge_CN},1);

foreach (@edge_key){
	until(exists ${$CN_pvalue}{$_}){${$CN_pvalue}{$_}="NA";}
	print EDGE "$_\t${$edge_CN}{$_}\t${$CN_pvalue}{$_}\t$constant_depth{$_}\n";
}
   
#============================================    
#===============sub function=================
#============================================

####get split reads mapping information#####
sub getnod {
	my ($read_site_m,$read_site_l,$loxp_region,$loxpnod)= @_;
	my $node_side = 0;
	foreach my $region_index(keys %{$loxp_region}){
       if ($read_site_m == ${${$loxp_region}{$region_index}}[2]){
       	  $node_side = ${${$loxpnod}{$region_index}}[1];      	  
       }elsif($read_site_m == ${${$loxp_region}{$region_index}}[1]){
       	  $node_side = ${${$loxpnod}{$region_index}}[0];    
          print "Warning! you make a mistake! be carefull the bug1...\n";
       }elsif($read_site_l==${${$loxp_region}{$region_index}}[1]){
       	  $node_side = ${${$loxpnod}{$region_index}}[0]; 
       }elsif($read_site_l==${${$loxp_region}{$region_index}}[2]){
       		$node_side = ${${$loxpnod}{$region_index}}[1];
       		print "Warning! you make a mistake! be carefull the bug2...\n";             
       }
	 }
         
 return  $node_side;	 
}


####get region index####

sub regionindex{
	 my ($coordinate,$pointindex)=@_;
	 my $loxpindex = 0;
	 if(exists ${$pointindex}{$coordinate}){
	 	 $loxpindex = ${$pointindex}{$coordinate};
	 }else{print "Warning! you make a mistake! be carefull the bug5...\n";}
	 return $loxpindex;
	}

###input cvg file###
sub readcvg{
	 my ($cvg,$chrid) = @_;
   open CVG,$cvg or die $!;   
   $/ =">";
    <CVG>;
   my %chrdepth;

   while(<CVG>){
	    chomp;
	    next if ($_ =~ /^#/);
	    my $depth = ();
	    my @line = split /\n/;
	    my $chrID = shift @line;
	    next if ($chrID ne $chrid);
	    $depth = join (" ",@line);       
 
      @{$chrdepth{$chrid}} = split /\s+/,$depth;
      my $deplength =  @{$chrdepth{$chrid}};      
      #print "$chrid $deplength\n";
      last;
   }
   return %chrdepth;
   $/ ="\n";
   close CVG;
}
###stat nomal mapping loxp reads on the breakpoint#######
sub statnomalnod{
  my ($chrID,$nomal_node,$soap,$soap2)=@_;
  my %nomal_node_stat;
  my %stat_region;
  foreach (keys %{$nomal_node}){
    $nomal_node_stat{$_} = 0;
	  my $site = (split /\t/,${$nomal_node}{$_})[0];
	  my $stat_region_left = $site - 67;
	  my $stat_region_right =$site - 31;
	  push @{$stat_region{$_}},($stat_region_left,$stat_region_right);
	}
	
	if($soap =~ /gz$/){open SOP,"gzip -dc $soap |" or die $!;}else{open SOP, "$soap" ||die $!;}
  
	if($soap2 =~ /gz$/){open SOP2,"gzip -dc $soap2 |" or die $!;}else{open SOP2, "$soap2" ||die $!;}
	
	while (<SOP>){
		 chomp;
	   my ($chrid,$site)= (split /\s+/,$_)[7,8];
	   next if ($chrid ne $chrID);
	    
	   foreach (keys %stat_region){
	   	$nomal_node_stat{$_} ++ if ($site >= ${$stat_region{$_}}[0]&&$site <= ${$stat_region{$_}}[1]);
	   }
	 }
	 
	while (<SOP2>){
		 chomp;
	   my ($chrid,$site)= (split /\s+/,$_)[7,8];
	   next if ($chrid ne $chrID);
	    
	   foreach (keys %stat_region){
	   	$nomal_node_stat{$_} ++ if ($site >= ${$stat_region{$_}}[0]&&$site <= ${$stat_region{$_}}[1]);
	   }
	 } 
	return (\%nomal_node_stat);	
}
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
  	print "$chrID\t$star\t$END\t0\n";
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
      	$findex = (split /_/,$_)[0];      	 
        push @{$newkey{$findex}},$_;
      }
   }elsif($opern==0){ 
  
   foreach (@keyid){
      	$findex = (split /\//,$_)[0];      	 
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
#####estimate CN by x minimum Masdist(Zscore)####
=head
sub Masdcn{
 	my ($cninfo,$opern) = @_;
 	my @expectation = (1.017,2.005,3.006,3.869);
 	my @stddev = (0.163200245097855,0.290479465711434,0.392868527627246,0.532309984877233);
 	my $masdist_min = (abs($cninfo-$expectation[0]))/($stddev[0]);
 	my $getcn = 1;
 	for(my $i=1;$i<4;$i++){
   	my $masdist = (abs($cninfo-$expectation[$i]))/($stddev[$i]);
   	if($masdist_min>$masdist){
   		$masdist_min = $masdist;
   		$getcn = $i + 1;  		
   	}
  }
  return($getcn);
}
=cut 	

#####get copy number(CN) and coefficient of variation(CV)#####
sub getCN_CV{
	my ($CNinfo,$opern) = @_;
	my ($region_CN,$CN_CV)=(0,0);
  my %eachCN;
  my ($totalCN,$n) = (0,0);
  my %tempnum = ();
  foreach(keys %{$CNinfo}){
    if(${$CNinfo}{$_} eq "--"){
    	  $eachCN{$_} = "--";
    	  next;}
    $eachCN{$_} = int (${$CNinfo}{$_} + 0.4); #Round off
#    if($cntemp>=5){
#    	$eachCN{$_} = int (${$CNinfo}{$_} + 0.5);
#    }else{$eachCN{$_} = Masdcn(${$CNinfo}{$_},1); }		
    $region_CN = $eachCN{$_} if ($eachCN{$_} > $region_CN);
    if (exists $tempnum{$eachCN{$_}}){$tempnum{$eachCN{$_}} ++;}else{$tempnum{$eachCN{$_}}=1;} 
    $totalCN += $eachCN{$_};
    $n++;    
  }
  if($opern==2){
  	my $surportnum=0;
  	foreach (keys %tempnum){
  		  if ($surportnum <$tempnum{$_}){
  	       $region_CN = $_;
  	       $surportnum = $tempnum{$_}; 
  	   }
  	   
  	}  
  }
   ###caculate the coefficient of variation(CV)###
  my $averCN = $totalCN/$n;
  my $square_deviation = 0;
    foreach (keys %eachCN){ 
       next if($eachCN{$_} eq "--");
       $square_deviation +=  ($eachCN{$_}-$averCN)*($eachCN{$_}-$averCN);
     }
  my $stdCN = sqrt($square_deviation/($n-1));
  $CN_CV =  $stdCN/$averCN;    
  return ($region_CN,$CN_CV,%eachCN);
 }
######get copy number(CN) and Masdist-test#####
sub CN_Mascheck{
	my ($CNinfo,$opern) = @_;
	
  my ($region_CN,$P_all,$CN_P);
  my @CNarray;
  
  foreach(keys %{$CNinfo}){
    if(${$CNinfo}{$_} eq "--"){
    	  next;
    }else{push @CNarray,${$CNinfo}{$_};}
  }
  my($best_CN,$min_p_value,@P_value) = Mascheck (@CNarray); #  Ttestcheck.pm only can test the CN <=4 precisely. 	
  return ($best_CN,$min_p_value,@P_value);
} 
######caculate copy number####
sub caculateCN{
	 #opern = 1 take the CN with the maximum value and estimate by CV.
	 #opern = 2 take the CN with the maximum surport and estimate by CV.
	 #opern = 3 take the CN with the minimum mahalanobis distance and estimate by misjudgement probability. 
   
   my ($testee_depths,$ref_depth,$testee_CN,$ref_CN,$regionkeys,$opern) = @_;
   
   my (%region_CN,%CN_CV,%CN_P);
   
   print "Edge: ";
   foreach (@{$regionkeys}){print "$_\t";}
   print "\nCN_F: ";
   foreach (@{$regionkeys}){print "${$testee_CN}{$_}\t";}
   
   if($opern != 3){ print "CN\tCV\n";
    }else{print "P(1)\tP(2)\tP(3)\tP(4)\tCN\tP_value\n"};
      
   foreach my $i(@{$regionkeys}){
   	  my %copy_number;
   	  foreach my $j(@{$regionkeys}){
   	  	  if ($i eq $j){
   	  	  $copy_number{$j} = "--"; next;}
   	      my $copy_num = (${$ref_depth}{$i}/${$testee_depths}{$j})*(${$testee_depths}{$i}/${$ref_depth}{$i})*(${$ref_CN}{$i}/${$ref_CN}{$j})*${$testee_CN}{$j};
   	      $copy_num =  sprintf("%.2f", $copy_num);
   	      $copy_number{$j} = $copy_num;  
      }
      #####get copy number#####
      my %eachCN; 
      
      if($opern != 3){
      	print "$i\t";
        ($region_CN{$i},$CN_CV{$i},%eachCN) = getCN_CV(\%copy_number,$opern);#opern = 2 take the CN with  maximum surport
        foreach(@{$regionkeys}){
      	   print "$eachCN{$_}\t"; 
         }
          print "$region_CN{$i}\t$CN_CV{$i}\n";  
      }else{
      	my @P_all;     	
       	($region_CN{$i},$CN_P{$i},@P_all) = CN_Mascheck(\%copy_number,3); 
        if(${$testee_CN}{$i} > 4){$region_CN{$i}=${$testee_CN}{$i};$CN_P{$i}="NA";}
       	print "$i\t";     	
        foreach(@{$regionkeys}){
      	   print "$copy_number{$_}\t"; 
         }         
          print "$P_all[0]\t$P_all[1]\t$P_all[2]\t$P_all[3]\t$region_CN{$i}\t$CN_P{$i}\n";        
      }
     
   }

   if($opern != 3){
     return (\%region_CN,\%CN_CV);       	
   }else{
   	 return (\%region_CN,\%CN_P);
   }
}


#####reckon copy number####

sub reckonCN{
	 my ($testee_depth,$ref_depth) = @_;
   my %edge_CN;
   
   foreach (keys %{$testee_depth}){
   	print "${$testee_depth}{$_}\n";
  }
   my @regionkeys = sortkey(\%{$testee_depth},1);	
   
   print "step0:Collect depth information...\n";
   print "Refseq: ";
   foreach (@regionkeys){print "$_\t";}
   print "\nAverD:  ";
   foreach (@regionkeys){print "${$ref_depth}{$_}\t";}
   print "\n##\n";
   print "Testee: ";
   foreach (@regionkeys){print "$_\t";}
   print "\nAverD:  ";   
   foreach (@regionkeys){
      my @infodepth = split /\t/,${$testee_depth}{$_};
      print "$infodepth[0]\t";      
      ###remove copy number(CN) = 0 region###  
      ###need refine....
      if ($infodepth[0] <= $mdep){
      	 $edge_CN{$_} = "0";
      	 delete ${$testee_depth}{$_};
      	 delete ${$ref_depth}{$_};
      	}
   }
   
   print "\n\nstep1:Remove the region of CN = 0...\n";
   my %testee_depths;
   @regionkeys=();
   @regionkeys = sortkey(\%{$testee_depth},1);	
   print "refseq: ";
   foreach (@regionkeys){print "$_\t";}
   print "\nAverD:  ";
   foreach (@regionkeys){print "${$ref_depth}{$_}\t";}
   print "\n##\n";
   print "Testee: ";
   foreach (@regionkeys){print "$_\t";}
   print "\nAverD:  "; 
   foreach (@regionkeys){
   	  my @infodepth = split /\t/,${$testee_depth}{$_};
   	  print "$infodepth[0]\t";   
      $testee_depths{$_} = $infodepth[0];        
   }
   print "\n\n";
   print "step2: CN initialization...\n";
   print " #Suppose the CN of each region is 1.\n"; 
   
   my (%ref_CN,%testee_CN);
   foreach (@regionkeys){
   	  $ref_CN{$_} = 1;
   	  $testee_CN{$_} = 1;
   	}
   #coefficient of variation(CV)	

   my ($region_CN1,$CN_CV1) = caculateCN(\%testee_depths,\%{$ref_depth},\%testee_CN,\%ref_CN,\@regionkeys,1);
  
   print "\nstep3: CN iterative estimate...\n";
   print "\n[1] CN estimate cycle 1...\n";
   my ($region_CN2,$CN_CV2) = caculateCN(\%testee_depths,\%{$ref_depth},\%{$region_CN1},\%ref_CN,\@regionkeys,2);
   
   my ($region_CNt,$CN_CVt) = ($region_CN2,$CN_CV2);
   my ($region_CNF,$CN_CVF);
   
   foreach my $cycle (2..$mcycle){         
      print "\n[$cycle] CN estimate cycle $cycle...\n";      
 
      ($region_CNF,$CN_CVF) = caculateCN(\%testee_depths,\%{$ref_depth},\%{$region_CNt},\%ref_CN,\@regionkeys,2);   
      
      my $ck = 0;
      foreach (keys %{$region_CNF}){
       $ck ++ unless ${$region_CNF}{$_} == ${$region_CNt}{$_}; 
      }
      last if $ck == 0; 
      
      undef %{$region_CNt};
      undef %{$CN_CVt};      
      ($region_CNt,$CN_CVt) = ($region_CNF,$CN_CVF);
    
    }
   
   print "\nstep4:check CN by mahalanobis distance...\n"; 
   my %edgeinfo;
   
   foreach (keys %{$region_CNF}){
   	$edgeinfo{$_} = "$_\t${$region_CNF}{$_}\t$constant_depth{$_}";
   }
   
   my @cns = cnstat ($chrID,$cutlen,\%edgeinfo,$coverage,$refcover,$loxpregion);
     
   my %finalcn;
   $/ = "\n";
   foreach (@cns){
    	chomp;
      my @CNstats = split /\s+/,$_;
     	my $edgeid = shift @CNstats;
   	  push @{$finalcn{$edgeid}},(@CNstats);
   }
  
   my (%region_CN,%CN_P,@P_all);
   my @finalkey;
   foreach (keys %{$region_CNF}){
   	 if(${$region_CNF}{$_} > 4){
   	 	$edge_CN{$_} =${$region_CNF}{$_};$CN_P{$_}="NA";
   	}else{push @finalkey,$_;}
   	}
   
   foreach(@finalkey){
      my($best_CN,$min_p_value,@P_value) = Mascheck (@{$finalcn{$_}});  
      print "$_ : ". join "\t",@P_value;
      print "\n";
      
      my $len = @{$finalcn{$_}};
      $edge_CN{$_} = $best_CN;
      $CN_P{$_} = sprintf("%.4f", $min_p_value);
   }
   return (\%edge_CN,\%CN_P);
 }

#########=============#######   
#########=============#######

=head1 Name
  
  extract_node_and_edge.pl--This script is used to stat the edges/nodes infomation
                             and parse the input for PATH . 
   
=head1 Description
   
  This script is only used for stat the edges/nodes infomation and parse the input 
  for PATH .      
  
  Theoretically, the copy number of each region in a reference genome can be 
  reflected from the sequencing depth of corresponding region. However, the fact 
  is copy number cannot be transformed from sequencing depth directly because of 
  random bias caused by sequencing. Therefore, we choose a iterative algorithm to 
  refine the copy number estimation. 
  
  Mahalanobis Distance(MD) analysis was introduced to evaluate our method and 
  calculate misjudgment value p. For each edge, with the estimated copy number, the
  edge was further divided into several sub regions and calculate the corresponding
  copy number of each sub region. For each sub region¡¯s copy number, The MD 
  analysis was used to evaluate whether it is correct or not. 
     
=head1 Version
  
  Author: Yun Wang, wangyun\@genomic.org.cn
  
  Version: 1.0, Date: 2012-11-05
   
=head1 Usage
  
  perl extract_node_and_edge.pl [option]
        -prefix     :   prefix of output file [default: parsepath] 
        -chrid      :   synthetic chromosome ID       
        -coverage   :   soapcoverage file of sample
        -refcover   :   soapcoverage file of reference sample(unscrambled sample)
        -abnorm     :   abnormal mapped split reads (*.sort.gvr)  
        -soap       :   pair_end soap result including normal mapped loxp reads
        -soap2      :   single_end soap result including normal mapped loxp reads
        -loxpregion :   loxpregion file
        -minread    :   minimal split reads for supporting a breakpoint [default: 2]
        -mdep       :   maximal average sequencing deth of copy number = 0 [default: 10] 
        -mcycle     :   maximal cycle for estimate copy number [default: 10] 
        -cutlen     :   cut each region into subregion for Mahalanobis Distance(MD) analysis [default: 500] #500 model only used
        -outdir     :   output directory [default: ./]
        -help       :   show help message
  
=head1 Example
  
  perl extract_node_and_edge.pl -prefix <str> -chrid <str> -coverage <coverage_file>
       -refcover <refcover_file> -abnorm <sort_gvr_file> -soap <soap_file> -soap2 
       <soap2_file> -loxpregion <loxp_regio_file> [option]
     
=cut       


