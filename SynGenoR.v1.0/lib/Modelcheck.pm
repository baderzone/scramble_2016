package Modelcheck;
# This module is used to test copy number by T-test.
# The program will estimate 1-10 copy by T-test and return the P value and best CN and p value.   

use strict;
use lib;
use lib "/ifs1/ST/SYNBIO/USER/wangyun/Bin/lib";#lib path for List::Util.pm
use List::Util qw /min max/;
use Statistics::Distributions qw /tprob uprob/;
use Statistics::Basic qw(:all);
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

@EXPORT_OK = qw(Ttestcheck Utestcheck Zscorecheck Mascheck); 

$VERSION     = 1.00;
@ISA         = qw(Exporter);


sub Ttestcheck {
	my @CN = @_;
	#expectation ...
	my @P_value;
	my $min_p_value=1;
	my $best_CN=0;
	my @expectation = (1.017,2.005,3.006,3.869,5.000,6.000,7.000,8.000,9.000,10.000); 	
	#caculate average and stddev;
	my $mean   = mean(@CN);
	my $stddev   = stddev(@CN);
	my $n = @CN;
	my $df = $n-1;
	my $i = 1;
	 
	foreach my $E(@expectation ){
   	my $t = abs ($mean-$E)/($stddev/sqrt($n));
   	my $uper_prob = tprob ($df,$t);
    my $prob = 1-(2*$uper_prob);
   	$prob = sprintf("%.5f", $prob);
   	push @P_value, $prob;   	
	  if ($min_p_value > $prob){
	  		$min_p_value = $prob;
	  		$best_CN =$i;
	  	}
	  $i ++;
	}		
 return ($best_CN,$min_p_value,@P_value);

}

sub Utestcheck {
	my @CN = @_;
	#expectation ...
	my $pi = 3.1415926536;
	my @P_value;
	my $min_p_value=1;
	my $best_CN=0;
	my @expectation = (1.017,2.005,3.006,3.869);
	#my @stddev = (3.469,6.197,8.749,11.530);
  my @stddev = (0.163200245097855,0.290479465711434,0.392868527627246,0.532309984877233);
	
	#caculate average and stddev;
	my $mean   = mean(@CN);
	my $sstddev = stddev(@CN);
	my $n = @CN;
	#print "###PM test number :$n\t mean $mean\t \n";
	my $df = $n-1;
	 
	foreach my $i(0..3){
   	my $u = (abs ($mean-$expectation[$i]))/($stddev[$i]/sqrt($n));
   	#print "###PM U:E $expectation[$i]\t##$stddev[$i]\t##$sstddev\n";      	
   	my $uper_prob = uprob ($u);
    my $prob = 1-(2*$uper_prob);
   	$prob = sprintf("%.5f", $prob);
   	#print "###PM U:$u\t $prob\n";
   	push @P_value, $prob;   	
	  if ($min_p_value > $prob){
	  		$min_p_value = $prob;
	  		$best_CN =$i+1;
	  	}
	}		
 return ($best_CN,$min_p_value,@P_value);
}

sub Zscorecheck{
	my @CN = @_;
	#expectation ...
	my $pi = 3.1415926536;
	my @P_value;
	my $min_p_value=1;
	my $best_CN=0;
	my @expectation = (1.017,2.005,3.006,3.869);
  my @stddev = (0.163200245097855,0.290479465711434,0.392868527627246,0.532309984877233);
	
	#caculate average and stddev;
	my $mean   = mean(@CN);
	my $sstddev = stddev(@CN);
	my $n = @CN;
	print "###PM test number :$n\t mean $mean\t \n";
	my $df = $n-1;
	
	foreach my $i(0..3){
   	my $Z_score = (abs ($mean-$expectation[$i]))/($stddev[$i]);
   	print "###PM Z score:E $expectation[$i]\t##S $stddev[$i]\t## Z score $Z_score\n";      	
   	my $uper_prob = uprob ($Z_score);
    my $prob = 1-(2*$uper_prob);
   	$prob = sprintf("%.5f", $prob);
   	print "###PM Z_score :$Z_score\t $prob\n";
   	push @P_value, $prob;   	
	  if ($min_p_value > $prob){
	  		$min_p_value = $prob;
	  		$best_CN =$i+1;
	  	}
	}		
 return ($best_CN,$min_p_value,@P_value);
}

sub Mascheck {
	my @CN = @_;
	#expectation ...
	my $pi = 3.1415926536;
	my @P_value;
	my $best_CN=0;
	my @expectation = (1.017,2.005,3.006,3.869);
  my @stddev = (0.163200245097855,0.290479465711434,0.392868527627246,0.532309984877233);
	
	#caculate average and stddev;
	my $mean   = mean(@CN);
	my $sstddev = stddev(@CN);
	my $n = @CN;
	print "MD test number :$n\t mean $mean\t \n";
	my $df = $n-1;
  
  my ($sumd1,$sumd2,$sumd3,$sumd4) = 0;
=head  
  foreach my $x(@CN){
  	$sumd1 += Masd($x,$expectation[0],$stddev[0]); 
    $sumd2 += Masd($x,$expectation[1],$stddev[1]); 
    $sumd3 += Masd($x,$expectation[2],$stddev[2]);
    $sumd4 += Masd($x,$expectation[3],$stddev[3]);
   }
   my $dist1 = sqrt ($sumd1);
   my $dist2 = sqrt ($sumd2);
   my $dist3 = sqrt ($sumd3);
   my $dist4 = sqrt ($sumd4);
=cut
   my $dist1 = abs($mean-$expectation[0])/$stddev[0];
   my $dist2 = abs($mean-$expectation[1])/$stddev[1];
   my $dist3 = abs($mean-$expectation[2])/$stddev[2];
   my $dist4 = abs($mean-$expectation[3])/$stddev[3];

   
  #$best_CN = GetbestCN($sumd1,$sumd2,$sumd3,$sumd4); 
  $best_CN = GetbestCN($dist1,$dist2,$dist3,$dist4); 

  
###回代估计误判概率P（P_value）#####
   
   my $best_p_value = Misjudgment_prob($best_CN,\@CN,\@expectation,\@stddev);
   my $p_value1 = Misjudgment_prob(1,\@CN,\@expectation,\@stddev);
   my $p_value2 = Misjudgment_prob(2,\@CN,\@expectation,\@stddev);
   my $p_value3 = Misjudgment_prob(3,\@CN,\@expectation,\@stddev);
   my $p_value4 = Misjudgment_prob(4,\@CN,\@expectation,\@stddev);
   
   push @P_value,($p_value1,$p_value2,$p_value3,$p_value4);
   return ($best_CN,$best_p_value,@P_value);
   
####Sub Routine####   
  sub Masd {
  	my ($x,$u,$s) = @_;
  	my $dist = ((abs($x-$u))**2/$s**2);
  	#my $dist = (abs($x-$u)/$s);
  	return $dist
  	}
  	
  sub GetbestCN{
  	my $bestcn = 0; 
  	my ($d1,$d2,$d3,$d4) = @_;
  	
  	my $mindist = min ($d1,$d2,$d3,$d4); 
  	if($d1 == $mindist){
  		$bestcn = 1;
  	}elsif($d2 == $mindist){
  		$bestcn = 2;
  	}elsif($d3 == $mindist){
  		$bestcn = 3;
  	}elsif($d4 == $mindist){
  		$bestcn = 4;
  	  }
  	return $bestcn;   	
  	} 
  ###回代估计误判概率P（P_value）#####	
  sub Misjudgment_prob{
    my ($jdgmcn,$cn,$expe,$stddev) = @_; 
    my $Total_cn = @{$cn};
    my $true_jdg = 0;
    my $false_jdg = 0;
    my $misjdg_prob = 0;
    
    
    foreach (@{$cn}){
    	my $masd_cn1 = sqrt (Masd ($_,${$expe}[0],${$stddev}[0]));
      my $masd_cn2 = sqrt (Masd ($_,${$expe}[1],${$stddev}[1]));
      my $masd_cn3 = sqrt (Masd ($_,${$expe}[2],${$stddev}[2]));
      my $masd_cn4 = sqrt (Masd ($_,${$expe}[3],${$stddev}[3]));
      my $best_cn = GetbestCN($masd_cn1,$masd_cn2,$masd_cn3,$masd_cn4); 
      
      if ($best_cn == $jdgmcn){
      	$true_jdg ++;
      }else{
      	$false_jdg ++;
      }
    }
   $misjdg_prob = $false_jdg/$Total_cn;
   return ($misjdg_prob);
  }
}
