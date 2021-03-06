#! /usr/bin/perl -w
use strict;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use List::Util qw/max min/;
use PerlIO::gzip;
use Cwd;

my $usage=<<USAGE;

Describe: filter reads generated by Hiseq; input file format: fq or fq.gz...
       
         
  Author: Yun Wang, wangyun\@genomic.org.cn
 Version: 1.0.0, Date: 2015-04-11

Usage:perl $0 <fq1> <fq2> -prefix <str> -o <output_path>
          
      -prefix     <str>  : Sample ID/prefix of output file     
      -outdir            : Directory of output [default: ./]   
      -help              : show this message
      
      -minlen            : minimal read lenght allowed [default 100]
      -tfPhred           : transformation value from ASCII to Phred score [default 64], for Miseq and Hisq2500, this parameter should be 33.   
      -lowPhred          : maximal phred score of base which is identified as a low quality base [default 10]
      -lowRate           : maximal ratio of low quality base in a read [default 0.01] 
      -Ns                : maximal number of N base arrowed [default 2]
      
      -trim              : trim to abtain a eligible read [default None]
      -seed              : seed length of trimed [default 50]
      
      -mvbase            : trim both END of read [default None]
      -star              : trim the number of left END [default 0]
      -end               : trim the number of right END [default 0]
      
      -deadapter         : remove adapter reads [operation default None]
      -adapterlist1      : list of adpterlist1
      -adapterlist2      : list of adpterlist2
      
      -outSE             : output for single-end reads [operation default None]     
      -outfilted         : output for abandoned reads [operation default None] #Not including reads with adapter sequence...
      
      
Example:
   1 filter reads ...
     perl $0 <*fq1> <*fq2> -prefix <str>  -minlen [int]  -lowPhred [int]  -lowRate [float] -Ns [int] -outdir <output_path>  
   1 filter reads by trim parameter...
     perl $0 <*fq1> <*fq2> -prefix <str>  -minlen [int]  -lowPhred [int]  -lowRate [float] -Ns [int] -trim -seed [int] -outdir <output_path>  
     
Note:  Only deal with pair-end reads... 

USAGE

my ($prefix,$minlen,$tfPhred,$lowPhred,$lowRate,$Ns,$trim,$seed,$mvbase,$mvbase_st,$mvbase_end,$deadapter,$adapterlist1,$adapterlist2,$outSE,$outfilted,$outdir,$help);
GetOptions(
   "prefix:s" => \$prefix,
   "minlen:s" => \$minlen,
   "tfPhred:s" => \$tfPhred,
   "lowPhred:s" => \$lowPhred,
   "lowRate:s" => \$lowRate,
   "Ns:s" => \$Ns,
   
   ##
   "trim" => \$trim,
   "seed:s" => \$seed,
   
   "mvbase" => \$mvbase,
   "star:s" => \$mvbase_st,
   "end:s" => \$mvbase_end,
   
   "deadapter" => \$deadapter,
   "adapterlist1:s" => \$adapterlist1,
   "adapterlist2:s" => \$adapterlist2,
   
   "outSE" => \$outSE,
   "outfilted" => \$outfilted, 
   
   "outdir:s"  => \$outdir,   
   "help"  => \$help,
);

die $usage if (@ARGV != 2 || ($deadapter && (!$adapterlist1 || !$adapterlist2) ||$help));

$prefix ||= basename ($ARGV[0]);
$minlen ||= 100;
$tfPhred ||= 64;
$lowPhred ||= 10;
$lowRate ||= 0.01;
$Ns ||= 2 unless ($Ns == 0);
$seed ||= 50;
$mvbase_st ||= 0;
$mvbase_end ||= 0;
$outdir ||= getcwd;

my (%mvrdsid1,%trmid1);
my (%mvrdsid2,%trmid2);

my %deadapterid;
if(defined $deadapter){
   foreach($adapterlist1,$adapterlist2){
  	 if ($_ =~ /.gz$/){
	      open IN, "gzip -dc $_ |" or die $!;
     }else{
        open IN, $_ or die $!;
     }
     
     while (<IN>){
        chomp;
        next if $_=~/^#/;
        my $adapterid = (split /\s+/,$_)[0];       
        $deadapterid{$adapterid} = 1
     }
     close IN;     
  }
}


if (defined $trim){
   (%mvrdsid1,%trmid1) = identify_rds($ARGV[0]); 
   (%mvrdsid2,%trmid2) = identify_rds($ARGV[1]);   
}else{
   (%mvrdsid1) = identify_rds($ARGV[0]); 
   (%mvrdsid2) = identify_rds($ARGV[1]);
}

open my $pe1, ">:gzip","$outdir/$prefix\_1.clean.fq.gz" ||die $!; 
open my $pe2, ">:gzip","$outdir/$prefix\_2.clean.fq.gz" ||die $!; 
open STT, ">$outdir/$prefix.clean.stat" ||die $!;

open SE, ">:gzip","$outdir/$prefix.SE.clean.fq.gz" ||die $! if defined $outSE; 
open AB, ">:gzip","$outdir/$prefix.AB.clean.fq.gz" ||die $! if defined $outfilted; 

my $sdstr;

if (defined $trim){
  $sdstr .= "1" foreach (1..$seed);
}

my (%stat1,%stat2);


%stat1 = filter_rds($ARGV[0],$pe1);
%stat2 = filter_rds($ARGV[1],$pe2);


my $raw_Q20_ratio1 = $stat1{raw_Q20}/$stat1{raw_base};
$raw_Q20_ratio1 = sprintf("%.3f",$raw_Q20_ratio1);

my $raw_Q30_ratio1 = $stat1{raw_Q30}/$stat1{raw_base};
$raw_Q30_ratio1 = sprintf("%.3f",$raw_Q30_ratio1);

my $clean_Q20_ratio1 = $stat1{clean_Q20}/$stat1{clean_base};
$clean_Q20_ratio1 = sprintf("%.3f",$clean_Q20_ratio1);

my $clean_Q30_ratio1 = $stat1{clean_Q30}/$stat1{clean_base};
$clean_Q30_ratio1 = sprintf("%.3f",$clean_Q30_ratio1);


my $raw_Q20_ratio2 = $stat2{raw_Q20}/$stat2{raw_base};
$raw_Q20_ratio2 = sprintf("%.3f",$raw_Q20_ratio2);

my $raw_Q30_ratio2 = $stat2{raw_Q30}/$stat2{raw_base};
$raw_Q30_ratio2 = sprintf("%.3f",$raw_Q30_ratio2);

my $clean_Q20_ratio2 = $stat2{clean_Q20}/$stat2{clean_base};
$clean_Q20_ratio2 = sprintf("%.3f",$clean_Q20_ratio2);

my $clean_Q30_ratio2 = $stat2{clean_Q30}/$stat2{clean_base};
$clean_Q30_ratio2 = sprintf("%.3f",$clean_Q30_ratio2);

print STT "#Raw Data\n";
print STT "Raw reads number of read1: $stat1{raw_rdsnum}\n";
print STT "Raw reads number of read2: $stat2{raw_rdsnum}\n";
print STT "Raw base pair of read1: $stat1{raw_base}\n";
print STT "Raw base pair of read2: $stat2{raw_base}\n";
print STT "Q20 base Ratio in raw read1 data: $raw_Q20_ratio1\n";
print STT "Q30 base Ratio in raw read1 data: $raw_Q30_ratio1\n";
print STT "Q20 base Ratio in raw read2 data: $raw_Q20_ratio2\n";
print STT "Q30 base Ratio in raw read2 data: $raw_Q30_ratio2\n";

print STT "\n#Clean Data\n";
print STT "Clean reads number of read1: $stat1{clean_rdsnum}\n";
print STT "Clean reads number of read2: $stat1{clean_rdsnum}\n";
print STT "Clean base pair of read1: $stat1{clean_base}\n";
print STT "Clean base pair of read2: $stat2{clean_base}\n";
print STT "Q20 base Ratio in clean read1 data: $clean_Q20_ratio1\n";
print STT "Q30 base Ratio in clean read1 data: $clean_Q30_ratio1\n";
print STT "Q20 base Ratio in clean read2 data: $clean_Q20_ratio2\n";
print STT "Q30 base Ratio in clean read2 data: $clean_Q30_ratio2\n";

my $filtRate_rds1 = 1 - $stat1{clean_rdsnum}/$stat1{raw_rdsnum};
$filtRate_rds1 = sprintf("%.3f",$filtRate_rds1);
my $filtRate_base1 = 1 - $stat1{clean_base}/$stat1{raw_base};
$filtRate_base1 = sprintf("%.3f",$filtRate_base1);

my $filtRate_rds2 = 1 - $stat1{clean_rdsnum}/$stat1{raw_rdsnum};
$filtRate_rds2 = sprintf("%.3f",$filtRate_rds2);
my $filtRate_base2 = 1 - $stat1{clean_base}/$stat1{raw_base};
$filtRate_base2 = sprintf("%.3f",$filtRate_base2);

print STT "\n#Filtering statistics\n";
print STT "Filtered read1 by read length < $minlen: $stat1{minlen}\n";
print STT "Filtered read1 by low-quality base[phred score < $lowPhred] >$lowRate: $stat1{lowquality}\n";
print STT "Filtered read1 by Ns >$Ns: $stat1{Ns}\n";
print STT "Filtered read1 by pair-end: $stat1{nopair}\n";
print STT "Filtered read1 by adapter: $stat1{adapter}\n";

print STT "Filtering ratio of read1 by read: $filtRate_rds1\n";
print STT "Filtering ratio of read1 by base: $filtRate_base1\n";
print STT "Trimed read of read1: $stat1{trim}\n" if defined $trim;
print STT "Trimed base of read1: $stat1{trimbase}\n" if defined $trim;

##
print STT "\nFiltered read2 by read length < $minlen: $stat2{minlen}\n";
print STT "Filtered read2 by low-quality base[phred score < $lowPhred] >$lowRate: $stat2{lowquality}\n";
print STT "Filtered read2 by Ns >$Ns: $stat2{Ns}\n";
print STT "Filtered read2 by pair-end: $stat2{nopair}\n";
print STT "Filtered read2 by adapter:$stat2{adapter}\n";

print STT "Filtering ratio of read2 by read: $filtRate_rds2\n";
print STT "Filtering ratio of read2 by base: $filtRate_base2\n";
print STT "Trimed read of read2: $stat2{trim}\n" if defined $trim;
print STT "Trimed base of read2: $stat2{trimbase}\n" if defined $trim;
##

sub filter_rds
{
  my ($fqfile,$outfd) = @_;
  
  if ($fqfile =~ /fq.gz$/){
     open IN, "gzip -dc $fqfile |" or die $!;
  }else{
     open IN, $fqfile or die $!;
  }
  
  my %stat; 
  $stat{raw_rdsnum} = 0;
  $stat{raw_base} = 0;
  $stat{raw_Q20} = 0;
  $stat{raw_Q30} = 0;

  $stat{clean_rdsnum} = 0;
  $stat{clean_base} = 0;
  $stat{clean_Q20} = 0;
  $stat{clean_Q30} = 0;

  $stat{minlen} = 0;
  $stat{lowquality} = 0;
  $stat{Ns} = 0;
  $stat{nopair} = 0;
  $stat{adapter} =0;
    
  $stat{trim} = 0;
  $stat{trimbase} = 0;

  my $rdslen;
  
  $/ = "\n@";
  while (<IN>){
   	my $rdsinfo = $_;
    chomp $rdsinfo;
    $stat{raw_rdsnum} ++;
    $rdsinfo =~ s/^@// if  $stat{raw_rdsnum} ==1;
    
    my @fqinfo = split /\n/,$rdsinfo;
    
    if (@fqinfo < 4) {
      my $addinfo = <IN>;
      chomp $addinfo;
      $addinfo = "\@".$addinfo;
      my @fqaddinfo = split /\n/,$addinfo;
      push @fqinfo, @fqaddinfo;
    }
  
    print "Warning:a error for read $fqinfo[0]\n" unless @fqinfo == 4;
     
    my @base = split //,$fqinfo[1];
    my @qual = split //,$fqinfo[3];
    
    if(defined $mvbase){
      if($mvbase_st){
        foreach (1..$mvbase_st){
          shift @base;
          shift @qual;
        }
      }
        
      if($mvbase_end){
        foreach (1..$mvbase_end){
          pop @base;
          pop @qual;
        }
      }
      $fqinfo[1] = join "",@base;
      $fqinfo[3] = join "",@qual;
    }
    
    $fqinfo[0] =~ /(.+)\/[12]/;

    my $rdsid1 = "$1"."/1";
    my $rdsid2 = "$1"."/2";
    
    ##check read length
    $rdslen = length($fqinfo[1]);
    $stat{raw_base} += $rdslen;
    if ($rdslen < $minlen){
      $stat{minlen} ++;
      next;
    }elsif(exists $deadapterid{$fqinfo[0]}){
       $stat{adapter} ++;
     	 next;
    }

    my ($Nsnum,$lowqual) = (0,0);
    my $qualstr = "";
    $Nsnum ++ while ($fqinfo[1] =~ /N/gi);

    my $i =0;
    foreach (@qual)
    {
       my $phred = ord($_) - $tfPhred;
       
       if ($phred >= 20){
        	$stat{raw_Q20} ++;
       }
       
       if ($phred >= 30){
          $stat{raw_Q30} ++;
       }
           
       if ($phred < $lowPhred){
          $lowqual ++;
       }
    }
    my $lowqualrate = $lowqual/$rdslen;
        
    if ($rdslen >= $minlen && $Nsnum <= $Ns && $lowqualrate  <= $lowRate){
       
       unless(exists $mvrdsid1{$rdsid1} || exists $mvrdsid2{$rdsid2}){
         print $outfd "@".join "\n",@fqinfo;
         print $outfd "\n";
      
         $stat{clean_base} += $rdslen;
         $stat{clean_rdsnum} ++;
  
         foreach (@qual)
         {
           my $phred = ord($_) - $tfPhred;
       
           if ($phred >= 20){
             $stat{clean_Q20} ++;
           }
           if ($phred >= 30){
             $stat{clean_Q30} ++;
           }
         }
      }else{
        $stat{nopair} ++;
        if(defined $outSE){
          print SE "@".join "\n",@fqinfo;
          print SE "\n";
        }else{
          
          if(defined $outfilted){
             print AB "@".join "\n",@fqinfo;
             print AB "\n";
           }
        }
     }
   }else{
      if (defined $trim){ 
         my $trminfo;
         if(exists $trmid1{$fqinfo[0]}){ 
         	  $trminfo = $trmid1{$fqinfo[0]};
         	}elsif(exists $trmid2{$fqinfo[0]}){
         		$trminfo = $trmid2{$fqinfo[0]};
         	}
          
         if(defined $trminfo){
      	    @fqinfo = split /\n/,$trminfo;      	    
      	    my @base = split //,$fqinfo[1];
            my @qual = split //,$fqinfo[3];
            
            unless(exists $mvrdsid1{$rdsid1} || exists $mvrdsid2{$rdsid2}){
               print $outfd "@".join "\n",@fqinfo;
               print $outfd "\n";
                
               my $trimrdslen = length ($fqinfo[0]);
               $stat{clean_base} += $trimrdslen;
               $stat{clean_rdsnum} ++;
               
               $stat{trim} ++;
               $stat{trimbase} += $rdslen - $trimrdslen;
               
                foreach (@qual)
                {
			             my $phred = ord($_) - $tfPhred;
			 
			             if ($phred >= 20){
			    	         $stat{clean_Q20} ++;
			             }
		  	           if ($phred >= 30){
		  	             $stat{clean_Q30} ++;
		  	           }
                }
             }else{
             	  $stat{nopair} ++;
        	      if(defined $outSE){
         	         print SE "@".join "\n",@fqinfo;
        	         print SE "\n";
        	      }else{
        	         if (defined $outfilted){
        	           print AB "@".join "\n",@fqinfo;
        	           print AB "\n";
        	         }
        	      }
             } 
      
         }else{
         	
         	  if ($Nsnum > $Ns){
         	  	$stat{Ns} ++;
         	  }elsif($lowqualrate  > $lowRate){
         	    $stat{lowquality} ++;
         	  }
         	  
            if (defined $outfilted){
        	     print AB "@".join "\n",@fqinfo;
        	     print AB "\n";
      	    } 
         }
      }else{
      	 if ($Nsnum > $Ns){
         	  	$stat{Ns} ++;
          }elsif($lowqualrate  > $lowRate){
         	    $stat{lowquality} ++;
         	}
         	  
         if (defined $outfilted){
        	   print AB "@".join "\n",@fqinfo;
        	   print AB "\n";
      	 }
      }
   } 
 }
 
 close IN;
 $/= "\n";
 return  %stat; 
      
}

##
sub identify_rds
{  my $fqfile = $_[0];
	 my (%mvrdsid,%trmrds);
	 
	 if (defined $trim){
      $sdstr .= "1" foreach (1..$seed);
   }
	 
	 if ($fqfile =~ /fq.gz$/){
	   open IN, "gzip -dc $fqfile |" or die $!;
   }else{
     open IN, $fqfile or die $!;
   }
   
   $/ = "\n@";
   my $n = 0;
   while (<IN>){
   	 my $rdsinfo = $_;
     chomp $rdsinfo;
     $n ++;
     $rdsinfo =~ s/^@// if $n ==1;
     
     my @fqinfo = split /\n/,$rdsinfo;
     
     if (@fqinfo < 4) {
       my $addinfo = <IN>;
       chomp $addinfo;
       $addinfo = "\@".$addinfo;
       my @fqaddinfo = split /\n/,$addinfo;
       push @fqinfo, @fqaddinfo;
     }
  
     print "Warning:a error for read $fqinfo[0]\n" unless @fqinfo == 4;
  
     my @base = split //,$fqinfo[1];
     my @qual = split //,$fqinfo[3];
     
     ##check read length
     my $rdslen = length($fqinfo[1]);

     if ($rdslen < $minlen){
	     $mvrdsid{$fqinfo[0]} = 1;
     	 next;
     }elsif(exists $deadapterid{$fqinfo[0]}){
       $mvrdsid{$fqinfo[0]} = 1;
     	 next;
     }
     
     my ($Nsnum,$lowqual) = (0,0);
     my $qualstr = "";

     my $i =0;
     foreach (@qual)
     {
			 my $phred = ord($_) - $tfPhred;
			 my $qualstr_tmp;
			 if ($phred < $lowPhred){
			 	  $lowqual ++;
			 	  $qualstr_tmp = "0"; 
			 }else{
			 	  $qualstr_tmp = "1"; 
			 }
			 if($base[$i] eq "N"){
			    $Nsnum ++;
			    $qualstr_tmp = "#"; 
			 }			 
			 $qualstr .= "$qualstr_tmp"; 
	     $i++;
    }

    my $lowqualrate = $lowqual/$rdslen;


   if($rdslen >= $minlen && $Nsnum <= $Ns && $lowqualrate  <= $lowRate){
     next;  
   }else{
     if(defined $trim){

       if ($qualstr =~ m/$sdstr/ig){
      	  my @str_qual = split //,$qualstr  ;
      	  my ($star,$ends);
          $star = index($qualstr,"$sdstr");
          ##elonge... 
          $ends = $star + $seed;
          
          while (1){
            if(defined $str_qual[$ends] && $str_qual[$ends] eq "1"){
            	 if($ends < $rdslen - 1){
            	 $ends ++;	
               }else{
        	     last;
               }    
            }else{
               $ends --;
               last;
            }
         }
      
         my ($star_tp,$ends_tp) = ($star,$ends);  
      
         while(1){
           	my ($star_tmp,$ends_tmp) = ($star_tp,$ends_tp);
      	    #Up elongation...
            if($star_tp > 1){
       	      $star = $star_tp - 2;
       	  
       	      while (1){      	 
                 if($str_qual[$star] eq "1"){
          	        if($star>0){
          		      $star --;
                    }else{
           	    	  last;
           	        }           	
                 }else{
             	      $star ++;
        		        last;
                 }
              }           	  
       	   }else{
       		     $star = $star_tp;
       	   }        
           my $chc_up_elong = checkstr($star,$ends_tp,$qualstr);
       
           #Down elongation...
           if ($ends_tp < $rdslen - 2){
        	    $ends = $ends_tp + 2;
       	      while (1){       	 
                 if($str_qual[$ends] eq "1"){
          	       if($ends< $rdslen - 1){
          		     $ends ++;
                   }else{
           		     last;
           	       }           	
                 }else{
             	     $ends --;
        		       last;
                 }
              }              
           }else{
              $ends = $ends_tp;
           }
           my $chc_dn_elong = checkstr($star_tp,$ends,$qualstr);
           my $checkout = checkstr($star,$ends,$qualstr);
    
           #print "#$star,$ends\t$checkout\t$chc_up_elong\t$chc_dn_elong\n";
           if($checkout){
             ($star_tp,$ends_tp) = ($star,$ends);  
           }elsif($chc_up_elong && $chc_dn_elong){
             my $len_up = $ends_tp - $star + 1;
             my $len_dn = $ends - $star_tp + 1;
             if($len_up >= $len_dn){
               ($star_tp,$ends_tp) = ($star,$ends_tp);            
             }else{
               ($star_tp,$ends_tp) = ($star_tp,$ends);  
             }         
          }elsif ($chc_up_elong){
            ($star_tp,$ends_tp) = ($star,$ends_tp);  
          }elsif ($chc_dn_elong){
            ($star_tp,$ends_tp) = ($star_tp,$ends);  
          }else{  
            last;
          }

          if($star_tp == $star_tmp && $ends_tp == $ends_tmp){
            last;
          }
        
         }
        
       if($ends_tp - $star_tp + 1 >= $minlen){
          ##trimed reads results...
          my $star_idx = shift;
	        my $ends_idx = shift;
          my $rdstrim = shift;
          my $len = $ends_tp - $star_tp + 1;
          my $trim_len_L = $star_tp;
          my $trim_len_R = $rdslen - $ends_tp - 1;
                    
          $fqinfo[0] = "$fqinfo[0]"."trim$trim_len_L"."_$trim_len_R";
          $fqinfo[1] = substr($fqinfo[1],$star_tp,$len);
          $fqinfo[3] = substr($fqinfo[3],$star_tp,$len);
          $rdsinfo = join "\n",@fqinfo;              
          $trmrds{$rdsinfo} = $rdsinfo;
        
       }else{
          $mvrdsid{$fqinfo[0]} = 1;
       } 
      
      }else{
        $mvrdsid{$fqinfo[0]} = 1;
     }
      
   }else{
     $mvrdsid{$fqinfo[0]} = 1;
   }   
  }
  } 
  close IN;
  $/= "\n";
  
  if (defined $trim){
     return (%mvrdsid,%trmrds);
  }else{
     return (%mvrdsid);
  }
             
}



sub checkstr
{
	my $star_idx = shift;
	my $ends_idx = shift;
  my $strim = shift;
  
  my $chc_len = $ends_idx - $star_idx + 1;
  my $chc_strim = substr($strim,$star_idx,$chc_len);
  my ($chc_Ns,$chc_low) = (0,0);
  $chc_Ns ++ while ($chc_strim =~ /#/gi);	
  $chc_low ++ while ($chc_strim =~ /0/gi);
  
  my $chc_low_Rate = $chc_low / $chc_len;
 
 if($chc_Ns <= $Ns && $chc_low_Rate  <= $lowRate){ 
   return 1;
 }else{
   return 0;
 }

}
##
