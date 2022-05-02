#! /usr/bin/perl -w
use strict;
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use List::Util qw/max min/;
use Cwd;

my $usage=<<USAGE;

Describe: Only used to check unmaped read ...

  Author: Yun Wang, wangyun\@genomic.org.cn
 Version: 1.0, Date: 2012-11-11

Usage:perl $0 <infile1> <infile2> -prefix <str> -o <output_path> 
          
      -prefix     <str>  : Sample ID/prefix of output file     
      -outdir            : Directory of output [default: ./]   
      -help              : show this message
        
      -bwachc            : run to get unmaped reads CHECK information by bwa 
         -tfphred  [int] : Quality transformed to Phred score
         -lowqual  [int] : maximal number of low quality base allowed 
         -Nnum     [int] : maximal number of N base allowed
         -mismatch [int] : maximal number of N mismatch allowed
         
      -bowtiechc         : run to get unmaped reads CHECK information by Bowtie 
   
Example:
   1 get unmaped reads CHECK information by bwa 
     perl $0 <infile1> <infile2> -prefix <str> -bwachc  -lowqual [int] -Nnum [int] -mismatch [int] -o <output_path>
   
   2 get unmaped reads CHECK information by bwa
     perl $0 <infile1> <infile2> <outfile> -prefix <str> -bowtiechc 
         
USAGE

my ($bwachc,$tfphred,$lowquality,$Nnum,$mismatch,$prefix,$bowtiechc,$outdir,$help);
GetOptions(
   "prefix:s" => \$prefix,
   "bwachc" => \$bwachc,
   "tfphred:s" => \$tfphred,
   "lowqual:s"  => \$lowquality,
   "Ns:s"  => \$Nnum,
   "mismatch:s"  => \$mismatch,
   "bowtiechc" => \$bowtiechc,
   "outdir:s"  => \$outdir,
   "help"  => \$help,
);

die $usage if ($help);
$tfphred ||= 64;
$lowquality ||= 10;
$Nnum ||= 4;
$mismatch ||= 2;
$prefix ||= "umapcheck";
$outdir ||= getcwd;

if ($bwachc){
  die $usage if (@ARGV != 2);
  open FQ, "$ARGV[0]" || die $!;
  open FAI, "$ARGV[1]" || die $!;
  open OUT, ">$outdir/$prefix.umc" || die $!;

  my %loxpread;
  $/ = "\n@";
  my $Fline = <FQ>;
  my @line= split /\n/,$Fline;
  $line[0] =~m /^@(.+)/;
  $line[0] = $1;
  #print "$line[0]\n";
  my $chunk1 = join "\n",@line;
  $loxpread{$line[0]} = $chunk1;
  my $Tnum =1;

  while (<FQ>){
	  chomp;
	  my $chunk = $_;
	  my $readid = (split /\n/,$_)[0];
	  $loxpread{$readid} = $chunk;
	  $Tnum++;
  }

  print "Total number of Loxp reads:$Tnum\n";
  close FQ;

  my $prinM ="";
  my $prinUM = "";	
  $/ = "\n";
  while(<FAI>){
	  chomp;
	  next if ($_ =~ /^\@SQ/);
	  my $readsID = "";  
	  my @line = split /\s+/,$_;
	  foreach my $key (keys %loxpread){
	    if ($key =~ /$line[0]/){
		    $readsID = $key;
		    delete $loxpread{$key};
		    last;
	    }
	  }
	  next unless ($readsID); 
	  my $umtype ="";
	  my $mismatch_num = 0;
	  my $Ns =0;
	  my $lowqual = 0;   # ord($Q) - 64 < 7
	  my $q20 = 0;       #phred quality >=20
	  if ($line[1]!=4){
		  my $strand;
		  if ($line[1] == 16) {
		  	$strand = "-";
		  }elsif($line[1] == 0){
			  $strand = "+";
		  }else{$strand = "*";}
	  
	    if ($line[5]=~ /[ID]/){
	  	  $umtype = "Indel";
#=head
      }elsif($line[5]=~ /s/i){
    	  $Ns ++ while ($line[9] =~ /N/gi);
		    my @qual = split //,$line[10];
		    foreach (@qual){
			    my $phred = ord($_) - $tfphred;
			    if ($phred < 7){
			 	    $lowqual ++;
			   	}
		      if ($phred >= 20){
		   	    $q20 ++;
	        }
	      }
	      if($Ns>= $Nnum){
	  	    $umtype = "Ns";
	      }elsif($lowqual >= $lowquality){
	  	    $umtype = "Low_quality";
	      }else{$umtype = "*";} 
	  
	      $prinUM .= sprintf ("%-50s%-14s%-8s%-8s%-8s","$readsID","$umtype","$Ns","$lowqual","$q20");
	      $prinUM .= "\n"; 
	      next;
#=cut
	    }else{
	   	  $mismatch_num  ++ while ($line[18]=~/[ATCGN]/gi);
	   	  $umtype = "Mismatch" if $mismatch_num > $mismatch;
	  	  }
	    if ($umtype eq ""){$umtype ="*";}  		 
        if(length($line[5])>=10){$line[5] = "$line[5]"." ";}	  
	    $prinM .= sprintf ("%-50s%-10s%-4s%-16s%-8s%-10s%-15s","$readsID","$umtype","$strand","$line[2]","$line[3]","$line[5]","$line[18]");
	    $prinM .= "\n";
	  }else{
		  $Ns ++ while ($line[9] =~ /N/gi);
		  my @qual = split //,$line[10];
		  foreach (@qual){
			   my $phred = ord($_) - $tfphred;
			   if ($phred < 7){
			 	    $lowqual ++;
			 	  }
		     if ($phred >= 20){
		   	    $q20 ++;
	       }
	     }
	    if($Ns>=$Nnum){
	  	  $umtype = "Ns";
	    }elsif($lowqual >= $lowquality){
	  	  $umtype = "Low_quality";
	    }else{$umtype = "*";} 
	  
	      $prinUM .= sprintf ("%-50s%-14s%-8s%-8s%-8s","$readsID","$umtype","$Ns","$lowqual","$q20");
	      $prinUM .= "\n"; 
	   }
  }

  print OUT "$prinM$prinUM";

  close OUT; 
}

if ($bowtiechc){

  die $usage if (@ARGV != 3);	
  open IN1, "$ARGV[0]" ||die $!; #JSnum.umc
  open IN2, "$ARGV[1]" ||die $!; #JSnum.splitted.map 
  open OUT, ">$ARGV[2]" ||die $!; #JSnum.split.umc


  my %splitmap;
  while(<IN2>){
	  next if ($_ =~ /^\@/);
    chomp;
    my @line = split /\s+/,$_;	
    my ($readsid,$strand,$chrid,$site,$macht,$md);
    if ($line[1]!=4){		
		  if ($line[1] == 16) {
			  $strand = "-";
		  }elsif($line[1] == 0){
			  $strand = "+";
		  }
	  
	    $readsid = $line[0];
	    $chrid =  $line[2];
	    $site = $line[3];
	    $macht =$line[5];
	    $md = $line[17];
#	    print "$readsid\t$chrid\t$strand\t$site\t$macht\t$md\n"; 
	    if(length($macht)>=10){$macht = "$macht"." ";}	  
	    $splitmap{$readsid} = sprintf ("%-50s%-4s%-16s%-8s%-10s%-15s","$readsid","$strand","$chrid",,"$site","$macht","$md");	
=head
	   }else{
		  $strand = "*";
		 ($readsid,$chrid,$site,$macht,$md) = ("$line[0]","*","*","*","*");
=cut
	   }   
  }

  while(<IN1>) {
	  chomp;
	
	  my @line = split /\s+/,$_;
	
	  if ($line[1] eq "*"){
		  my $spltid1 = "$line[0]"."_1";
		  my $spltid2 = "$line[0]"."_2";
		  if(exists $splitmap{$spltid1}){
		    print OUT "$splitmap{$spltid1}\n";
		    delete $splitmap{$spltid1};
		  }
		  if(exists $splitmap{$spltid2}){
			  print OUT "$splitmap{$spltid2}\n";
			  delete $splitmap{$spltid2};
		  }
	  }else{next;}
  }

  close IN1;
  close IN2;
  close OUT;
}
##===end===##
