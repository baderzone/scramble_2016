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
Describe: dynamic programming path(Dynpath) to reconstruct sequence 
          and get index order...
                              
  Author: Yun Wang, wangyun\@genomic.org.cn
 Version: 1.0, Date: 2012-11-05

Usage :perl $0 <edg_file> <node_file> [option]
     -spid   <str> : sample name that shown in the prefix of output file [default: path_result] 
     -mcf    <int> : maximal cumulate confict in one path allowed
     -outdir       : output directory [default: ./]
     -help         : show this message
     
Note : Only Dynamic programming based  method is available. Euler path is in development. 

 
USAGE

my ($samplename,$max_conflict,$outdir,$help);
GetOptions(  
   "spid:s" => \$samplename,
   "mcf:s" => \$max_conflict,
   "outdir:s"  => \$outdir,
   "help"  => \$help,
);

die $usage if (@ARGV != 2 ||$help);

my ($ire_file,$vis_file)= @ARGV;

$samplename ||= "path_result";
$outdir ||= getcwd;
unless($max_conflict==0){$max_conflict ||= 3;}
open IRE,"$ire_file" ||die $!;
open VIS,"$vis_file" ||die $!;
open CODE,">$outdir/$samplename.encod" ||die $!;
open INFO,">$outdir/$samplename.info" ||die $!;
open LOGG,">$outdir/$samplename.log" ||die $!;

##loading edge and vertice info
my (%edgecopy,%regind);
my %vertice;
while (<IRE>){
	print LOGG "$_";
	next if ($_=~ /^#/);
  chomp;
  my @line = (split /\s+/,$_);
  my $len = @line;
  $edgecopy{$line[0]} = $line[1]; 
  $regind{$line[0]} = $line[4];
  #print LOGG "##$line[0]\t$line[1]\t$line[2]\n";    
}

my $edge_number = keys %edgecopy;
print LOGG "Total number of edge:$edge_number\n";

while (<VIS>){
	print LOGG "$_";
  next if ($_=~m/^#/);
  chomp;
##change##  
  my @line = (split /\s+/);
  #$vertice{$line[0]} = $line[1];
  $vertice{$line[0]} = 1;  
}
my $vertice_number = keys %vertice;
print LOGG "Total number of vertice:$vertice_number\n";

close IRE;
close VIS;

##Path by dynamic programming
my %path;
my @edgecopykeys = sortkey(\%edgecopy,1);

#=head
foreach (@edgecopykeys){
  if($edgecopy{$_}>0){
    push @{$path{1}},($_,0); 
    last;
  }
}
#=cut

#push @{$path{1}},("3_2",0);

#print "${$path{1}}[0]\n"; 
my $perfpathnum = 1;
my $illepathnum = 1;

while(1){
 #print  LOGG "start....\n";
 my %tempath =();
 %tempath= %path;
  
 foreach my $ind(keys %tempath){	      
   
   my ($conflictinfo,$cumulate_conflict) = lastpath(${$tempath{$ind}}[0],\%edgecopy,\%vertice);
   my ($indexseq,$control_error);
   my $conflict = ${$tempath{$ind}}[1];
   #print  LOGG "#############$conflict\n"; 
   #Judge the extension permission 
   if ($conflict >$max_conflict){
   	  ($indexseq,$control_error) = indexorder(${$tempath{$ind}}[0],\%regind);
   	  print INFO "##Find illegal path$illepathnum:\nconflict =$conflict\ncumulate conflict=$cumulate_conflict\nPath= ${$tempath{$ind}}[0]\nIndex sequence: $indexseq\n\n" if $cumulate_conflict <=$max_conflict+3;
      delete $path{$ind};
      $illepathnum++ if $cumulate_conflict <=$max_conflict+3;
  	  next;
   }
   
   #print LOGG "${$tempath{$ind}}[0],${$tempath{$ind}}[1],$cumulate_conflict\t";
   #print LOGG "$conflictinfo\n";
   my $ck = 0;
   
   #Judge the perfect Path..
   if ($cumulate_conflict == 0){
      my ($indexseq,$control_error) = indexorder(${$tempath{$ind}}[0],\%regind);
      if($control_error==0){
      	  $ck == 0?print CODE "$samplename\t$indexseq\n":print CODE "$samplename\_$ck\t$indexseq\n";
        }
      print INFO "##Find a perfect path$perfpathnum:\ncumulate conflict=$cumulate_conflict\nPath= ${$tempath{$ind}}[0]\nIndex sequence: $indexseq\n\n";
      $ck++;
      delete $path{$ind};
      $perfpathnum++;
  	  next; 
   }
      
   my %edgenum;
   my %dotnum;
   my @edge = split /\//,${$tempath{$ind}}[0];
   my @dot = split /_/,${$tempath{$ind}}[0];
   my $last_side= pop @dot;
   shift @dot;
   
      
   #choose a possible vertice and edge##add a vertice and a edge.
   my $i= 1;
   my $temp_path = ${$tempath{$ind}}[0]; 
   my $new_ind;
   my $newpath = "";
   $conflict = 0; 
   foreach my $verticekey (keys %vertice){
   	  my ($vertice_L,$vertice_R)=split /\//,$verticekey;
   	  my $dotright = ""; #record the right side of possible vertice
   	  my $edgeadd = "";  #record the possible edge
   	     if($last_side == $vertice_L || $last_side == $vertice_R){   	        
   	         $dotright = $vertice_R if ($last_side == $vertice_L);
   	         $dotright = $vertice_L if ($last_side == $vertice_R); 
   	         
   	         foreach my $edgekey(keys %edgecopy){
   	         	  my ($edgekey_L,$edgekey_R) = (split /_/,$edgekey)[0,1];
   	         	  if ($edgekey_L == $dotright||$edgekey_R == $dotright){
   	         	      $edgeadd= $edgekey if $edgekey_L == $dotright;
   	         	      $edgeadd= "$edgekey_R"."_$edgekey_L" if $edgekey_R == $dotright;	         	  
   	         	   
   	         	      $new_ind = "$ind"."_$i";   	  	             	  	           
   	         	      $newpath = "$temp_path"."/$edgeadd";  	         	      
   	         	      $conflict = conflictpath($newpath,\%edgecopy,\%vertice,1);	         	      
   	         	      push @{$path{$new_ind}},($newpath,$conflict);
   	         	      $i ++;
   	         	      
   	         	   }
   	           }
   	        }   	  	 
 
         }
     my $testind = "$ind"."_1";		 
     unless(exists $path{$testind}){
     ($indexseq,$control_error) = indexorder(${$path{$ind}}[0],\%regind);
	 print INFO "##Find illegal path$illepathnum:\nconflict =$conflict\ncumulate conflict=$cumulate_conflict\nPath= ${$tempath{$ind}}[0]\nIndex sequence: $indexseq\n\n" if $cumulate_conflict <=$max_conflict+3;
	 $illepathnum++ if $cumulate_conflict <=$max_conflict+3;
	 }
     delete $path{$ind} if (exists $path{$ind});    
         
    }
  last if (keys %path == 0); 
 }

     
###   	  
############sub function###############   

##calculate the edge and vertice conflict in the Path.   

#0-1/4-5/7   -->operation 0
#0-1/4-5/7-6 -->operation 1

sub conflictpath{
   	my ($pathway,$edgeN,$vert,$operation) = @_;
   	my $conflict=0;
   	my %edgenum;
    my %dotnum; 
    
    my @edge = split /\//,$pathway;
    pop @edge if $operation==0;
#    my @dot = split /_/,$pathway;
#    pop @dot if $operation==1;
#    shift @dot;
    
    foreach (0..@edge-1){
    	 next if (exists ${$edgeN}{$edge[$_]});
    	 my ($edgeL,$edgeR)= split /_/,$edge[$_];
    	 $edge[$_]="$edgeR"."_$edgeL" if $edgeL>$edgeR;
    }
#    foreach (0..@dot-1){
#    	 next if (exists ${$vert}{$dot[$_]});
#    	 my ($dotL,$dotR)= split /\//,$dot[$_];
#    	 $dot[$_]="$dotR"."/dotL" if $dotL>$dotR;
#    }
    #calculate the edge and vertice in the subpath.
    foreach (@edge){
   	  #(exists $edgenum{$_})?$edgenum{$_} = $edgenum{$_} + 1:$edgenum{$_} = 1;
      if(exists $edgenum{$_}){$edgenum{$_} = $edgenum{$_} + 1;}else{$edgenum{$_} = 1;}
    }
    
##change##    
#    foreach (@dot){
      #(exists $dotnum{$_})?$dotnum{$_} = $dotnum{$_} + 1:$dotnum{$_} = 1;
#      if(exists $dotnum{$_}){$dotnum{$_} = $dotnum{$_} + 1;}else{$dotnum{$_} = 1;} 		
#   }
       
    #calculate the edge conflict.
    foreach (keys %{$edgeN}){
    	my $n = 0;
    	if (exists $edgenum{$_}){
    		$n =  $edgenum{$_}-${$edgeN}{$_};
    		$conflict += $n if $n>0;
     }
   }
##change##   
   #calculate the vertice conflict.
#   foreach(keys %{$vert}){
#   	    my $n = 0;
#    	if (exists $dotnum{$_}){
#    		$n =  $dotnum{$_}-${$vert}{$_};
#    		$conflict += $n if $n>0; 
#       }
#     }
     return ($conflict);
  }
      
   
##calculate the edge and vertice conflict of Last Path.    
sub lastpath{
   	my ($pathway,$edgeN,$vert) = @_;
   	my $len1= keys %{$edgeN};
   	my $len2= keys %{$vert};
   	my $cumulate_conflict=0;
   	my $conflictinfo="";
   	my %edgenum;
    my %dotnum; 
    
    my @edge = split /\//,$pathway;
    my @dot = split /_/,$pathway;
    pop @dot;
    shift @dot;

    foreach (0..@edge-1){
       next if (exists ${$edgeN}{$edge[$_]});          	  
    	 my ($edgeL,$edgeR)= split /_/,$edge[$_];
    	 my $temp = $edge[$_];
    	 $edge[$_]="$edgeR"."_$edgeL" if $edgeL>$edgeR;
 
    }
   
    foreach (0..@dot-1){
    	 next if (exists ${$vert}{$dot[$_]});
    	 my ($dotL,$dotR)= split /\//,$dot[$_];
    	 $dot[$_]="$dotR"."/$dotL" if $dotL>$dotR;    
    }
    
    #calculate the edge and vertice in the subpath.
    foreach (@edge){
   	  #(exists $edgenum{$_})?$edgenum{$_} = $edgenum{$_} + 1:$edgenum{$_} = 1;
   	  if(exists $edgenum{$_}){$edgenum{$_} = $edgenum{$_} + 1;}else{$edgenum{$_} = 1;}
    }
    foreach (@dot){
      #(exists $dotnum{$_})?$dotnum{$_} = $dotnum{$_} + 1:$dotnum{$_} = 1; 		
      if(exists $dotnum{$_}){$dotnum{$_} = $dotnum{$_} + 1;}else{$dotnum{$_} = 1;}
    }
       
    #calculate the edge conflict.
    foreach (keys %{$edgeN}){
    	my $n = 0;
    	if (exists $edgenum{$_}){
    		$n =  $edgenum{$_}-${$edgeN}{$_};
    		$conflictinfo .= "$_:"."$n ";
    		$cumulate_conflict += abs($n);
    	}else{
    		$n = ${$edgeN}{$_};
    		$conflictinfo .= "$_:"."-$n ";
    		$cumulate_conflict += $n;    
     }
   }
   
   #calculate the vertice conflict.
   $conflictinfo .= "\t";
   my $n = 0;
   foreach(keys %{$vert}){
   	    
##change##  
#   my $n = 0;  	
#    	if (exists $dotnum{$_}){
#    		$n =  $dotnum{$_}-${$vert}{$_};
#    		$conflictinfo .= "$_:"."$n ";
#    		$cumulate_conflict += abs($n);
#    	}else{
#    		$n = ${$vert}{$_};
#    		$conflictinfo .= "$_:"."-$n ";
#    		$cumulate_conflict += $n;    
#      }
      $n++ unless (exists $dotnum{$_});     
     }
     $cumulate_conflict += $n; 
     return ($conflictinfo,$cumulate_conflict);
  }
  
##get the index order of the scrambled sequence
sub indexorder{
	  my ($pathway,$regind) = @_;
	  my @edge = split /\//,$pathway;
	  my $indexseq = "";
	  my $control_error=0;
	  my $i=0;  
	  foreach (@edge){
	  	 $i ++;
	     my($edgeL,$edgeR) = split /_/,$_;	     
	     my $reverse_edge = "$edgeR"."_$edgeL";
	     	if(exists ${$regind}{$_}){
	     	    #$i==1?$indexseq = ${$regind}{$_}:$indexseq .= ","."${$regind}{$_}";
	     	    if ($i==1){
	     	    	$indexseq = ${$regind}{$_};
	     	    }else{$indexseq .= ","."${$regind}{$_}";}
	     	    	
	     	         	     	     	   
	     	}elsif(exists ${$regind}{$reverse_edge}){
	     		 my @temparray = split /,/,${$regind}{$reverse_edge};
	     		 @temparray = reverse @temparray;
	     		 my $idex = @temparray -1;	     		 
	     		 foreach my $i(0..$idex){
	     		    $temparray[$i] = "-$temparray[$i]";
	     		    }
	     		 my $reverse_regind = join ",",@temparray;	 
	     		 #$i==1?$indexseq = $reverse_regind:$indexseq .= ","."$reverse_regind"; 
	     	   if ($i==1){
	     	    	$indexseq = $reverse_regind;
	     	   }else{$indexseq .= ","."$reverse_regind";}
	     	    		     		   			     		 
	     	}else{
           print "Warning:find a illegal edge in the path --$_\n";
           $control_error ++; 
         }
        
       }
    return ($indexseq,$control_error);
   }      
            
######sort keys########
sub sortkey{
	 my ($hash,$opern) = @_;
   my @keyid = keys %{$hash};
   my $findex;
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

##
    
  
   
 
