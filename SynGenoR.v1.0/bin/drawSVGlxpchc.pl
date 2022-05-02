#! /usr/bin/perl -w

use strict;
use FindBin qw ($Bin);
use lib "$Bin/../lib";
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use List::Util qw/max min/;
use Statistics::Basic qw(:all);
use Cwd;

my $usage=<<USAGE;
Describe: This script is used to draw a loxp reads mapping svg figure for 
          validation of the reconstructed reference. .
          In addition, add split reads mapping negtive check ...
          
  Author: Yun Wang, wangyun\@genomic.org.cn
 Version: 1.0, Date: 2012-11-08

Usage :perl $0 [option]
        -prefix      <str> :  prefix of output file/sample name  
        -insize      [int] :  insert size [default:500]
        -DNAtype     <str> :  DNA type,Linear|Circular [[default:Circular]         
        -refid       <str> :  id of synthetic reference
        -refcode     <str> :  segment index order of synthetic sequence
        -refseq      <str> :  synthetic reference sequence
        -reflen      [int] :  length of synthetic sequence
        -loxpftr     <str> :  loxp feature information of synthetic sequence (gff format)
        -lxpfq       <str> :  reads containing Loxp site sequence (fq format)
        -lxpspr      <str> :  pair_end loxp reads info from normal and abnormal SOAP results (mapping on synthetic reference)
        -lxpspid     <str> :  loxp reads id in SOAP mapping result (mapping on synthetic reference)
        -lxpspstat   <str> :  statistic of loxp reads at each loxp site in synthetic reference
        -rencode     <str> :  segment index order of reconstructed sequence 
        -reseq       <str> :  reconstructed sequence file (scrambled genome)
        -relxpftr    <str> :  loxp feature information of reconstructed sequence (gff format)
        -relxpspr    <str> :  pair_end loxp reads info from normal and abnormal SOAP results (mapping on reconstructed sequence)
        -relxpspid   <str> :  loxp reads id in SOAP mapping result (mapping on reconstructed reference)
        -relxpspstat <str> :  statistic of loxp reads at each loxp site in reconstructed reference
        -reumc       <str> :  ummapping loxp reads CHECK information
        -splitumc    <str> :  amphibolous loxp reads CHECK information
        -outdir            :  output directory [default: ./]
        -help              :  show this help message

USAGE


my ($samplename,$insertsize,$DNAtype,$refid,$refencode,$refseq,$reflen,$loxpftr,$lxpfq,$loxp_spr,$loxp_spid,$lxpsp_stat,$reseq_loxpftr,$rencode,$scbseq,$reseq_loxp_spr,$reseq_loxp_spid,$relxpsp_stat,$reseq_umc,$split_umc,$outdir,$help);

GetOptions(
    "prefix:s"     => \$samplename,
    "insize:s"     => \$insertsize,
    "DNAtype:s"    => \$DNAtype,
    "refid:s"      => \$refid, 
    "refcode:s"    => \$refencode,
    "refseq:s"     => \$refseq,
    "reflen:s"     => \$reflen,
    "loxpftr:s"    => \$loxpftr,
    "lxpfq:s"      => \$lxpfq,
    "lxpspr:s"     => \$loxp_spr,
    "lxpspid:s"    => \$loxp_spid,
    "lxpspstat:s"  => \$lxpsp_stat,
    "rencode:s"    => \$rencode,
    "reseq:s"      => \$scbseq, 
    "relxpftr:s"   => \$reseq_loxpftr, 
    "relxpspr:s"   => \$reseq_loxp_spr,
    "relxpspid:s"  => \$reseq_loxp_spid,
    "relxpspstat:s"=> \$relxpsp_stat,
    "reumc:s"      => \$reseq_umc,
    "splitumc:s"   => \$split_umc,
    "outdir:s"     => \$outdir,
    "help"       => \$help,
);

die $usage if (!$refid||!$refencode||!$refseq||!$loxpftr||!$lxpfq||!$loxp_spr||!$loxp_spid||!$lxpsp_stat);
die $usage if (!$reseq_loxpftr||!$rencode||!$scbseq||!$reseq_loxp_spr||!$reseq_loxp_spid||!$relxpsp_stat||!$reseq_umc||!$split_umc||$help);

print "#$samplename\n";
print "#drawing SVG\n";

$reflen ||= 100371;
$refid ||= "Refseq";
$insertsize ||= 500;
my $version = "v1.0";
my $software = basename ($0); 



open RENC,$refencode ||die $!;
open FSEQ,$refseq ||die $!;
open LXF,$lxpfq ||die $!;
open SPR,$loxp_spr ||die $!;
open SPI,$loxp_spid ||die $!;
open LXS,$lxpsp_stat ||die $!; 

#open RFTR,$reseq_loxpftr ||die $!;
open ENC,$rencode ||die $!;
open RSEQ,$scbseq ||die $!;
open RSPR,$reseq_loxp_spr ||die $!;
open SSPI,$reseq_loxp_spid ||die $!;
open RUMC,$reseq_umc ||die $!;
open SUMC,$split_umc ||die $!;
open RLXS,$relxpsp_stat ||die $!; 

$outdir ||= getcwd;
open OUT, ">$outdir/$samplename.lxpchc.svg" ||die $!;
open XLS, ">$outdir/$samplename.lxpchc.xls" ||die $!;

my $height = estimateY($reseq_umc,$split_umc,410);
$height = $height + 100;

###loading scrambled sequence information###
<RSEQ>;
my $reseq = <RSEQ>;
chomp $reseq;
my $relen = length($reseq);
print "The scrambled reference length: $relen\n";

my $width;

#$reseqlen < $reflen?$width = 10240:$width = int($reseqlen/10)+400;
if ($relen < $reflen){
	$width = estimateX($reflen,$loxpftr);
}else{
	$width = estimateX($relen,$reseq_loxpftr);
}

print OUT "<svg height=\"$height\" width=\"$width\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";

Drabg($width,$height);

###loading loxp fq information####
my %loxpread;
$/ = "\n@";
my $Fline = <LXF>;
my @line= split /\n/,$Fline;
$line[0] =~m /^@(.+)/;
$line[0] = $1;
my $chunk1 = join "\n",@line;
$loxpread{$line[0]} = $chunk1;
my $Tnum =1;

while (<LXF>){
	chomp;
	my $chunk = $_;
	my $readid = (split /\n/,$_)[0];
	$loxpread{$readid} = $chunk;
	$Tnum++;
}

print "Total number of Loxp reads:$Tnum\n";
close LXF;
$/ = "\n";

###loading encode of reference information###
my $refcodeseq = <RENC>;
chomp $refcodeseq;
$refcodeseq = (split /\s+/,$refcodeseq)[1]; 
my @temparray = (split /,/,$refcodeseq);
my $lastcode = pop @temparray;

###loading reference sequence information###
<FSEQ>;
my $seq_ref = <FSEQ>;
chomp $seq_ref;
print "The reference length: $reflen\n";

###start drawing mapping svg of reference loxp reads ####
my $posy = 100;
my $ordinfo = Drawaxis($posy,$reflen,$loxpftr,$refcodeseq);

my $axipos = $posy  - 8;
my $scripty = $axipos - 70; 
Drawscript(5,$scripty,"Reads with loxp site that can regularly map to reference before Scrambling");
$scripty = $axipos - 4;
Drawscript(16,$scripty ,$refid);
$scripty = $posy  + 35;
Drawscript(18,$scripty ,"Segment_ID");


while(<SPR>){
  chomp;
  my ($site1,$ori1,$site2,$ori2) = (split /\s+/,$_)[2,3,6,7];
  next if abs($site2-$site1) > ($insertsize*1.3);
  my $ry = 50; 
  Drawsprds($site1,$ori1,$site2,$ori2,$axipos,$ry,\@{$ordinfo});
}
while (<LXS>){
	chomp;
	next if $_ =~ /^#/;
	my ($lxpstar,$lxpend,$number)= (split /\s+/,$_)[1,2,3];
  $scripty = $posy  + 27;
  Drawlxprnum($lxpstar,$lxpend,$number,\@{$ordinfo},$scripty);
 }

my $mapnum_ref = 0;
my %mapid_ref;
while(<SPI>){
	chomp;
	$mapnum_ref ++;	
	$mapid_ref{$_} = 1;
 }
print "Total number of reference mapped reads: $mapnum_ref\n";
 
###############################################################
###loading encode of scrambled sequence information####
my $codeseq = <ENC>;
chomp $codeseq;
$codeseq = (split /\s+/,$codeseq)[1];
$codeseq =~ s/,0//g;

if($DNAtype eq "Circular"){
  $codeseq =~ s/(-*1),(-*$lastcode)|(-*$lastcode),(-*1)/$1 $2/g;
  #$codeseq =~ s/$lastcode,1/$lastcode 1/g;
 }
 
###loading new construction seq mapped read ID ########
my $mapnum_scb = 0;
my %mapid_scb;
while(<SSPI>){
	chomp;
	$mapnum_scb ++;	
	$mapid_scb{$_} = 1;
 }
print "Total number of scrambled sequence mapped reads: $mapnum_scb\n";

###start drawing loxp reads mapping svg of scrambled sequence ####
$posy = $posy +150;
my $ordinfo2 = Drawaxis($posy,$relen,$reseq_loxpftr,$codeseq);

$axipos = $posy  - 8;

$scripty = $axipos - 70; 
Drawscript(5,$scripty,"Reads with loxp site that can regularly map to reconstructed sequence");
$scripty = $axipos - 4;
Drawscript(16,$scripty ,"Reconstructed_Seq");

$scripty = $posy  + 35;
Drawscript(16,$scripty ,"Segment_ID");

while (<RLXS>){
	chomp;
	next if $_ =~ /^#/;
	my ($lxpstar,$lxpend,$number)= (split /\s+/,$_)[1,2,3];
  $scripty = $posy  + 27;
  Drawlxprnum($lxpstar,$lxpend,$number,\@{$ordinfo2},$scripty);
 }
 
while(<RSPR>){	
  chomp;
  my ($id1,$site1,$ori1,$id2,$site2,$ori2) = (split /\s+/,$_)[0,2,3,4,6,7];
  next if abs($site2-$site1) > 650;
  my $ry = 50;
  if(exists $mapid_ref{$id1}||exists $mapid_ref{$id2}){ 
    Drawsprds($site1,$ori1,$site2,$ori2,$axipos,$ry,\@{$ordinfo2});
  }else{
  	my $rgb = "165,42,42";
  	Drawsprds($site1,$ori1,$site2,$ori2,$axipos,$ry,\@{$ordinfo2},$rgb);
  	}
}


my $tmp_posy = $posy + 70;
###Output the details of ummapped reads(loxp unmap check)#####
$posy += 160; 
Drawscript(5,$posy,"## Unmapping loxp reads check");

my $i = 1;
my $noise = 0;
my @dubious;
while (<RUMC>){
	 chomp; 
	 my @line = split /\s+/;
	 $noise ++ if $line[1] ne "*";
   if ($line[1] eq "*"){
   	push @dubious,[@line];
   	next;}
   $posy += 14;
	 if (@line == 7){
	 my ($id,$umtype,$strand,$chr,$site,$match,$MD) = @line;
	 $chr =~ s/IXR_BACseq_scb/Reconstructed_IXR/;
   Drawscript(8,$posy,$i);
   Drawscript(32,$posy,$id);
   Drawscript(420,$posy,"BWA");
   Drawscript(500,$posy,$umtype);
   Drawscript(580,$posy,$strand);
   Drawscript(600,$posy,$chr);
   Drawscript(750,$posy,$site);
   if(length($match) >=16){
     	$match = substr($match,0,13);
     	$match = "$match"."...";
   	}
   Drawscript(810,$posy,$match);   
   Drawscript(950,$posy,$MD);      
   }else{
   my ($id,$umtype,$Ns,$lowqual,$q20) = @line;
   if ($umtype eq "Ns"){
   	$umtype = "N_base";
   }
   Drawscript(8,$posy,$i);
   Drawscript(32,$posy,$id);
   Drawscript(420,$posy,"-");
   Drawscript(500,$posy,$umtype);
   Drawscript(600,$posy,$Ns);
   Drawscript(670,$posy,$lowqual);
   Drawscript(750,$posy,$q20);
   }
   $i ++;
 }

my @sumc;
if (-s $split_umc){
#  $posy += 30;
#  Drawscript(5,$posy,"## Final split reads check");
#  $i = 1;
   
   my ($id,$strand,$chr,$site,$match,$MD,$mptype);
   
   my $filine = <SUMC>;
   chomp $filine;	 
	 my @filine = split /\s+/,$filine;
	 my ($idtmp,$strandtmp,$chrtmp,$sitetmp,$matchtmp,$MDtmp) = @filine;
   my $readidtmp = $idtmp;
   $readidtmp =~ s/_\d+//;
   my $control = 0;
   
  while(<SUMC>){
	 chomp;
	 $control ++;
	 #$posy += 14; 
	 my @line = split /\s+/;
	 	 
	 ($id,$strand,$chr,$site,$match,$MD) = @line;	 
	 my $readid = $id;
	 $readid =~ s/_\d+//;
	 $mptype = ""; 
	 
	 if($readid ne $readidtmp){
	  $mptype = "SE";
	  }else{$mptype = "PE";}
   push @sumc,[$idtmp,"Bowtie",$mptype,$strandtmp,$chrtmp,$sitetmp,$matchtmp,$MDtmp];
   if ($mptype eq "PE"){
    push @sumc,[$id,"Bowtie",$mptype,$strand,$chr,$site,$match,$MD];
    my $mline = <SUMC>;
    next if defined $readid;
    chomp $mline;
    $control ++;
    my @line = split /\s+/,$mline;
    ($id,$strand,$chr,$site,$match,$MD) = @line;	
     $readid = $id;
	   $readid =~ s/_\d+//;
    ($idtmp,$strandtmp,$chrtmp,$sitetmp,$matchtmp,$MDtmp) = ($id,$strand,$chr,$site,$match,$MD);
    }	   
   ($idtmp,$strandtmp,$chrtmp,$sitetmp,$matchtmp,$MDtmp) = ($id,$strand,$chr,$site,$match,$MD);
   $readidtmp = $readid;
  }
  # $posy += 14; 
   if($control == 0){
    $mptype = "SE";
    push @sumc,[$idtmp,"Bowtie",$mptype,$strandtmp,$chrtmp,$sitetmp,$matchtmp,$MDtmp];
   }elsif($mptype ne "PE"){
   push @sumc,[$idtmp,"Bowtie",$mptype,$strandtmp,$chrtmp,$sitetmp,$matchtmp,$MDtmp];  	
  }
}

my $m = @dubious - 1;
foreach my $p(0..$m){
	 $posy += 14;
	 Drawscript(8,$posy,$i);
   Drawscript(32,$posy,$dubious[$p][0]);
   Drawscript(420,$posy,"-");
   Drawscript(500,$posy,$dubious[$p][1]);
   Drawscript(600,$posy,$dubious[$p][2]);
   Drawscript(670,$posy,$dubious[$p][3]);
   Drawscript(750,$posy,$dubious[$p][4]);
	 my $n = @sumc -1;
	 my $ck = 0;
	 foreach my $q(0..$n){
	   my $readid = $sumc[$q][0]; 
	   $readid =~ s/_\d+//;
	   if ($dubious[$p][0] eq  $readid){
	   	  $ck ++;   	  
	   	  last if $ck > 2;
	   	  $posy += 14;
	   	  #Drawscript(8,$posy,$i);
        Drawscript(32,$posy,$sumc[$q][0]);
        Drawscript(420,$posy,$sumc[$q][1]);
        Drawscript(500,$posy,$sumc[$q][2]);
        Drawscript(580,$posy,$sumc[$q][3]);
        Drawscript(600,$posy,$sumc[$q][4]);
        Drawscript(750,$posy,$sumc[$q][5]);
          if(length($sumc[$q][6]) >=16){
     	      $sumc[$q][6] = substr($sumc[$q][6],0,13);
     	      $sumc[$q][6] = "$sumc[$q][6]"."...";
      	}
        Drawscript(810,$posy,$sumc[$q][6]);   
        Drawscript(950,$posy,$sumc[$q][7]); 
	   	  $sumc[$q][0] = "XX";
	   	}
	  }
	 
	 $i++;
	}
	
##Output summary of loxp reads check##
  my $Tumap_ref = $Tnum - $mapnum_ref;
  my $Tumap_scb = $Tnum - $mapnum_scb;
  my $tfreads =  $Tumap_ref - $Tumap_scb ;
  my $scbreads = $Tumap_ref - $noise;
  
  my $un_tfreads = $scbreads - $tfreads;
  my $Rate_tf =($tfreads/$scbreads)*100;
  my $Rem_Rate_tf  = 100 - $Rate_tf;
  $Rate_tf = sprintf("%.2f%%",$Rate_tf);
  $Rem_Rate_tf = sprintf("%.2f%%",$Rem_Rate_tf);

$posy = $tmp_posy; 
Drawscript(5,$posy,"##Summary of loxp reads mapping");
$posy += 16; 
Drawscript(20,$posy,"[1] $mapnum_ref of $Tnum reads with loxp site can regularly map back to reference before scrambling.");
$posy += 16; 
Drawscript(20,$posy,"[2] For the rest $Tumap_ref reads, $scbreads reads can NOT map back to reference before scrambling because of recombination event.");
$posy += 16;
Drawscript(50,$posy,"After sequence reconstruction,  $Rate_tf of the $scbreads reads can regularly map back to new reconstructed reference and the rest $Rem_Rate_tf and the remaining $noise of $mapnum_ref reads still can NOT map back to the new reconstructed reference.");
$posy += 16; 
Drawscript(50,$posy,"For all these unmapping reads, addtional analysis is performed with following results:");

print XLS "$samplename\t$Tnum\t$mapnum_ref\t$Tumap_ref\t$mapnum_scb\t$Tumap_scb\t$tfreads\t$noise\t$scbreads\t$Rate_tf\t$un_tfreads\n";  

$scripty = $height - 20;
my $loctime = scalar (localtime); 
Drawscript(25,$scripty,"## Created by $software $version at $loctime.");

print OUT "</svg>";


####======sub function=========##### 

### draw description #####
sub Drawscript{
	 my ($posx,$posy,$script,$size) = @_;
	 $size ||= 14;
	 print OUT "<text fill=\"rgb(47,79,79)\" font=\"Arial\" font-size=\"$size\" font-weight=\"bold\" x=\"$posx\" y=\"$posy\">$script </text>\n";    
}

##### draw background#####
sub Drabg{
	my ($width,$height)= @_;
	my $stary = 0;
	my $endy = $height;
  my $starx = 0;
  my $endx = $width;
  my $len = $endx - $starx + 1;
 
	my $n = int ($len/100); 

	foreach my $i(0..$n){
		 my $starxm = $i*100;
		 print OUT "<line style=\"fill: rgb(211,211,211); fill-opacity: 1.0; stroke: rgb(211,211,211); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxm\" x2=\"$starxm\" y1=\"$stary\" y2=\"$endy\" />\n" ;

	       last if ($i== $n);
	       foreach my $j(1..9){
	       	  my $starxs = (($i*1000)+ $j*100)/10 ;
		        print OUT "<line style=\"fill: rgb(224,255,255); fill-opacity: 1.0; stroke: rgb(224,255,255); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxs\" x2=\"$starxs\" y1=\"$stary\" y2=\"$endy\" />\n";
	         }
	   } 
	my $m = int (($len%100)/10) + 1;   
	foreach my $k (1..$m){
	   my $starxs = (($n*1000)+ $k*100)/10;
		 print OUT "<line style=\"fill: rgb(224,255,255); fill-opacity: 1.0; stroke: rgb(224,255,255); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxs\" x2=\"$starxs\" y1=\"$stary\" y2=\"$endy\" />\n";
	 }       	
}
	
#### Draw  separator #####
sub separate{
	my ($x,$y)=@_;
	my $sprtup = $y - 6;
	my $sprtdn = $y + 6;	
	my $sprtl = $x - 4;
	my $sprtr = $x + 4;	
	print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$sprtr\" x2=\"$sprtl\" y1=\"$sprtup\" y2=\"$sprtdn\" />\n";
  $sprtl += 10;
  $sprtr += 10;
  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$sprtr\" x2=\"$sprtl\" y1=\"$sprtup\" y2=\"$sprtdn\" />\n";
	}
	
##### Draw code seq #####
sub Drawcode{
	my($pos_code,$loxpos1,$loxpos2,$code) = @_;
	my $posup = $pos_code - 4;
	my $posdn = $pos_code + 4;
	my $marklen = length($code);
	my $chrwidth = $marklen*8 + 4;
	my $xlen = (($loxpos2 - $loxpos1) - $chrwidth)/2;
	my $x1 = $loxpos1 + $xlen;
  my $x2 = $loxpos2 - $xlen;
     if($marklen>5){$x2=$x2 - 10;}
  my $txtX = $x1 + 4;
  my $txtY = $pos_code + 4;
  
  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$loxpos1\" x2=\"$loxpos1\" y1=\"$posup\" y2=\"$posdn\" />\n"; 
	print OUT "<text fill=\"rgb(0,0,0)\" font=\"Helvetica\" font-size=\"11\" font-weight=\"normal\" x=\"$txtX\" y=\"$txtY\">$code </text>\n";
  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$loxpos2\" x2=\"$loxpos2\" y1=\"$posup\" y2=\"$posdn\" />\n";
  
  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$loxpos1\" x2=\"$x1\" y1=\"$pos_code\" y2=\"$pos_code\" />\n";
  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$x2\" x2=\"$loxpos2\" y1=\"$pos_code\" y2=\"$pos_code\" />\n"; 
 }
    	
##### Draw coordinate axis #####
sub Drawaxis{
	my($pos,$seqlen,$loxpftr,$codeseq) = @_; 
	my $movex = 0;                                    #modify the X axis coordination..
	my $dist = 650;
	open FTR,$loxpftr ||die $!;
	#my ($pos,$seqlen,$height,$oper)= @_;
	my $starym = $pos - 4;
	my $endym = $pos + 4;
	my $starys = $pos - 2;
	my $endys = $pos + 2;
	my $sprtup = $pos - 6;
	my $sprtdn = $pos + 6;
	my $codin = 0;
	my $mx = int($seqlen/10) + 120;
	my $mxa = $mx - 4;
	
	my @tford;
	my @lxpinfo;
  my $i = 0;
  my $lastmean = 0;
  my $meantt;
	my @code= split /,/,$codeseq;
	my $fftr = <FTR>;
	my ($tloxps,$tloxpe) = (split /\s+/,$fftr)[3,4];
	push @lxpinfo,[$tloxps,$tloxpe];
	#justfy if the length of the first edge >= 750;
	
	my $tmean = mean($tloxps,$tloxpe);
	if ($tmean >=$dist + 100){
		my $starsite = $tmean - $dist;
		$movex = 120;
		my $x_start =  $movex;
		my $x_end = int($dist/10) +  $movex;
	  push @{$tford[$i]},($starsite,$tmean,$x_start,$x_end);
	 
	  $i++;
	  $movex = $x_end;	  
	  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"80\" x2=\"110\" y1=\"$pos\" y2=\"$pos\" />\n";
	  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"80\" x2=\"84\" y1=\"$pos\" y2=\"$starym\" />\n";
	  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"80\" x2=\"84\" y1=\"$pos\" y2=\"$endym\" />\n";
	  
	  ###0k coordinate###
	  my $TextX = 95;
	  my $TextY = $pos + 16;
	  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"100\" x2=\"100\" y1=\"$starym\" y2=\"$endym\" />\n"; 
	  print OUT "<text fill=\"rgb(0,0,0)\" font=\"Helvetica\" font-size=\"11\" font-weight=\"normal\" x=\"$TextX\" y=\"$TextY\">0k </text>\n";
	  separate(110,$pos);
#	  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"114\" x2=\"106\" y1=\"$sprtup \" y2=\"$sprtdn\" />\n";
#	  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"124\" x2=\"116\" y1=\"$sprtup\" y2=\"$sprtdn\" />\n";	
	}else{
		my $starsite = 0;
		$movex = 100;
		my $x_start = int($starsite/10) +  $movex;
		my $x_end = int($tmean/10) +  $movex;
	  push @{$tford[$i]},($starsite,$tmean,$x_start,$x_end);
	  
	  $i++;
	  $movex = $x_end;
	  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"80\" x2=\"110\" y1=\"$pos\" y2=\"$pos\" />\n";
	  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"80\" x2=\"84\" y1=\"$pos\" y2=\"$starym\" />\n";
	  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"80\" x2=\"84\" y1=\"$pos\" y2=\"$endym\" />\n";
		}
	
	while (<FTR>){
	 my ($loxps,$loxpe) = (split /\s+/,$_)[3,4];
	 push @lxpinfo,[$loxps,$loxpe];
	 my $mean = mean($loxps,$loxpe);
	 $lastmean = $mean;
	 #my $meantt;
	 #justfy if the length of the first edge >= 1400;
   if ($mean - $tmean >= $dist*2 + 200){
      my $site1 = $tmean + $dist;
		  my $x1 = $movex;
	 	  my $x2 = int($dist/10) +  $movex;
	    push @{$tford[$i]},($tmean,$site1,$x1,$x2);
	    
	    $i++;
	    $movex = $x2 + 10;
	    
	    separate($x2,$pos);

	    my $site2 = $mean - $dist;
	    $x1 = $movex;
	    $x2 = int($dist/10) +  $movex;
	    push @{$tford[$i]},($site2,$mean,$x1,$x2);
      $i++;
	   
	    $movex = $x2;
	    $tmean = $mean;	     
	  }else{
	  	my $tmeantt = $mean; 
	  
	  	while(<FTR>){
	  		chomp;
	  		my ($loxpstt,$loxpett) = (split /\s+/,$_)[3,4];
	  		push @lxpinfo,[$loxpstt,$loxpett];
	      $meantt = mean($loxpstt,$loxpett);
	      
	      $lastmean = $meantt;
	      if ($meantt - $tmeantt >= $dist*2 + 200){
	     	  last;
	     	}else{$tmeantt = $meantt;}
	    }

	   if($tmeantt == $meantt||$tmeantt == $lastmean ){           #justify if the last second edge <= $dist*2 + 200
	    next;
	   }else{	   	
	    my $site1 = $tmeantt + $dist;
		  my $x1 = $movex;
	 	  my $x2 = int(($site1-$tmean)/10) +  $movex;
	    push @{$tford[$i]},($tmean,$site1,$x1,$x2);
	    
	    $i++;
	    $movex = $x2 + 10;
	    
	    separate($x2,$pos);

	    my $site2 = $meantt - $dist;
	    $x1 = $movex;
	    $x2 = int($dist/10) +  $movex;
	    push @{$tford[$i]},($site2,$meantt,$x1,$x2);
	  	    
	    $i++;
	    $movex = $x2;
	    $tmean = $meantt;
	    }	 
	  }
  }

##lastregion....
  
  if ($seqlen - $lastmean > $dist + 100){
  	my $lastsite = $lastmean + $dist; 
  	my $x_start =  $movex;
		my $x_end = int(($lastsite-$tmean)/10) +  $movex;
	  push @{$tford[$i]},($tmean,$lastsite,$x_start,$x_end);
	 	
	  $i++;
	  $movex = $x_end + 10;	  
	  separate($x_end,$pos);
	  
	  
	  my $lastpos_star = $seqlen - 400;
	  my $lastpos_end = $seqlen + 200;
	  $x_start = $movex;
	  $x_end =  $movex + 60;
    push @{$tford[$i]}, ($lastpos_star,$lastpos_end,$x_start,$x_end);
     
	  $i ++;
	  $movex = $x_end ;	 
	  my $arrowx = $x_end - 4;
	  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$arrowx\" x2=\"$x_end\" y1=\"$starym\" y2=\"$pos\" />\n";
	  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$arrowx\" x2=\"$x_end\" y1=\"$endym\" y2=\"$pos\" />\n";	  
   }else{
    my $x_start =  $movex;
    my $x_end = int(($seqlen -$tmean)/10) + $movex + 20;
    push @{$tford[$i]},($tmean,$seqlen,$x_start,$x_end);
   
    my $arrowx = $x_end - 4;
	  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$arrowx\" x2=\"$x_end\" y1=\"$starym\" y2=\"$pos\" />\n";
	  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$arrowx\" x2=\"$x_end\" y1=\"$endym\" y2=\"$pos\" />\n";	  

    }
    
##merge tford####
  my @mergord;
  my $k = @tford - 1;

  $i = 0;                  # index of @tford
  my $j = 0;               # index of @mergord
  while($i < $k){
  	if($tford[$i][1]==$tford[$i+1][0]){
  	 
  	 push @{$mergord[$j]},($tford[$i][0],$tford[$i+1][1],$tford[$i][2],$tford[$i+1][3]); 	   	 
  	 
  	 print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$mergord[$j][2]\" x2=\"$mergord[$j][3]\" y1=\"$pos\" y2=\"$pos\" />\n";
  	 $j ++;	
  	}else{
  	 print "Warning,you make a mistake in a..!!\n";
  	 
  		}
  	$i +=2;	
  }
#  if($k%2 ==0){
  	push @{$mergord[$j]},($tford[$k][0],$tford[$k][1],$tford[$k][2],$tford[$k][3]);
    print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$mergord[$j][2]\" x2=\"$mergord[$j][3]\" y1=\"$pos\" y2=\"$pos\" />\n";

#   }
   

  my $n = int ($seqlen/1000); 
  my $p = $j;
	foreach my $i(0..$n){
		 my $starxm = $i*1000;		 
		 my $TextX;
		 my $TextY = $pos + 16;
		 my $ik ="$i"."k";
		 foreach my $ind(0..$p){ 
		 	
		 if($starxm >= $mergord[$ind][0] && $starxm <= $mergord[$ind][1]){
		 	 $starxm =  $mergord[$ind][2] + ($starxm -$mergord[$ind][0])/10;
		 	 $TextX = $starxm - 5;
		   print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxm\" x2=\"$starxm\" y1=\"$starym\" y2=\"$endym\" />\n"; 
	     print OUT "<text fill=\"rgb(0,0,0)\" font=\"Helvetica\" font-size=\"11\" font-weight=\"normal\" x=\"$TextX\" y=\"$TextY\">$ik </text>\n";
	     last;
	    }
	  }
	   #print "ind =15 : $mergord[15][0],,$mergord[15][1],,$mergord[15][2],,$mergord[15][3]\n";
	   last if ($i== $n);
	   foreach my $j(1..9){
	     my $starxs = (($i*1000)+ $j*100);
	     foreach my $ind(0..$p){ 
	     	 if($starxs >= $mergord[$ind][0] && $starxs <= $mergord[$ind][1]){   	  
	          $starxs =  $mergord[$ind][2] + ($starxs -$mergord[$ind][0])/10;
	          print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxs\" x2=\"$starxs\" y1=\"$starys\" y2=\"$endys\" />\n"; 
	          last;
	        }
	      }
	    }
 }   
  my $m = int (($seqlen%1000)/100) + 1;
  foreach my $k (1..$m){
	   my $starxs = ($n*1000)+ $k*100;
	   foreach my $ind(0..$p){
	   	 if($starxs >= $mergord[$ind][0] && $starxs <= $mergord[$ind][1]){   	 
	        $starxs =  $mergord[$ind][2] + ($starxs -$mergord[$ind][0])/10;
	        print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxs\" x2=\"$starxs\" y1=\"$starys\" y2=\"$endys\" />\n"; 
	        last;
	     }
	   } 
	 }       	
###add the loxp site add codeseq##
  my $lpn = @lxpinfo;
  my $pos_code = $pos + 33;
  @code = ();
  @code = split /,/,$codeseq;
  my $tloxpos = 100;
  
  for (my $i=0;$i<$lpn;$i++){
    my $loxpos = ($lxpinfo[$i][0]+$lxpinfo[$i][1])/2;
    foreach my $ind(0..$p){
	  if($loxpos >= $mergord[$ind][0] && $loxpos <= $mergord[$ind][1]){   	 
	     $loxpos =  $mergord[$ind][2] + ($loxpos -$mergord[$ind][0])/10;
	     Drawcode($pos_code,$tloxpos,$loxpos,$code[$i]);
	     Drawloxp($loxpos,$pos);
	     last;
    }
   }
   $tloxpos = $loxpos;
 }
 my $endpos = $mergord[-1][3] - 20;
 Drawcode($pos_code,$tloxpos,$endpos,$code[-1]);
 
return (\@mergord);
}	 
##### Draw lxp mapping reads number #####
sub Drawlxprnum{
  my ($lxpstar,$lxpend,$number,$ordinfo,$posy) = @_;
  my $posx = 0;
  my $p = @{$ordinfo} - 1;
	foreach my $ind(0..$p){ 
		
	   if($lxpstar >= $$ordinfo[$ind][0] && $lxpend <= $$ordinfo[$ind][1]){
	   
	   	  my $mid = ($lxpstar+$lxpend)/2;   	  
	      $posx =  $$ordinfo[$ind][2] + ( $mid -$$ordinfo[$ind][0])/10 - 6;
	      print OUT "<text fill=\"rgb(153,50,204)\" font=\"Helvetica\" font-size=\"11\" font-weight=\"bold\" x=\"$posx\" y=\"$posy\">$number </text>\n";

	      
	      last;
	     }
	 }
  
}
    
#### Draw Loxp feature ####
sub Drawloxp {
	 my ($posx,$posy) = @_;	 
	 $posy = $posy - 4;
	 print OUT "<rect style=\"fill: rgb(107,142,35); fill-opacity: 1.0; stroke: rgb(107,142,35); stroke-opacity: 1.0; stroke-width: 1\" height=\"8\" width=\"2.5\" x=\"$posx\" y=\"$posy\" />\n";
 }
 
#### Draw spr reads #####
sub Drawsprds{
	my ($site1,$ori1,$site2,$ori2,$posy,$ry,$ordinfo,$rgb) = @_;
	
	
	my $rgba;
	if($ori1 eq $ori2){
	 	$rgba = "165,42,42";
	}else{$rgba = "106,90,205";}
   
  $rgb ||=	$rgba;

		
	my ($starx,$endx,$rx)=(0,0,0); 
	my $p = @{$ordinfo} - 1;
	foreach my $ind(0..$p){ 
	   if($site1 >= $$ordinfo[$ind][0] && $site1 <= $$ordinfo[$ind][1]){   	  
	      $starx =  $$ordinfo[$ind][2] + ($site1 -$$ordinfo[$ind][0])/10;
	      last;
	     }
	 }
	foreach my $ind(0..$p){ 
	   if($site2 >= $$ordinfo[$ind][0] && $site2 <= $$ordinfo[$ind][1]){   	  
	      $endx =  $$ordinfo[$ind][2] + ($site2 -$$ordinfo[$ind][0])/10;
	      last;
	     }
	 }        
  my $max = max($starx,$endx);
  my $min = min($starx,$endx);
  $starx =  $min;
  $endx = $max;
	$rx = ($endx-$starx)/2;
	print OUT "<path d=\"M$starx,$posy,A$rx,$ry,0,0,1,$endx,$posy\" style=\"fill: none; fill-opacity: 1.0; stroke: rgb($rgb); stroke-opacity: 1.0; stroke-width: 0.8; stroke-linecap: square\"/>\n";
}

sub estimateY{
	my ($reseq_umc,$split_umc,$height) = @_;	
	open UMC, "$reseq_umc" ||die $!;
	open SPC, "$split_umc" ||die $!;
	while(<UMC>){
	 	$height += 14;
	  }
	if(-s  $split_umc){ 
   while(<SPC>){
	 	$height += 14;
	  }
	}	  
	  return $height;
}
sub estimateX{
	my ($seqlen,$loxpftr) = @_;
	my $lastmean = 0;
  my $movex = 0;                                    #modify the X axis coordination..
	my $dist = 650;
	open FTR,$loxpftr ||die $!;
		
	my $fftr = <FTR>;
	my ($tloxps,$tloxpe) = (split /\s+/,$fftr)[3,4];
	#justfy if the length of the first edge >= 750;	
	my $tmean = mean($tloxps,$tloxpe);
	if ($tmean >=$dist + 100){
		my $starsite = $tmean - $dist;
		$movex = 120;
		my $x_start =  $movex;
		my $x_end = int($dist/10) +  $movex;
	  $movex = $x_end;	  
	}else{
		my $starsite = 0;
		$movex = 100;
		my $x_start = int($starsite/10) +  $movex;
		my $x_end = int($tmean/10) +  $movex; 
	  $movex = $x_end;
	}
	
	while (<FTR>){
	 my ($loxps,$loxpe) = (split /\s+/,$_)[3,4];
	 my $mean = mean($loxps,$loxpe);
	 $lastmean = $mean;
	 #justfy if the length of the first edge >= 1400;
   if ($mean - $tmean >= $dist*2 + 200){
      my $site1 = $tmean + $dist;
		  my $x1 = $movex;
	 	  my $x2 = int($dist/10) +  $movex;
	    $movex = $x2 + 10;
	    
	    my $site2 = $mean - $dist;
	    $x1 = $movex;
	    $x2 = int($dist/10) +  $movex;
	    $movex = $x2;
	    $tmean = $mean;	     
	  }else{
	  	my $tmeantt = $mean; 
	  	my $meantt;
	  	while(<FTR>){
	  		chomp;
	  		my ($loxpstt,$loxpett) = (split /\s+/,$_)[3,4];
	      $meantt = mean($loxpstt,$loxpett);
	      $lastmean = $meantt;
	      if ($meantt - $tmeantt >= $dist*2 + 200){
	     	  last;
	     	}else{$tmeantt = $meantt;}
	    }
	   if($tmeantt == $meantt){           #justify if the last second edge <= $dist*2 + 200
	    next;
	   }else{
	    my $site1 = $tmeantt + $dist;
		  my $x1 = $movex;
	 	  my $x2 = int(($site1-$tmean)/10) +  $movex;
	    $movex = $x2 + 10;
	    my $site2 = $meantt - $dist;
	    $x1 = $movex;
	    $x2 = int($dist/10) +  $movex;  
	    $movex = $x2;
	    $tmean = $meantt;
	    }	 
	  }
  }
##lastregion....
   
  if ($seqlen - $lastmean > $dist + 100){
  	my $lastsite = $lastmean + $dist; 
  	my $x_start =  $movex;
		my $x_end = int(($lastsite-$tmean)/10) +  $movex;
	  $movex = $x_end + 10;	  
	  
	  my $lastpos_star = $seqlen - 400;
	  my $lastpos_end = $seqlen + 200;
	  $x_start = $movex;
	  $x_end =  $movex + 60;
	  $movex = $x_end ;	 
	}else{
    my $x_start =  $movex;
    my $x_end = int(($seqlen -$tmean)/10) + $movex + 20;
    $movex = $x_end ;	
    }
  my $width = $movex +100;
  return $width;   
}
