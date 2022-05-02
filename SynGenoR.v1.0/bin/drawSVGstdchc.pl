#! /usr/bin/perl -w

use strict;
use FindBin qw ($Bin);
use lib "$Bin/../lib";
use File::Basename;
use Data::Dumper;
use Getopt::Long;
use List::Util qw/max min/;
use Cwd;

my $usage=<<USAGE;
Describe: This script is used to draw a standard reads mapping svg figure for 
          validation of the reconstructed reference. .
   
  Author: Yun Wang, wangyun\@genomic.org.cn
 Version: 1.0, Date: 2012-11-08

Usage :perl $0 [option]
        -prefix   <str> :  prefix of output file 
        -insize   [int] :  insert size [default:500] 
        -refid    <str> :  id of synthetic reference
        -reflen   [int] :  length of synthetic chromosome
        -reseq    <str> :  reconstructed sequence /scrambled genome
        -loxpftr  <str> :  loxp feature information (gff format)
        -cvg      <str> :  average sequence coverage 
        -lxps     <str> :  mapping stat of loxp reads 
        -gvr      <str> :  abnormal mapped split reads (*.sort.gvr)
        -spr      <str> :  pair_end mapping info from normal and abnormal SOAP results --mapping on synthetic reference
        -varid    <str> :  tansformation Loxp reads id  (from ummapped on synthetic reference to mapped on reconstructed referenece) 
        -val2     <str> :  pair_end mapping info result from abnormal SOAP result (mapping on reconstructed sequence)
        -val      <str> :  pair_end mapping info result from normal SOAP result (mapping on reconstructed sequence)
        -lsv      <str> :  structure variation information  
        -outdir         :  output directory [default: ./]
        -help           :  show this help message


USAGE

my ($samplename,$insertsize,$refid,$reflen,$reseq,$loxpftr,$cvg,$lxs,$lsv,$gvr,$spr,$varid,$val2,$val,$outdir,$help);
GetOptions(
    "prefix:s"     => \$samplename,
    "insize:s"     => \$insertsize,
    "refid:s"      => \$refid,
    "reflen:s"     => \$reflen,
    "reseq:s"      => \$reseq,
    "loxpftr:s"    => \$loxpftr,
    "cvg:s"        => \$cvg,
    "lxps:s"        => \$lxs,    
    "gvr:s"        => \$gvr,
    "spr:s"        => \$spr,
    "varid:s"      => \$varid,
    "val2:s"       => \$val2,
    "val:s"        => \$val,
    "lsv:s"        => \$lsv,
    "outdir:s"     => \$outdir,
    "help"       => \$help,
);


die $usage if (!$refid||!$reseq||!$loxpftr||!$cvg||!$lxs||!$gvr||!$spr||!$varid||!$val2||!$val||!$reflen||$help);

print "#$samplename\n";
print "#drawing SVG\n";

$reflen ||= 100371;
$refid ||= "Refseq";
$insertsize ||= 500;
my $version = "v1.0";
my $software = basename ($0); 

open RSEQ, $reseq ||die $!;
open FTR, $loxpftr ||die $!;
open CVG, $cvg ||die $!;
open LXS, $lxs ||die $!;
open SPR, $spr ||die $!;
open VID, $varid ||die $!;
open VAL2, $val2 ||die $!;
open VAL, $val ||die $!;

my $height = estimateY($gvr,$lsv,80);
#print "###$height\n";
$height = $height + 500;

<RSEQ>;
my $re_seq = <RSEQ>;
chomp $re_seq;
my $reseqlen= length($re_seq);
my $width;

#$reseqlen < $reflen?$width = 10240:$width = int($reseqlen/10)+400;
if ($reseqlen < $reflen){$width = int($reflen/10) + 200;
}else{$width = int($reseqlen/10)+200;}

open GVR, $gvr ||die $!;
$outdir ||= getcwd;
open OUT, ">$outdir/$samplename.stdchc.svg" ||die $!;
print OUT "<svg height=\"$height\" width=\"$width\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">\n";

Drawscript(5,35,"Coverage");    
while (<CVG>){
	chomp;
	my ($site,$depth) = (split /\s+/,$_)[1,3];
  Drawdpline(50,$site,$depth);
}
my $axipos = 58;
Drawscript(5,54,"$refid");  

Drawaxis ($axipos,$reflen,$height,1);
=head
if($reseqlen < $reflen){
  Drawaxis ($axipos,$reflen,$height,1);
}else{
	Drawaxis ($axipos,$reflen,$height,0);
	}
=cut

my $loxpY = 80;
Drawscript(5,80,"Loxp site");

while (<FTR>){
	 next if $_ =~ /^#/;
	 my ($star,$END) = (split /\s+/,$_)[3,4];
	 
	 Drawloxp($star,$loxpY);
	}

my %lsv = ();
my $posysv = 100;

if($lsv && -s $lsv){
   open LSV, $lsv ||die $!;
  while (<LSV>){
	  next if $_ =~ /^#/;
	  my ($star,$END,$tstar,$tend) = (0,0,0,0);
	  my ($type,$direc);
	  $type = (split /\s+/,$_)[0];
  	if ($type eq "DEL" || $type eq "INV"){
		   ($type,$star,$END) = (split /\s+/,$_)[0,1,2];
	    	push @{$lsv{$type}},[$star,$END];	
  	}else {
		   ($type,$star,$END,$tstar,$tend,$direc) = (split /\s+/,$_)[0,1,2,3,4,5];
		   push @{$lsv{$type}},[$star,$END,$tstar,$tend,$direc];	
   }
  }



  if (exists $lsv{"DEL"}){
 	  Drawscript(5,100,"Deletion");
 	  my @temp = [0,0];
 	  my $i = 0;
 	  foreach my $del(@{$lsv{"DEL"}}){
 	  	  foreach my $tempa(@temp){
 	  	  	if(${$tempa}[0]>${$del}[1] ||${$tempa}[1]<${$del}[0]){$i = 0;
 	  	  	}else{$i = 1;last;}
 	  	  }
 	  	   if ($i == 1){$posysv += 20;}	
 	  	   $i = 0; 	  		    
 	  	   Drawvariate ($posysv,"DEL",${$del}[0],${$del}[1]);
 	  	   push @temp,[${$del}[0],${$del}[1]];
 	  }
 	 }

  if (exists $lsv{"INV"}){
	  $posysv += 20;
    Drawscript(5,$posysv,"Inversion");
 	  my @temp = [0,0];
 	  my $i = 0;
 	  foreach my $inv(@{$lsv{"INV"}}){
 	  	  foreach my $tempa(@temp){
 	  	  	  if(${$tempa}[0]>${$inv}[1] ||${$tempa}[1]<${$inv}[0]){$i = 0;
 	  	  	  }else{$i = 1;last;}
 	  	  	}
 	  	 if ($i == 1){$posysv += 20;}	  	  		    
 	  	 $i = 0;
 	  	 Drawvariate ($posysv,"INV",${$inv}[0],${$inv}[1]);
 	  	 push @temp,[${$inv}[0],${$inv}[1]];
 	 
 	  	}
 	 }
	  	 
  if (exists $lsv{"TDUP"}){
	  $posysv += 20;
    Drawscript(5,$posysv,"Tandem duplication");
 	  my @temp = [0,0];
 	  my $i = 0;
 	  foreach my $tx(@{$lsv{"TDUP"}}){
 	  	  my ($txmax,$txmin);
 	  	  foreach my $tempa(@temp){
 	  	  	  $txmax = max (${$tx}[0],${$tx}[1],${$tx}[2],${$tx}[3]);
 	  	  	  $txmin = min (${$tx}[0],${$tx}[1],${$tx}[2],${$tx}[3]);
 	  	  	  if(${$tempa}[0]>$txmax ||${$tempa}[1]<$txmin){$i = 0;
 	  	  	  }else{$i = 1;last;}
 	  	  	}
 	  	   if ($i == 1){$posysv += 20;}	  	  		    
 	  	   $i = 0;
 	  	   Drawvariate ($posysv,"TDUP",${$tx}[0],${$tx}[1],${$tx}[2],${$tx}[3],${$tx}[4]);
 	  	   push @temp,[$txmin,$txmax];
 	  	}
 	  }
  if (exists $lsv{"NTDUP"}){
	  $posysv += 20;
    Drawscript(5,$posysv,"Nontandem duplication");
 	  my @temp = [0,0];
 	  my $i = 0;
 	  foreach my $tx(@{$lsv{"NTDUP"}}){
 	  	  my ($txmax,$txmin);
 	  	  foreach my $tempa(@temp){
 	  	  	  $txmax = max (${$tx}[0],${$tx}[1],${$tx}[2],${$tx}[3]);
 	  	  	  $txmin = min (${$tx}[0],${$tx}[1],${$tx}[2],${$tx}[3]);
 	  	      if(${$tempa}[0]>$txmax ||${$tempa}[1]<$txmin){$i = 0;
 	  	  	  }else{$i = 1;last;}
 	  	  	}
 	  	   if ($i == 1){$posysv += 20;}	  	  		    
 	  	   $i = 0;
 	  	   Drawvariate ($posysv,"NTDUP",${$tx}[0],${$tx}[1],${$tx}[2],${$tx}[3],${$tx}[4]);
 	  	   push @temp,[$txmin,$txmax];
 	  	}
 	  }
  if (exists $lsv{"D&I"}){
	  $posysv += 20;
    Drawscript(5,$posysv,"Deletion and insertion");
 	  my @temp = [0,0];
 	  my $i = 0;
 	  foreach my $tx(@{$lsv{"D&I"}}){
 	  	  my ($txmax,$txmin);
 	  	  foreach my $tempa(@temp){
 	  	  	  $txmax = max (${$tx}[0],${$tx}[1],${$tx}[2],${$tx}[3]);
 	  	  	  $txmin = min (${$tx}[0],${$tx}[1],${$tx}[2],${$tx}[3]);
 	  	      if(${$tempa}[0]>$txmax ||${$tempa}[1]<$txmin){$i = 0;
 	  	  	  }else{$i = 1;last;}
 	  	  	}
 	  	   if ($i == 1){$posysv += 20;}	  	  		    
 	  	   $i = 0;
 	  	   Drawvariate ($posysv,"D&I",${$tx}[0],${$tx}[1],${$tx}[2],${$tx}[3],${$tx}[4]);
 	  	   push @temp,[$txmin,$txmax];
 	  	}
 	 }
}	   	   	  
 	  	
my $posylxp = $posysv + 10;
my $scripty = $posylxp + 12;
Drawscript(5,$scripty,"Loxp mapping reads");
	
while (<LXS>){
	 next if $_ =~ /^#/;  
	 chomp;
	 my ($site,$fowrds,$revrds) = (split /\s+/,$_)[1,4,5];
	 Drawloxprds ($site,$posylxp,$fowrds,$revrds);
	}

$scripty = $posysv + 50;
Drawscript(5,$scripty,"Split reads");
my @recpos;
my $lingvr = <GVR>;
chomp $lingvr;
my ($direct,$lstar,$lend,$lori,$rstar,$rend,$rori)= (split /\s+/,$lingvr)[0,3,4,5,7,8,9];
my $posy = $posylxp + 30;
my $suposy = 0;
Drawsplitrds($direct,$lstar,$lend,$lori,$rstar,$rend,$rori,$posy);
my $max = max($lend,$rstar);
push @recpos, [$max,$posy];	

while (<GVR>){
	$posy += 16;
	chomp;
	($direct,$lstar,$lend,$lori,$rstar,$rend,$rori)= (split /\s+/,$_)[0,3,4,5,7,8,9];
  my $minpos = min ($lend,$rstar);
  my $maxpos = max ($lend,$rstar);
  if($minpos <$max){
      Drawsplitrds($direct,$lstar,$lend,$lori,$rstar,$rend,$rori,$posy);
      push @recpos,[$maxpos,$posy];
      $max = $maxpos;	
  }else{
  	 my @recposR = reverse @recpos;
  	 foreach my $record (@recposR){
  	 	 ${$record}[0] < $minpos?$posy= ${$record}[1]:last; 
  	 	 }
  	 @recpos = ();	 
  	 Drawsplitrds($direct,$lstar,$lend,$lori,$rstar,$rend,$rori,$posy);
  	 push @recpos,[$maxpos,$posy];
  	 $max = $maxpos;
   }
  $suposy = $posy if $posy > $suposy; }
print "The maximal y-axis coordinate of splitted reads:$suposy\n";

#$axipos = ($height-450)+80;
$axipos = $suposy + 100;
Drawaxis ($axipos,$reflen,$height,0);

open FTR, $loxpftr ||die $!;
$loxpY = $axipos - 4; 
while (<FTR>){
	 next if $_ =~ /^#/;
	 my ($star,$END) = (split /\s+/,$_)[3,4];
	 
	 Drawloxp($star,$loxpY);
	}


$posy = $axipos - 8;
$scripty = $axipos - 70;
Drawscript(5,$scripty,"Reads irregular mapping to reference before scrambling");
$scripty = $axipos - 30;
#Drawscript(5,$scripty,"Irregular mapping");
$scripty = $axipos - 4;
#Drawscript(5,$scripty,"UnSCRaMbLEd");

while(<SPR>){
  chomp;
  my ($site1,$ori1,$site2,$ori2) = (split /\s+/,$_)[2,3,6,7];
  my $ry = 50;
  Drawsprds($site1,$ori1,$site2,$ori2,$posy,$ry);
}

####################################
my %vrid;
my $Tnum=0;
while(<VID>){
	chomp;
	next if $_=~ /^#/;
	my ($id,$ctgy) = (split /\s+/)[0,1];
	$vrid{$id} = $ctgy;
	$Tnum++;
}

close VID;

####################################
$axipos = $axipos+100;

Drawaxis ($axipos,$reseqlen,$height,0);
=head
if ($reflen <= $reseqlen){
	Drawaxis ($axipos,$reseqlen,$height,1);
 }else{
	Drawaxis ($axipos,$reseqlen,$height,0);
}
=cut

$posy = $axipos - 8;
$scripty = $axipos - 70;
Drawscript(5,$scripty,"Validation by regular(up) and irregular (down) reads mapping to new reference");
#$scripty = $axipos - 40;
#Drawscript(5,$scripty,"Regular mapping");
$scripty = $axipos - 4;
#Drawscript(5,$scripty,"SCRaMbLEd");

my $tfN=0;
while(<VAL>){
  chomp;
  my ($id,$site1,$ori1,$site2,$ori2) = (split /\s+/,$_)[0,2,3,6,7];
  my $ry = 50;
  my $rgb;
  if ($id =~ /(.+)\/\d+$/){$id = $1;}
  
  if(exists $vrid{$id}){
  	$tfN++;
  	if($vrid{$id}==0){
  	  $rgb=	"34,139,34";
  	}elsif($vrid{$id}==1){
  		$rgb = "165,42,42";
  	}elsif($vrid{$id}==2){
  	  $rgb = "106,90,205";
  	}  
  Drawsprds($site1,$ori1,$site2,$ori2,$posy,$ry,$rgb);
 }else{
 	Drawsprds($site1,$ori1,$site2,$ori2,$posy,$ry);
 	}
}

close VAL;
####################################
$axipos = $axipos+80;
Drawaxis ($axipos,$reseqlen,$height,0);

$posy = $axipos - 8;
$scripty = $axipos - 40;
#Drawscript(5,$scripty,"Irregular mapping");
$scripty = $axipos - 4;
#Drawscript(5,$scripty,"SCRaMbLEd");

while(<VAL2>){
  chomp;
  my ($site1,$ori1,$site2,$ori2) = (split /\s+/,$_)[2,3,6,7];
  my $ry = 50; 
  Drawsprds($site1,$ori1,$site2,$ori2,$posy,$ry) if (abs($site2-$site1)>$insertsize||$ori1 eq $ori2);
}
close VAL2;

my $tfrate = int(($tfN/$Tnum)*1000)/10;

$scripty = $axipos + 50;
Drawscript(25,$scripty,"In total, $Tnum reads can not map back to reference before scrambling.");
#$scripty = $scripty + 25;
#Drawscript(25,$scripty,"Total $tfN of reads change into regular mapping");
$scripty = $scripty + 25;
Drawscript(25,$scripty,"After construction of actual reference sequence, $tfrate% of these irregular reads can map back to new constructed reference.");
$scripty = $scripty + 25;
if ($tfrate>90){
   #Drawscript(25,$scripty,"Perfect validation!!");
}

$scripty = $height - 20;
my $loctime = scalar (localtime); 
Drawscript(25,$scripty,"## Created by $software $version at $loctime.");
print OUT "</svg>";

###======sub function=========##### 

### draw description #####
sub Drawscript{
	 my ($posx,$posy,$script) = @_;
	 print OUT "<text fill=\"rgb(47,79,79)\" font=\"Arial\" font-size=\"14\" font-weight=\"bold\" x=\"$posx\" y=\"$posy\">$script </text>\n";    
}
 
### draw sequence coverage #### 
sub Drawdpline{
	my ($starY,$starX,$depth) = @_;
	$starX = ($starX/10)+100;
	my $endY =$starY*(1-$depth);  
	print OUT "<line style=\"fill: rgb(0,0,255); fill-opacity: 1.0; stroke: rgb(0,0,255); stroke-opacity: 1.0; stroke-width: 0.5; stroke-linecap: square\" x1=\"$starX\" x2=\"$starX\" y1=\"$starY\" y2=\"$endY\" />\n";
}

##### draw coordinate axis and background#####
sub Drawaxis{
	my ($pos,$seqlen,$height,$oper)= @_;
	my $starym = $pos - 4;
	my $endym = $pos + 4;
	my $starys = $pos - 2;
	my $endys = $pos + 2;
	my $codin = 0;
	my $mx = int($seqlen/10) + 120;
	my $mxa = $mx - 4;
	
	print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"80\" x2=\"$mx\" y1=\"$pos\" y2=\"$pos\" />\n";
	print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"80\" x2=\"84\" y1=\"$pos\" y2=\"$starym\" />\n";
	print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"80\" x2=\"84\" y1=\"$pos\" y2=\"$endym\" />\n";
	print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$mx\" x2=\"$mxa\" y1=\"$pos\" y2=\"$starym\" />\n";
	print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$mx\" x2=\"$mxa\" y1=\"$pos\" y2=\"$endym\" />\n";
			
 
	my $n = int ($seqlen/1000); 

	foreach my $i(0..$n){
		 my $starxm = $i*100+100;
		 my $TextX = $starxm-5;
		 my $TextY = $pos +16;
		 my $ik ="$i"."k";
		   
		 print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxm\" x2=\"$starxm\" y1=\"$starym\" y2=\"$endym\" />\n"; 
		 print OUT "<line style=\"fill: rgb(211,211,211); fill-opacity: 1.0; stroke: rgb(211,211,211); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxm\" x2=\"$starxm\" y1=\"$pos\" y2=\"$height\" />\n" if $oper == 1;
	   print OUT "<text fill=\"rgb(0,0,0)\" font=\"Helvetica\" font-size=\"11\" font-weight=\"normal\" x=\"$TextX\" y=\"$TextY\">$ik </text>\n";
	       last if ($i== $n);
	       foreach my $j(1..9){
	       	  my $starxs = (($i*1000)+ $j*100)/10 + 100;
	       	  print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxs\" x2=\"$starxs\" y1=\"$starys\" y2=\"$endys\" />\n"; 
		        print OUT "<line style=\"fill: rgb(224,255,255); fill-opacity: 1.0; stroke: rgb(224,255,255); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxs\" x2=\"$starxs\" y1=\"$pos\" y2=\"$height\" />\n"if $oper == 1;
	         }
	   }
  my $m = int (($seqlen%1000)/100) + 1;
  foreach my $k (1..$m){
	   my $starxs = (($n*1000)+ $k*100)/10 + 100;
	   print OUT "<line style=\"fill: rgb(0,0,0); fill-opacity: 1.0; stroke: rgb(0,0,0); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxs\" x2=\"$starxs\" y1=\"$starys\" y2=\"$endys\" />\n"; 
		 print OUT "<line style=\"fill: rgb(224,255,255); fill-opacity: 1.0; stroke: rgb(224,255,255); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starxs\" x2=\"$starxs\" y1=\"$pos\" y2=\"$height\" />\n"if $oper == 1;
	 }       	
}	

#### Draw Loxp feature ####
sub Drawloxp {
	 my ($site,$posy) = @_;
	 my $posx = $site/10+100;
	 print OUT "<rect style=\"fill: rgb(107,142,35); fill-opacity: 1.0; stroke: rgb(107,142,35); stroke-opacity: 1.0; stroke-width: 1\" height=\"8\" width=\"2.5\" x=\"$posx\" y=\"$posy\" />\n";
 }

### Draw structure variation ###  
sub Drawvariate {
	 my @para = @_;
	 my ($star,$END,$tstar,$tend) = (0,0,0,0);
	 my ($type,$direc); 
	 if (@para == 4){
	 	  ($posy,$type,$star,$END) = @para;
	 	  my $posx = $star/10+100;
	 	  my $width = ($END - $star + 1)/10;	 	  
	 	  if ($type eq "DEL"){	 	  	
	 	  	print OUT "<rect style=\"fill: rgb(105,105,105); fill-opacity: 1.0; stroke: rgb(105,105,105); stroke-opacity: 1.0; stroke-width: 1\" height=\"8\" width=\"$width\" x=\"$posx\" y=\"$posy\" />\n";
	 	  }elsif($type eq "INV"){
	 	  	print OUT "<rect style=\"fill: rgb(165,42,42); fill-opacity: 1.0; stroke: rgb(165,42,42); stroke-opacity: 1.0; stroke-width: 1\" height=\"8\" width=\"$width\" x=\"$posx\" y=\"$posy\" />\n";
      }else{print "lsv format error\n";}
   }else{
   	  ($posy,$type,$star,$END,$tstar,$tend,$direc) = @para;
   	  
   	  my $posx = $star/10+100;
	 	  my $width = ($END - $star + 1)/10;
	 	  my ($pathsx,$pathex) = (0,0);
	 	  
      if($END <= $tstar){
	 	   	 $pathsx = $END/10+100;
	 	     $pathex = $tstar/10+100;
	 	  }else{	 	     
	 	   	 $pathsx = ($star/10)+100;
	 	     $pathex = ($tstar/10)+100;
	 	   	}
	 	   	
	 	  
	 	  my $script = "$direc";
	 	  my $arrowl = $pathex-3;
	 	  my $arrowr = $pathex+3;
	 	  my $posym = $posy + 4;
	 	  my $posymd = $posy + 2;
	 	  my $posyd = $posy + 8;
	 	  my $posyu = $posy - 2;
	 	  my $rgb;
	 	  
	 	  if($type eq "TDUP"){
	 	   $rgb = "25,25,112";#midnightblue
	 	   
	 	  }elsif($type eq "NTDUP"){
	 	   $rgb = "255,165,0";# orange
	 	  }elsif($type eq "D&I"){
	 	   $rgb = "148,0,211";#purple
	 			   		 	      	
      }else{
      print "Weird structure variationn !!\n"	
     }
      print OUT "<rect style=\"fill: rgb($rgb); fill-opacity: 1.0; stroke: rgb($rgb); stroke-opacity: 1.0; stroke-width: 1\" height=\"8\" width=\"$width\" x=\"$posx\" y=\"$posy\" />\n";
      print OUT "<line style=\"stroke-dasharray: 2; fill: none; fill-opacity: 1.0; stroke: rgb(119,136,153); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$pathsx\" x2=\"$pathex\" y1=\"$posym\" y2=\"$posym\" />\n";
      print OUT "<line style=\"fill: rgb(25,25,112); fill-opacity: 1.0; stroke: rgb(119,136,153); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$pathex\" x2=\"$pathex\" y1=\"$posyu\" y2=\"$posyd\" />\n";
     	print OUT "<line style=\"fill: rgb(25,25,112); fill-opacity: 1.0; stroke: rgb(119,136,153); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$pathex\" x2=\"$arrowl\" y1=\"$posyd\" y2=\"$posymd\" />\n";
      print OUT "<line style=\"fill: rgb(25,25,112); fill-opacity: 1.0; stroke: rgb(119,136,153); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$pathex\" x2=\"$arrowr\" y1=\"$posyd\" y2=\"$posymd\" />\n";
      print OUT "<text fill=\"rgb(119,136,153)\" font=\"Arial\" font-size=\"15\" font-weight=\"bold\" x=\"$arrowl\" y=\"$posyu\">$script </text>\n";    
     
  }
}

##### Draw loxp Reads ######
 
sub Drawloxprds{
	  my ($site,$posy,$fowrds,$revrds) = @_;
	  my $posx = $site/10+100;
	  
	  if($fowrds > 0){
	  my $posxs = $posx - 6;
	  my $posxe = $posx + 8;
	  my $posxa = $posxe - 3;
	  my $posyu = $posy + 1;
	  my $posyl = $posy + 4;
		my $posyd = $posy + 7;
	  my $posxt = $posx + 9;
	  my $posyt = $posy + 8;
	  
	  print OUT "<line style=\"fill: rgb(0,139,139); fill-opacity: 1.0; stroke: rgb(0,139,139); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$posxs\" x2=\"$posxe\" y1=\"$posyl\" y2=\"$posyl\" />\n";
	  print OUT "<line style=\"fill: rgb(0,139,139); fill-opacity: 1.0; stroke: rgb(0,139,139); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$posxe\" x2=\"$posxa\" y1=\"$posyl\" y2=\"$posyu\" />\n";
    print OUT "<line style=\"fill: rgb(0,139,139); fill-opacity: 1.0; stroke: rgb(0,139,139); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$posxe\" x2=\"$posxa\" y1=\"$posyl\" y2=\"$posyd\" />\n";
    print OUT "<rect style=\"fill: rgb(107,142,35); fill-opacity: 1.0; stroke: rgb(107,142,35); stroke-opacity: 1.0; stroke-width: 1\" height=\"8\" width=\"2.5\" x=\"$posx\" y=\"$posy\" />\n";
    print OUT "<text fill=\"rgb(0,0,0)\" font=\"Helvetica\" font-size=\"10\" font-weight=\"normal\" x=\"$posxt\" y=\"$posyt\">$fowrds </text>\n";    
  }
###    
    if ($revrds>0){
    $posy = $posy + 12;	
    my $posxs = $posx - 8;
    my $posxa = $posxs + 3;
	  my $posxe = $posx + 6;
	  my $posyu = $posy + 1;
	  my $posyl = $posy + 4;
	  my $posyd = $posy + 7;
	  my $posxt = $posx + 9;
	  my $posyt = $posy + 8;
    print OUT "<line style=\"fill: rgb(0,139,139); fill-opacity: 1.0; stroke: rgb(0,139,139); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$posxs\" x2=\"$posxe\" y1=\"$posyl\" y2=\"$posyl\" />\n";
	  print OUT "<line style=\"fill: rgb(0,139,139); fill-opacity: 1.0; stroke: rgb(0,139,139); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$posxs\" x2=\"$posxa\" y1=\"$posyl\" y2=\"$posyu\" />\n";
    print OUT "<line style=\"fill: rgb(0,139,139); fill-opacity: 1.0; stroke: rgb(0,139,139); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$posxs\" x2=\"$posxa\" y1=\"$posyl\" y2=\"$posyd\" />\n";
    print OUT "<rect style=\"fill: rgb(107,142,35); fill-opacity: 1.0; stroke: rgb(107,142,35); stroke-opacity: 1.0; stroke-width: 1\" height=\"8\" width=\"2.5\" x=\"$posx\" y=\"$posy\" />\n";
    print OUT "<text fill=\"rgb(0,0,0)\" font=\"Helvetica\" font-size=\"10\" font-weight=\"normal\" x=\"$posxt\" y=\"$posyt\">$revrds </text>\n";    
   }
}

##### Draw splitted loxp reads ###
sub Drawrds {
	my ($star,$loxpend,$ori,$posy) = @_;
  my ($starx,$endx,$arxo,$arx,$aruy,$ardy); 
  
  $aruy = $posy - 3;
  $ardy = $posy + 3;	
	if ($star>$loxpend){
		  $starx = ($loxpend/10) + 100;
		  $endx = ($star/10 +100) + 10;
	}else{       
		  $starx = ($star/10+100) - 10;
		  $endx = ($loxpend/10) + 100;
	}
	
	if ($ori eq "+"){
		  $arxo = $endx;
		  $arx = $arxo - 3;
		}else{
		  $arxo = $starx;
	    $arx = $arxo + 3;
	 }

#Only uesd for regulate picture!!
	 if ($star>$loxpend && $ori eq "-"){
	 	   $arxo = $arxo + 3;
	 	   $arx = $arxo + 3;
	 	   }	 	  
##END
	 
	print OUT "<line style=\"fill: rgb(75,0,130); fill-opacity: 1.0; stroke: rgb(75,0,130); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$starx\" x2=\"$endx\" y1=\"$posy\" y2=\"$posy\" />\n";
  print OUT "<line style=\"fill: rgb(75,0,130); fill-opacity: 1.0; stroke: rgb(75,0,130); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$arxo\" x2=\"$arx\" y1=\"$posy\" y2=\"$aruy\" />\n";
  print OUT "<line style=\"fill: rgb(75,0,130); fill-opacity: 1.0; stroke: rgb(75,0,130); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$arxo\" x2=\"$arx\" y1=\"$posy\" y2=\"$ardy\" />\n";
}

	  
sub Drawsplitrds{
	 my ($direct,$lstar,$lend,$lori,$rstar,$rend,$rori,$posy) = @_;
	 my $rgb  = "0,0,0";
	 if($direct eq "R"){
	 	  $rgb = "165,42,42";
	 	}else {$rgb = "106,90,205";}
	 
	 my $loxpl = $lend/10 +100;
	 my $loxpr = $rstar/10 +100; 
	 
	 print OUT "<rect style=\"fill: rgb($rgb); fill-opacity: 1.0; stroke: rgb($rgb); stroke-opacity: 1.0; stroke-width: 1\" height=\"8\" width=\"1.5\" x=\"$loxpl\" y=\"$posy\" />\n";
	 print OUT "<rect style=\"fill: rgb($rgb); fill-opacity: 1.0; stroke: rgb($rgb); stroke-opacity: 1.0; stroke-width: 1\" height=\"8\" width=\"1.5\" x=\"$loxpr\" y=\"$posy\" />\n";
	 $posy = $posy + 4;
	 if ($direct eq "F"){
	 print OUT "<line style=\"stroke-dasharray: 2; fill: none; fill-opacity: 1.0; stroke: rgb(119,136,153); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$loxpl\" x2=\"$loxpr\" y1=\"$posy\" y2=\"$posy\" />\n";
	}else{
	 print OUT "<line style=\"stroke-dasharray: 2; fill: none; fill-opacity: 1.0; stroke: rgb(255,192,203); stroke-opacity: 1.0; stroke-width: 1; stroke-linecap: square\" x1=\"$loxpl\" x2=\"$loxpr\" y1=\"$posy\" y2=\"$posy\" />\n";	
		}
	 Drawrds ($lstar,$lend,$lori,$posy);
	 Drawrds ($rend,$rstar,$rori,$posy);
	} 
	
### estimate the maximal y-axis coordinate ####	     
sub estimateY{
  my ($gvr,$lsv,$posy) = @_;
  #print "$posy\n";  
  if(-s $lsv){
  open LSV, $lsv ||die $!;
  while (<LSV>){
	next if $_ =~ /^#/;
	my ($star,$END,$tstar,$tend) = (0,0,0,0);
	my ($type,$direc);
	$type = (split /\s+/,$_)[0];
	if ($type eq "DEL" || $type eq "INV"){
		 ($type,$star,$END) = (split /\s+/,$_)[0,1,2];
	   	push @{$lsv{$type}},[$star,$END];	
	}else {
		 ($type,$star,$END,$tstar,$tend,$direc) = (split /\s+/,$_)[0,1,2,3,4,5];
		 push @{$lsv{$type}},[$star,$END,$tstar,$tend,$direc];	
    }
  }
 close LSV;
}
 my $posysv = 0;

 if (exists $lsv{"DEL"}){
 	  my @temp = [0,0];
 	  my $i = 0;
 	  foreach my $del(@{$lsv{"DEL"}}){
 	  	  foreach my $tempa(@temp){
 	  	  	if(${$tempa}[0]>${$del}[1] ||${$tempa}[1]<${$del}[0]){$i = 0;
 	  	  	}else{$i = 1;last;}
 	  	  }
 	  	   if ($i == 1){$posysv += 20;}	
 	  	   $i = 0; 	  		    
 	  	   push @temp,[${$del}[0],${$del}[1]];
 	  }
 	 }

if (exists $lsv{"INV"}){
	  $posysv += 20;
 	  my @temp = [0,0];
 	  my $i = 0;
 	  foreach my $inv(@{$lsv{"INV"}}){
 	  	  foreach my $tempa(@temp){
 	  	  	  if(${$tempa}[0]>${$inv}[1] ||${$tempa}[1]<${$inv}[0]){$i = 0;
 	  	  	  }else{$i = 1;last;}
 	  	  	}
 	  	 if ($i == 1){$posysv += 20;}	  	  		    
 	  	 $i = 0;
 	  	 push @temp,[${$inv}[0],${$inv}[1]];
 	 
 	  	}
 	 }
	  	 
if (exists $lsv{"TDUP"}){
	  $posysv += 20;
 	  my @temp = [0,0];
 	  my $i = 0;
 	  foreach my $tx(@{$lsv{"TDUP"}}){
 	  	  my ($txmax,$txmin);
 	  	  foreach my $tempa(@temp){
 	  	  	  $txmax = max (${$tx}[0],${$tx}[1],${$tx}[2],${$tx}[3]);
 	  	  	  $txmin = min (${$tx}[0],${$tx}[1],${$tx}[2],${$tx}[3]);
 	  	  	  if(${$tempa}[0]>$txmax ||${$tempa}[1]<$txmin){$i = 0;
 	  	  	  }else{$i = 1;last;}
 	  	  	}
 	  	   if ($i == 1){$posysv += 20;}	  	  		    
 	  	   $i = 0;
 	  	   push @temp,[$txmin,$txmax];
 	  	}
 	  }
if (exists $lsv{"NTDUP"}){
	  $posysv += 20;
 	  my @temp = [0,0];
 	  my $i = 0;
 	  foreach my $tx(@{$lsv{"NTDUP"}}){
 	  	  my ($txmax,$txmin);
 	  	  foreach my $tempa(@temp){
 	  	  	  $txmax = max (${$tx}[0],${$tx}[1],${$tx}[2],${$tx}[3]);
 	  	  	  $txmin = min (${$tx}[0],${$tx}[1],${$tx}[2],${$tx}[3]);
 	  	      if(${$tempa}[0]>$txmax ||${$tempa}[1]<$txmin){$i = 0;
 	  	  	  }else{$i = 1;last;}
 	  	  	}
 	  	   if ($i == 1){$posysv += 20;}	  	  		    
 	  	   $i = 0;
 	  	   push @temp,[$txmin,$txmax];
 	  	}
 	  }
if (exists $lsv{"D&I"}){
	  $posysv += 20;
 	  my @temp = [0,0];
 	  my $i = 0;
 	  foreach my $tx(@{$lsv{"D&I"}}){
 	  	  my ($txmax,$txmin);
 	  	  foreach my $tempa(@temp){
 	  	  	  $txmax = max (${$tx}[0],${$tx}[1],${$tx}[2],${$tx}[3]);
 	  	  	  $txmin = min (${$tx}[0],${$tx}[1],${$tx}[2],${$tx}[3]);
 	  	      if(${$tempa}[0]>$txmax ||${$tempa}[1]<$txmin){$i = 0;
 	  	  	  }else{$i = 1;last;}
 	  	  	}
 	  	   if ($i == 1){$posysv += 20;}	  	  		    
 	  	   $i = 0;
 	  	   push @temp,[$txmin,$txmax];
 	  	}
 	  } 	   
  $posy = $posy + $posysv;


  open GVR, $gvr ||die $!;
  my @recpos;
  my $lingvr = <GVR>;
  chomp $lingvr;
  my ($direct,$lstar,$lend,$lori,$rstar,$rend,$rori)= (split /\s+/,$lingvr)[0,3,4,5,7,8,9];
  my $suposy = 0;

  my $max = max($lend,$rstar);
  push @recpos, [$max,$posy];	

  while (<GVR>){
	  $posy += 16;
	  chomp;
	  ($direct,$lstar,$lend,$lori,$rstar,$rend,$rori)= (split /\s+/,$_)[0,3,4,5,7,8,9];
    my $minpos = min ($lend,$rstar);
    my $maxpos = max ($lend,$rstar);
      if($minpos <$max){      
         push @recpos,[$maxpos,$posy];
         $max = $maxpos;	
     }else{
  	     my @recposR = reverse @recpos;
  	     foreach my $record (@recposR){
  	 	   ${$record}[0] < $minpos?$posy= ${$record}[1]:last; 
  	 	 }
  	 @recpos = ();	 
  	 push @recpos,[$maxpos,$posy];
  	 $max = $maxpos;
   }
  $suposy = $posy if $posy > $suposy;
 }
 close GVR;
 return $suposy;
}
    	
sub Drawsprds{
	my ($site1,$ori1,$site2,$ori2,$posy,$ry,$rgb) = @_;
	my $min = min($site1,$site2);
	my $max = max($site1,$site2);
	my $rgba;
	if($ori1 eq $ori2){
	 	$rgba = "165,42,42";
	}elsif(($max-$min)>=500){$rgba = "106,90,205";
	}else{$rgba = "128,128,128";}
   
  $rgb ||=	$rgba;
	
	my ($starx,$endx,$rx); 
	 $starx = ($min/10) + 100;
	 $endx = ($max/10) + 100;
	 $rx = ($endx-$starx)/2;
	 print OUT "<path d=\"M$starx,$posy,A$rx,$ry,0,0,1,$endx,$posy\" style=\"fill: none; fill-opacity: 1.0; stroke: rgb($rgb); stroke-opacity: 1.0; stroke-width: 0.8; stroke-linecap: square\"/>\n";
}
