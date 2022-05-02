#!/usr/bin/perl -w

=head1 Name

   run_syngenor.pl -- a pipeline for synthetic genome restruction analysis 
                      or restructured genome analysis
     
=head1 Description
   
   This pipeline is only used for synthetic genome restruction analysis or 
   restructured genome analysis. The recombination event is induced by a 
   inducible DNA elements,such as loxp and rox,in the synthetic chromosome. 
   The inducible DNA elements are the potential breakpoint. Therefore we 
   based on split reads mapping, pair end relationship and coverage 
   information to reconstruct the restructured genome sequence.
     
   This progress consists of 6 stages: quality control(QC),reads mapping, preprocessing of 
   mapping result files, node & edge information extracting, sequence 
   reconstructing and validation. 
   
   Stage 1: Quality control
   Quality control of input reads.
   
   Stage 2: Mapping
   We map the reads to the reference (synthetic chromosome sequence) by 
   applying SOAP alignment with refined settings. Three mapping result 
   files are generated and used for subsequent relative information 
   extraction in the following analysis: pair-end mapping, single-end 
   mapping and unmapping result files. 
   
   Stage 3: Preprocessing
   We first estimate the depth distribution of generated reads with SOAP 
   Coverage(http://soap.genomics.org.cn/). We screen out all Loxp reads 
   from the three mapping result files. With all Loxp reads, we further 
   identify unmapping reads for read splitting. Then we map the split reads
   to the reference by applying bowtie2-2.0.0(http://bowtie-bio.sourceforge.
   net/index.shtml)and cluster these split reads with similar positions and
   orientations for further utilization in stage3. Pair-end Loxp reads are
   also sorted by the positions and utilized in the following stages.
   
   Stage 4: Extracting Node & Edge information
   This stage is to extract out all existing node information fromgenerated
   reads. Pair-end Loxp reads clusters are used to identify nodes without 
   recombination event happened, as supporting evidence for structure 
   variation calling in stage 5. Clusters of split reads areused to 
   identify new-formed nodes, which we define as breakpoints. We calculate 
   the number of split reads for each breakpoint and use it as a certainty 
   factor for following structure variation calling. 
 
   Each breakpoint is a signal for the connection relation between two 
   regions. So we use the node information extracted out from last stage to
   identify all existing segments, called as edges. The copy number will be
   crucial information if the variant is a deletionor duplication. So based
   on the SOAP coverage result, we further calculate the average sequencing
   depth to identify the copy number of each edge. 
   
   As our average sequencing depth is about 30x, based on our experience, 
   edges with average depth less than 10 can be identified as noise and the
   corresponding copy number value is 0.However, it is hard to directly 
   transform the average depth to copy number for random bias caused by 
   sequencing. To solve this issue, we used an iterative algorithm and 
   Mahalanobis Distance (MD) analysis to precisely estimate each edge��s copy
   number. 
   
   Stage 5: Reconstructing sequence & calling structure variation
   With all stored node and edge information, we simply treat the sequence 
   reconstruction as a path analysis. Nodes are all possible connections 
   between edges, while copy number indicates the number of each edge. So we
   use two methods, Dynamic programming based and Euler path based, to 
   reconstruct the sequence and compare the two results as mutual 
   authentication. Then we compare the reconstruction with reference to 
   obtain all structure 
   variations.
   
   Stage 6:validation
   Using the SOAPaligner2, the clean sequencing reads are aligned to the 
   restructured reference sequence. The regular mapping and irregular 
   mapping reads on new sequence based on the alignments.In theory, a perfect 
   restructured sequence, which was refactored by total variation,is mapped 
   by the sequencing reads with no fragmentation in regular mapping and no
   any clusters in irregular mapping. And all irregular pair end mapping 
   reads with variation information and split reads with loxp sequence on 
   synthetic reference can map to it normally.
      
=head1 Version
	
   Author: Wang Yun, wangyun@genomics.org.cn
           Xue Yaxin, xueyaxin@genomics.org.cn
           Shen Yue, shenyue@genomics.org.cn
            
   Version: 2.1  Date: 2016-05-31

=head1 Usage

   perl run_syngenor.pl <config_file> [-optin]
   <config>  <file>  configuration file for synthetic genome restruction analysis 
                     or restructured genome analysis.
   -show_sfware      show the software use at the process to screen.
   -step     [str]   run_syngenor step as follow, default=123456.
                     1 run quality control
                     2 run mapping by SOAP
                     3 run preprocessing of mapping result files, 
                       SOAP Coverage and split reads mapping by Bowtie
                     4 run Extracting Node & Edge information                 
                     5 run Reconstructing sequence 
                     6 run validation                             
   -run              the step for running shell, default equal -step set.
   -shell    [str]   output main shell name, default=run_syngenor.sh.
   -help             output help information to screen.  

=head1 Note
  
  1 You can use -run no, then just write the shell but not run it.
  2 At config_file options end with "\n", see the form at $Bin/config.cfg.
  3 You may not be able to get a reconstructed sequence,due to CN evaluation 
    error. You can CHECK the *.edg and change the copy number,then run step45
    again.
 
=head1 Example
 
 1 run all the progress.
   perl run_syngenor.pl config.cfg &
 2 just write the shell and then run run_syngenor.sh
   perl run_syngenor.pl config.cfg -run no
   nohup run_syngenor.sh &

=cut 
        
use strict;
use Getopt::Long;
use Cwd;
use File::Path;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin/../lib";

my  ($show_sfware,$cfg_examp,$step,$run,$shname,$help);

GetOptions(
	"show_sfware" => \$show_sfware,
	"cfg_examp"   => \$cfg_examp,
	"step:s"      => \$step,
	"run:s"       => \$run,
	"shell:s"     => \$shname,
	"help"        => \$help
);


if (@ARGV != 1 ||$help)
{
	die "
Describe: a pipeline for synthetic genome restruction analysis 
                      or restructured genome analysis
  Author: Wang Yun, wangyun\@genomics.org.cn
          Xue Yaxin, xueyaxin\@genomics.org.cn
          Shen Yue, shenyue\@genomics.org.cn
            
 Version: 1.0  Date: 2012-11-05

 Usage:perl run_syngenor.pl <config_file> [-option]
       <config>  <file>  configuration file for synthetic genome restruction analysis 
                         or restructured genome analysis.
       -show_sfware      show the software use at the process to screen.
       -cfg_examp        show config_file example, example file at $Bin/../config.cfg
       -step     [str]   run_syngenor step as follow, [default=12345 ]
                         1 run mapping by SOAP
                         2 run preprocessing of mapping result files, 
                           SOAP Coverage and split reads mapping by Bowtie
                         3 run Extracting Node & Edge information                        
                         4 run Reconstructing sequence 
                         5 run validation                             
       -run      [str]   the step for running shell, [default equal -step set ]
       -shell    [str]   output main shell name, [default=run_syngenor.sh ]
       -help             output help information to screen.  \n
       ";
}

### option value ######
$step ||= 123456;
$run ||= $step;
my $cfgf = shift;
(-s $cfgf) || die "error: can't find available config file 
   $cfgf $!";
$shname ||= "run_syngenor.sh";

my %sfware = (
          "QCfilter"    => "perl $Bin/filter_gz.pl",
	  "statreads"   => "perl $Bin/statreads.pl",
 	  "pickunmap"   => "perl $Bin/pickunmap.pl",
 	  "splitreads"  => "perl $Bin/splitreads.pl",
 	  "getvars"     => "perl $Bin/getvariareads.pl",
 	  "sortool"     => "perl $Bin/sortool.pl",
 	  "extNaE"      => "perl $Bin/extnodedge.pl",
 	  "dynpath"     => "perl $Bin/dynpath2recon.pl",
 	  "restructseq" => "perl $Bin/restructseq.pl",
 	  "SVGstdchc"   => "perl $Bin/drawSVGstdchc.pl",
 	  "SVGlxpchc"   => "perl $Bin/drawSVGlxpchc.pl",
 	  "shcontrol"   => "perl $Bin/shcontrol.pl",
 	  "tfcoverage"  => "perl $Bin/tfcoverage.pl",
 	  "getpereads"  => "perl $Bin/getpereads.pl",
 	  "classifyspr" => "perl $Bin/classifyspr.pl",
 	  "lxpspstat"   => "perl $Bin/lxpspstat.pl",
 	  "unmapcheck"  => "perl $Bin/unmapcheck.pl",
 	  "fdloxpftr"   => "perl $Bin/fdloxpftr.pl",
 	);

if($cfg_examp){system"less $Bin/../config.cfg";exit;} 
  
##step0 prepare check####
my %cfgh = get_cfg($cfgf,'cfg');

foreach (keys %cfgh){
	next if $cfgh{$_} =~ /qsub/;
  if($cfgh{$_} =~ /\/\w+\//){
    if($_ =~ /tools/){
    (-e $cfgh{$_}) || die "error:No executable alignment tools:
  $cfgh{$_} at \"$_\" option in config file\n";
    }else{
    (-d $cfgh{$_} ||-s $cfgh{$_}) || die "error:$! :
  $cfgh{$_} at \"$_\" option in config file\n";   
    }
   }
 }

foreach(keys %sfware){
  my $program = (split /\s+/,$sfware{$_})[1];	
	(-s $program)|| die "error:absence of such software:
	$program \n";
 }

if($show_sfware){
   foreach(keys %sfware){
    print "$sfware{$_}\n" ;
   }
  print 
"#alignment tools:
SOAP2 $cfgh{soaptools}
Bowtie2 $cfgh{bowtietools}
BWA $cfgh{bwatools}
Soap.coverage $cfgh{soapcvgtools}\n";
}


##mkdir for each step###
my $temdir = "$cfgh{outdir}/$cfgh{sample_ID}";
(-d $temdir) || mkdir ($temdir);

my ($dir1,$dir2,$dir31,$dir32,$dir33,$dir34,$dir4,$dir5,$dir61,$dir62,$dir63,$dir64)=
qw(01.cleandata 02.soap 03.prep/01.coverage 03.prep/02.readstat 03.prep/03.bowtie 03.prep/04.sortgvr 04.extNaE 05.dynpath 06.valid/01.valsoap 06.valid/02.stdchc 06.valid/03.umapBWA 06.valid/04.lxpchc);
foreach($dir1,$dir2,$dir31,$dir32,$dir33,$dir34,$dir4,$dir5,$dir61,$dir62,$dir63,$dir64){$_ = "$temdir/$_";}


##======================================================##
##         To make directory and write shell            ##
##======================================================##

open OPG,">$temdir/$shname" || die $!;

my %script;
my $sucess = "echo \"work completed successfully !!\"";

##step1 Quality QC##
if($step =~ /1/ && $cfgh{stage1}){
	 mdir($dir1);
   my $QCshell;
   $QCshell = "$sfware{QCfilter} $cfgh{raw_fq_1} $cfgh{raw_fq_2} -prefix $cfgh{sample_ID} $cfgh{filter_option} ";
   
   $QCshell .= "-trim -seed $cfgh{seed} " if $cfgh{trim};
   $QCshell .= "-deadapter -adapterlist1 $cfgh{adapterlist1} -adapterlist2 $cfgh{adapterlist2} " if $cfgh{deadapter};
   $QCshell .= "-outSE " if $cfgh{outSE};
   $QCshell .= "-outfilted " if $cfgh{outfilted};
   $QCshell .= "-outdir $dir1\n\n";
   if ($cfgh{dedup}){
      $QCshell .= "gzip -d $cfgh{sample_ID}\_1.clean.fq.gz\n";
      $QCshell .= "gzip -d $cfgh{sample_ID}\_2.clean.fq.gz\n\n";
      $QCshell .= "$cfgh{deduptools} $cfgh{sample_ID}\_1.clean.fq $cfgh{sample_ID}\_2.clean.fq $cfgh{sample_ID}\_1.clean.dedup.fq $cfgh{sample_ID}\_2.clean.dedup.fq $cfgh{sample_ID}\.clean.dedup.stat\n\n";
      $QCshell .= "gzip $cfgh{sample_ID}\_1.clean.dedup.fq\n";
      $QCshell .= "gzip $cfgh{sample_ID}\_2.clean.dedup.fq\n\n";   
   }
   
   $QCshell .= "$sucess >QCfilter.rec\n";	
   mfile("$dir1/QCfilter.sh","$QCshell");
   mfile("$dir1/work.sh","$cfgh{QCfilter_qsub} QCfilter.sh\n$sfware{shcontrol} QCfilter.rec\n");  
   
   shell_bos($dir1,"step1 Quality control (filtering)","work.sh",);
}

  
##step2 mapping##
#my $timemaker = "date | awk '{print \$0\" -> star $progress\"}'";
my ($cleanfq1,$cleanfq2);

if ($cfgh{stage1} && $cfgh{dedup}){
   $cleanfq1 = "$dir1/$cfgh{sample_ID}\_1.clean.dedup.fq.gz";
   $cleanfq2 = "$dir1/$cfgh{sample_ID}\_2.clean.dedup.fq.gz";
}elsif($cfgh{stage1}){
   $cleanfq1 = "$dir1/$cfgh{sample_ID}\_1.clean.fq.gz";
   $cleanfq2 = "$dir1/$cfgh{sample_ID}\_2.clean.fq.gz";  
}elsif(defined $cfgh{fq_1} &&  defined $cfgh{fq_2}){
   $cleanfq1 = $cfgh{fq_1};
   $cleanfq2 = $cfgh{fq_2};
}else{
   die "Not Fould input files $! \n";
}

if($step =~ /2/ && $cfgh{stage2}){
  mdir($dir2);
  my $soap_2bwt = dirname($cfgh{soaptools});  
  my $soapshell;
  unless ($cfgh{soap_index}==0 || (-s "$cfgh{soap_ref}.index.bwt")){$soapshell = "$soap_2bwt/2bwt-builder $cfgh{soap_ref}\n";}
  $soapshell .= "$cfgh{soaptools} -a $cleanfq1 -b $cleanfq2 $cfgh{soap_option} -o $dir2/$cfgh{sample_ID}.soap -2 $dir2/$cfgh{sample_ID}.soap2 -D $cfgh{soap_ref}.index\n"; 
  $soapshell .= "gzip $dir2/$cfgh{sample_ID}.soap\n";
  $soapshell .= "gzip -c $dir2/$cfgh{sample_ID}.soap2 >$dir2/$cfgh{sample_ID}.soap2.gz\n";
  $soapshell .= "$sucess >soap.rec\n";		
  mfile("$dir2/soap.sh","$soapshell");
  mfile("$dir2/work.sh","$cfgh{soap_qsub} soap.sh\n$sfware{shcontrol} soap.rec\n");  
  
  shell_bos($dir2,"step2 mapping","work.sh");
}

##step2 Preprocessing##

if($step =~ /3/ && $cfgh{stage3}){
	my $prepdir = "$temdir/03.prep";
	
	##[s1 soapcoverage]
	mdir($prepdir);
	mdir($dir31);
  my $cvgshell = "$cfgh{soapcvgtools} $cfgh{soap_coverage_option} -refsingle $cfgh{soap_ref} -i $dir2/$cfgh{sample_ID}.soap.gz $dir2/$cfgh{sample_ID}.soap2.gz -depthsingle $dir31/$cfgh{sample_ID}.coverage.depthsingle -o $dir31/soappaired.out\n";		
	$cvgshell .= "$sucess >soapcvg.rec\n";	
	mfile("$dir31/soapcvg.sh","$cvgshell");
	
	#[s2 read statistic]
	mdir($dir32);
  my $readstatshell = "$sfware{statreads} $cleanfq1 $cleanfq2 $dir2/$cfgh{sample_ID}.soap.gz $dir2/$cfgh{sample_ID}.soap2.gz -spid $cfgh{sample_ID}  -syntype $cfgh{syn_type} -chrsyn $cfgh{chr_syn} -synlen $cfgh{syn_len} -chrwt $cfgh{chr_wt} -wtlen $cfgh{wt_len} -reflen $cfgh{ref_len} -DNAe $cfgh{DNAe} -o $dir32\n";
  $readstatshell .= "$sucess >readstat.rec\n";
  mfile("$dir32/readstat.sh","$readstatshell");
  
  #[s3 splitting reads & Bowtie alignment ]
	mdir($dir33);
  my $bowtie_ref_dir = dirname($cfgh{syn_seq});
  my $bowtie_ref_name= basename($cfgh{syn_seq});
  $bowtie_ref_name = (split /\./,$bowtie_ref_name)[0];
  #print "check ::bowtie_ref_dir--$bowtie_ref_dir\nbowtie_ref_name--$bowtie_ref_name\n";
  my $bowtieshell = "## s3.1 split reads ##
$sfware{pickunmap} $dir32/$cfgh{sample_ID}.loxp.fq $dir32/$cfgh{sample_ID}.loxp.sp -spid $cfgh{sample_ID} -outdir $dir33\n";
  $bowtieshell .= "$sfware{splitreads} $dir33/$cfgh{sample_ID}.umap.fq -spid $cfgh{sample_ID} -o $dir33\n";
  $bowtieshell .= "## s3.2 Bowtie alignment ##\n";
  unless ($cfgh{bowtie_index}== 0 || (-s "$bowtie_ref_dir/$bowtie_ref_name.1.bt2")){$bowtieshell .= "$cfgh{bowtietools}\-build $cfgh{syn_seq} $bowtie_ref_dir/$bowtie_ref_name\n";}  
  $bowtieshell .= "$cfgh{bowtietools} -x $bowtie_ref_dir/$bowtie_ref_name -q $dir33/$cfgh{sample_ID}.splitted.fq -S $dir33/$cfgh{sample_ID}.splitted.map\n";
  $bowtieshell .= "$sucess >Bowtie.rec\n";
  mfile("$dir33/Bowtie.sh","$bowtieshell");
  
  #[s4 sorting split reads mapping]
  mdir($dir34);
  my $sortgvrshell = "$sfware{getvars} $dir33/$cfgh{sample_ID}.splitted.map -spid $cfgh{sample_ID} -outdir $dir34\n";
  $sortgvrshell .= "sort -k 5g,5 -k 9gr,9 $dir34/$cfgh{sample_ID}.gvr > $dir34/$cfgh{sample_ID}.gvr.temp\n";
  $sortgvrshell .= "$sfware{sortool} $dir34/$cfgh{sample_ID}.gvr.temp -spid $cfgh{sample_ID} -type gvr -outdir $dir34\n";
  $sortgvrshell .= "$sucess >sortgvr.rec\n";
  mfile("$dir34/sortgvr.sh","$sortgvrshell");
  ##
  mfile("$prepdir/work.sh","##[s1 soapcoverage]
cd $dir31\n$cfgh{soapcvg_qsub} soapcvg.sh\n$sfware{shcontrol} soapcvg.rec
#[s2 read statistic]
cd $dir32\nsh readstat.sh >readstat.log
#[s3 split reads & Bowtie alignment ]
cd $dir33\n$cfgh{bowtie_qsub} Bowtie.sh\n$sfware{shcontrol} Bowtie.rec
#[sorting split reads mapping]  
cd $dir34\nsh sortgvr.sh > sortgvr.log
");

  shell_bos($prepdir,"step3 Preprocessing","work.sh");
}

##step3 Extracting Node & Edge (extNaE) information## 
if($step =~ /4/ && $cfgh{stage4}){
	mdir($dir4);
	my $extNaEshell = "$sfware{extNaE} -prefix $cfgh{sample_ID} -chrid $cfgh{chr_syn} -coverage $dir31/$cfgh{sample_coverage} -refcover $cfgh{ref_cover} -abnorm $dir34/$cfgh{sample_ID}.sort.gvr -soap $dir2/$cfgh{soap}.gz -soap2 $dir2/$cfgh{soap2}.gz -loxpregion $cfgh{loxp_region} $cfgh{extNaE_option} -outdir $dir4 >$dir4/$cfgh{sample_ID}.log\n";
	$extNaEshell .= "$sucess >extNaE.rec\n"; 
  mfile("$dir4/extNaE.sh","$extNaEshell");	  
	mfile("$dir4/work.sh","cd $dir4\n$cfgh{extNaE_qsub} extNaE.sh\n$sfware{shcontrol} extNaE.rec\n");  
  shell_bos($dir4,"step4 Extracting Node & Edge (extNaE) information","work.sh");  
}
	
##step5 dynamic programming path(Dynpath) to reconstruct sequence ##
if($step =~ /5/ && $cfgh{stage5}){
  mdir($dir5);
  my $dynpathshell = "$sfware{dynpath} $dir4/$cfgh{sample_ID}.edg $dir4/$cfgh{sample_ID}.nod  -spid $cfgh{sample_ID} -mcf $cfgh{max_confict} -outdir $dir5\n";
  $dynpathshell .= "$sucess >Dynpath.rec\n"; 
  mfile("$dir5/Dynpath.sh","$dynpathshell");  
  mfile("$dir5/work.sh","cd $dir5\n$cfgh{dynpath_qsub} Dynpath.sh\n$sfware{shcontrol} Dynpath.rec\n");  
  shell_bos($dir5, "step5 Dynpath to reconstruct sequence", "work.sh"); 	
}

##step6 validation ##
if($step =~ /6/ && $cfgh{stage6}){
	my $validir =  "$temdir/06.valid";
  
  #[s1 validation SOAP ]
  mdir($validir);
  mdir($dir61);  
  $cfgh{val_scale} =~ /first/i || print "\"val_scale = all\" option is not available.The first path is used for validation\n";
  
  my $bkgname = basename($cfgh{bkgseq});
  $bkgname = (split /\./,$bkgname)[0];
  my $soap_2bwt = dirname($cfgh{soaptools}); 
  mdir ("$dir61/index"); 
  my $workshell = "cd $dir61\n";
  $workshell .= "$sfware{restructseq} $cfgh{syn_seq} $cfgh{loxp_region} $dir5/$cfgh{sample_ID}.encod -spid $cfgh{sample_ID} -reseqid $cfgh{reseqid} -outdir $dir61\n";
  $workshell .= "cat $cfgh{bkgseq} $dir61/$cfgh{sample_ID}.reseq >./index/$bkgname\_scb.fa\n";
  $workshell .= "cd $dir61/index\n";
  $workshell .= "$soap_2bwt/2bwt-builder $bkgname\_scb.fa\n";
  
  my $valspshell = "$cfgh{soaptools} -a $cleanfq1 -b $cleanfq2 $cfgh{soap_option} -o $dir61/$cfgh{valsp} -2 $dir61/$cfgh{valsp2} -D $dir61/index/$bkgname\_scb.fa.index\n";
  $valspshell .= "$sucess >valsoap.rec\n";	
  mfile("$dir61/valsoap.sh","$valspshell");
    
  $workshell .= "cd $dir61\n$cfgh{soap_qsub} valsoap.sh\n$sfware{shcontrol} valsoap.rec\n";  
  
  #[s2  standard mapping check ]
  if($cfgh{val_type} =~ /stdchc|both/i){
     mdir($dir62);
     my $stdchcshell = "##data preparation\n";
     $stdchcshell .= "$sfware{tfcoverage} $dir31/$cfgh{sample_coverage} -spid $cfgh{sample_ID} -chrid $cfgh{chr_syn} -wind $cfgh{tf_wind} -outdir $dir62\n";
     $stdchcshell .= "awk '\$8==\"$cfgh{chr_syn}\" {print \$0}' $dir2/$cfgh{soap2} > $dir62/$cfgh{sample_ID}.soap2\n";
     $stdchcshell .= "$sfware{getpereads} $dir62/$cfgh{sample_ID}.soap2 -spid $cfgh{sample_ID} $cfgh{getpe_option} -refseq $cfgh{syn_seq} -outfile $dir62/$cfgh{sample_ID}.spr\n";
     $stdchcshell .= "$sfware{sortool} $dir62/$cfgh{sample_ID}.spr -spid $cfgh{sample_ID} -type spr -outdir $dir62\n";
     $stdchcshell .= "$sfware{classifyspr} $dir34/$cfgh{sample_ID}.mpi $dir62/$cfgh{sample_ID}.sort.spr -spid $cfgh{sample_ID} -insertsize $cfgh{insertsize} -outdir $dir62\n";
     $stdchcshell .= "awk '\$8==\"$cfgh{chr_syn}\" {print \$0}' $dir2/$cfgh{soap2} $dir2/$cfgh{soap} > $dir62/$cfgh{sample_ID}.sp.temp\n";
     $stdchcshell .= "$sfware{lxpspstat} $dir62/$cfgh{sample_ID}.sp.temp $cfgh{loxp_ftr} -prefix $cfgh{sample_ID} -DNAe $cfgh{DNAe} -readlen $cfgh{readlen} -outdir $dir62\n";
     $stdchcshell .= "rm $dir62/$cfgh{sample_ID}.sp.temp\n";
     $stdchcshell .= "##validation soap aligment data preparation\n";
     $stdchcshell .= "awk '\$8==\"$cfgh{reseqid}\" {print \$0}' $dir61/$cfgh{valsp2} > $dir62/$cfgh{sample_ID}.val.soap2\n";
     $stdchcshell .= "$sfware{getpereads} $dir62/$cfgh{sample_ID}.val.soap2 -spid $cfgh{sample_ID} $cfgh{getpe_option} -refseq $dir61/$cfgh{sample_ID}.reseq -outfile $dir62/$cfgh{sample_ID}.val2.spr\n";
     $stdchcshell .= "$sfware{sortool} $dir62/$cfgh{sample_ID}.val2.spr -spid $cfgh{sample_ID} -type spr -outfile $dir62/$cfgh{sample_ID}.sort.val2.spr\n";
     $stdchcshell .= "##\n";
     $stdchcshell .= "awk '\$8==\"$cfgh{reseqid}\" {print \$0}' $dir61/$cfgh{valsp} > $dir62/$cfgh{sample_ID}.val.soap\n";
     $stdchcshell .= "$sfware{getpereads} $dir62/$cfgh{sample_ID}.val.soap -spid $cfgh{sample_ID} -refseq $dir61/$cfgh{sample_ID}.reseq -outfile $dir62/$cfgh{sample_ID}.val.spr\n";
     $stdchcshell .= "$sfware{sortool} $dir62/$cfgh{sample_ID}.val.spr -spid $cfgh{sample_ID} -type spr -outfile $dir62/$cfgh{sample_ID}.sort.val.spr\n";
     $stdchcshell .= "## draw std svg file\n";
     $stdchcshell .= "$sfware{SVGstdchc} -prefix $cfgh{sample_ID} -insize $cfgh{insertsize} -refid $cfgh{chr_syn} -reflen $cfgh{syn_len} -reseq  $dir61/$cfgh{sample_ID}.reseq -loxpftr $cfgh{loxp_ftr} -cvg $dir62/$cfgh{sample_ID}.depth -lxps $dir62/$cfgh{sample_ID}.lxps -gvr $dir34/$cfgh{sample_ID}.sort.gvr -spr $dir62/$cfgh{sample_ID}.sort.spr -varid $dir62/$cfgh{sample_ID}.varid -val2 $dir62/$cfgh{sample_ID}.sort.val2.spr -val $dir62/$cfgh{sample_ID}.sort.val.spr -lsv $dir62/$cfgh{sample_ID}.lsv -outdir $dir62\n";  
     
     mfile("$dir62/stdchc.sh","$stdchcshell");
     my $svgshell = "$sfware{SVGstdchc} -prefix $cfgh{sample_ID} -insize $cfgh{insertsize} -refid $cfgh{chr_syn} -reflen $cfgh{syn_len} -reseq  $dir61/$cfgh{sample_ID}.reseq -loxpftr $cfgh{loxp_ftr} -cvg $dir62/$cfgh{sample_ID}.depth  -lxps $dir62/$cfgh{sample_ID}.lxps -gvr $dir34/$cfgh{sample_ID}.sort.gvr -spr $dir62/$cfgh{sample_ID}.sort.spr -varid $dir62/$cfgh{sample_ID}.varid -val2 $dir62/$cfgh{sample_ID}.sort.val2.spr -val $dir62/$cfgh{sample_ID}.sort.val.spr -lsv $dir62/$cfgh{sample_ID}.lsv -outdir $dir62\n";  
     mfile("$dir62/drawsvg.sh","$svgshell");
     
     $workshell .= "cd $dir62\n sh stdchc.sh >stdchc.log\n";       
   }
  if($cfgh{val_type} =~ /lxpchc|both/i){
     mdir($dir63);
     mdir($dir64);
     mdir("$dir63/index");
   	 my $lxpchcshell = "##data preparation\n";
   	 $lxpchcshell .= "$sfware{getpereads} $dir32/$cfgh{sample_ID}.loxp.sp -spid $cfgh{sample_ID} -outfile $dir64/$cfgh{sample_ID}.loxp.spr\n";
     $lxpchcshell .= "$sfware{sortool} $dir64/$cfgh{sample_ID}.loxp.spr -spid $cfgh{sample_ID} -type spr -outfile $dir64/$cfgh{sample_ID}.loxp.sort.spr\n";
     $lxpchcshell .= "$sfware{lxpspstat} $dir32/$cfgh{sample_ID}.loxp.sp  $cfgh{loxp_ftr} -prefix $cfgh{sample_ID} -DNAe $cfgh{DNAe} -readlen $cfgh{readlen} -splxpid  -outdir $dir64\n";
     $lxpchcshell .= "##..\n";
     $lxpchcshell .= "$sfware{statreads} $cleanfq1 $cleanfq2 $dir61/$cfgh{valsp} $dir61/$cfgh{valsp2} -spid $cfgh{sample_ID}  -syntype $cfgh{syn_type} -chrsyn $cfgh{chr_syn} -synlen $cfgh{syn_len} -chrwt $cfgh{chr_wt} -wtlen $cfgh{wt_len} -reflen $cfgh{ref_len} -DNAe $cfgh{DNAe} -o $dir64\n";
     $lxpchcshell .= "mv $dir64/$cfgh{sample_ID}.loxp.sp $dir64/$cfgh{sample_ID}.val.loxp.sp\n";
     $lxpchcshell .= "$sfware{pickunmap} $dir64/$cfgh{sample_ID}.loxp.fq  $dir64/$cfgh{sample_ID}.val.loxp.sp -spid $cfgh{sample_ID}.val -outdir $dir64\n";     
     
     my $umBWAshell = "ln -s $dir61/index/$bkgname\_scb.fa $dir63/index/$bkgname\_scb.fa\n";
     $umBWAshell .= "$cfgh{bwatools} index $dir63/index/$bkgname\_scb.fa\n";
     $umBWAshell .= "$cfgh{bwatools} aln $cfgh{bwa_option} $dir63/index/$bkgname\_scb.fa $dir64/$cfgh{sample_ID}.val.umap.fq >$dir63/$cfgh{sample_ID}.val.umap.sai\n";
     $umBWAshell .= "$cfgh{bwatools} samse $dir63/index/$bkgname\_scb.fa $dir63/$cfgh{sample_ID}.val.umap.sai $dir64/$cfgh{sample_ID}.val.umap.fq >$dir63/$cfgh{sample_ID}.val.umap.fai\n";
     $umBWAshell .= "sort -k 3,3 -k 4g,4 $dir63/$cfgh{sample_ID}.val.umap.fai >$dir63/$cfgh{sample_ID}.val.umap.sort.fai\n";
     mfile("$dir63/ummapbwa.sh","$umBWAshell");
     
     $lxpchcshell .= "sh $dir63/ummapbwa.sh\n$sucess >ummapbwa.rec\n";
     $lxpchcshell .= "$sfware{unmapcheck}  $dir64/$cfgh{sample_ID}.val.umap.fq $dir63/$cfgh{sample_ID}.val.umap.sort.fai -prefix $cfgh{sample_ID} -bwachc  $cfgh{unmapcheck} -outdir $dir64\n";
     $lxpchcshell .= "$sfware{fdloxpftr} $dir61/$cfgh{sample_ID}.reseq $dir64/$cfgh{sample_ID}.lxpftr.lst -DNAe loxp\n";
     $lxpchcshell .= "$sfware{lxpspstat} $dir64/$cfgh{sample_ID}.val.loxp.sp  $dir64/$cfgh{sample_ID}.lxpftr.lst -prefix $cfgh{sample_ID}.val -DNAe $cfgh{DNAe} -readlen $cfgh{readlen} -splxpid  -outdir $dir64\n";
     $lxpchcshell .= "$sfware{getpereads} $dir64/$cfgh{sample_ID}.val.loxp.sp -spid $cfgh{sample_ID}.val -outfile $dir64/$cfgh{sample_ID}.val.loxp.spr\n";
     $lxpchcshell .= "$sfware{sortool} $dir64/$cfgh{sample_ID}.val.loxp.spr -spid $cfgh{sample_ID}.val -type spr -outfile $dir64/$cfgh{sample_ID}.val.loxp.sort.spr\n";
     $lxpchcshell .= "$sfware{unmapcheck}  $dir64/$cfgh{sample_ID}.umc $dir63/$cfgh{sample_ID}.splitted.map $dir64/$cfgh{sample_ID}.splitted.umc -prefix $cfgh{sample_ID} -bowtiechc  \n";
     $lxpchcshell .= "$sfware{SVGlxpchc} -prefix $cfgh{sample_ID} -insize $cfgh{insertsize} -DNAtype $cfgh{DNAtype} -refid $cfgh{chr_syn} -refcode $cfgh{syn_code} -refseq $cfgh{syn_seq} -reflen $cfgh{syn_len}  -loxpftr  $cfgh{loxp_ftr} -lxpfq $dir64/$cfgh{sample_ID}.loxp.fq -lxpspr $dir64/$cfgh{sample_ID}.loxp.sort.spr -lxpspid $dir64/$cfgh{sample_ID}.loxp.spid -lxpspstat $dir64/$cfgh{sample_ID}.lxps -rencode  $dir61/$cfgh{sample_ID}.encod -reseq $dir61/$cfgh{sample_ID}.reseq -relxpftr $dir64/$cfgh{sample_ID}.lxpftr.lst -relxpspr $dir64/$cfgh{sample_ID}.val.loxp.sort.spr -relxpspid $dir64/$cfgh{sample_ID}.val.loxp.spid -relxpspstat $dir64/$cfgh{sample_ID}.val.lxps -reumc $dir64/$cfgh{sample_ID}.umc -splitumc $dir64/$cfgh{sample_ID}.splitted.umc -outdir  $dir64\n";
     mfile ("$dir64/lxpchc.sh","$lxpchcshell");
     my $svgshell = "$sfware{SVGlxpchc} -prefix $cfgh{sample_ID} -insize $cfgh{insertsize} -DNAtype $cfgh{DNAtype} -refid $cfgh{chr_syn} -refcode $cfgh{syn_code} -refseq $cfgh{syn_seq} -reflen $cfgh{syn_len}  -loxpftr  $cfgh{loxp_ftr} -lxpfq $dir64/$cfgh{sample_ID}.loxp.fq -lxpspr $dir64/$cfgh{sample_ID}.loxp.sort.spr -lxpspid $dir64/$cfgh{sample_ID}.loxp.spid -lxpspstat $dir64/$cfgh{sample_ID}.lxps -rencode  $dir61/$cfgh{sample_ID}.encod -reseq $dir61/$cfgh{sample_ID}.reseq -relxpftr $dir64/$cfgh{sample_ID}.lxpftr.lst -relxpspr $dir64/$cfgh{sample_ID}.val.loxp.sort.spr -relxpspid $dir64/$cfgh{sample_ID}.val.loxp.spid -relxpspstat $dir64/$cfgh{sample_ID}.val.lxps -reumc $dir64/$cfgh{sample_ID}.umc -splitumc $dir64/$cfgh{sample_ID}.splitted.umc -outdir  $dir64\n";
     mfile ("$dir64/drawsvg.sh","$svgshell");
     
     $workshell .= "cd $dir64\n sh lxpchc.sh >lxpchc.log\n";      
   }
     
   mfile("$validir/work.sh","$workshell");
   shell_bos($validir, "step6 Validation", "work.sh"); 	   
  }


##======================================================##
##                   To run the shell                   ##
##======================================================##
if($run !~ /\D/){
	if($run eq $step){
		system"sh $temdir/$shname";
	}else{
		open SESH, ">$temdir/$shname.se.sh"||die $!;
		foreach (1..5){	
			if($run =~ /$_/){
		   my $key = "step$_";		  
		   print SESH "$script{$key}";	 
		  }
		}
			system"sh $temdir/$shname.se.sh";
	}
}  
  
#######################################################  
##     ==       ==   SUB FUNCTION   ==       ==      ##
#######################################################

##================##
sub get_cfg
##================##
{
	my $contxt = shift;
	my $get_line = 'awk \'(!/^#/)\'';
  my @arry = split /\s*\n+/,`$get_line $contxt`;
  my %path;
	shift @arry;
	foreach(@arry) {
	  my ($a,$b)	= split/\s+=\s+/,$_;
		$path{$a} = $b;
	}
  return %path
}

##================##
sub mdir
##================##
{
   (-d $_[0]) || `mkdir "$_[0]"`;
}

##================##
sub mfile
##================##
{
    $_[2] ? (open FILE,">>$_[0]" || die$!) : (open FILE,">$_[0]" || die$!);
    print FILE $_[1];
    $_[2] && print FILE $_[2];
    close FILE;
}

##================##
sub shell_bos
##================##
{
  my ($dir,$progress,$shell) = @_;
  my $stage = (split /\s+/,$progress)[0];
  $stage =~ s/(\d+)/$1/; 
  $script{$stage} = "##$progress
date | awk '{print \$0\" -> star $progress\"}'
cd $dir\nsh $shell
date | awk '{print \$0\" -> finish $progress\"}'\n\n";

  print OPG "$script{$stage}";
}
