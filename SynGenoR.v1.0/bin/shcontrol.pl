#!/usr/bin/perl -w
(@ARGV<1) && die"Usage: perl sh_contral.pl <record> [lnum,1] [shell] [sleept,100]\n";
my($rcf,$lnum,$shelf,$sleept) = @ARGV;
$lnum ||= 1;
$sleept ||= 100;
my $tem_lnum = ((-f $rcf) ? (split/\s+/,`wc -l $rcf`)[0] : 0);
until($tem_lnum>=$lnum){
	sleep($sleept);
	$tem_lnum = ((-f $rcf) ? (split/\s+/,`wc -l $rcf`)[0] : 0);
}
$shelf && (-s $shelf) && `sh $shelf`;
