#!/usr/bin/perl

use strict;
use warnings;
use lib qw(/share/project001/yanhuang/zhenghch/perl5lib/lib64/perl5/site_perl/5.8.5/x86_64-linux-thread-multi/);
use Math::Cephes qw(&chdtrc);

package Distribute;

my $version = 0.02; my $date = '2009-4-17';

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = (); 
our @EXPORT_OK = qw( 
	chisqure normal_stat
);




sub chisqure{
	my $data = shift;

	unless (@$data){
		warn "Empty  data set, nothing calculat\n";
		return;
	}

	my $row = @$data;
	my $col = @{$data->[0]};
	
	my ($total, @r_total, @c_total);

	###calculate total
	for (my $i = 0; $i < $row; $i++){
		for (my $j = 0; $j < $col; $j++){
			$r_total[$i] += $data->[$i][$j];
			$c_total[$j] += $data->[$i][$j];
			$total += $data->[$i][$j];
		}
	}
	
	###return if total == 0 ; calculate chi
	if ($total == 0){
		warn "Total = 0\n";
		return;
	}

	my $chi = 0;
	my $flag = $row == 2 && $col == 2;
	map{
		my $line = $_;
		map{
			my $exp = $r_total[$line] * $c_total[$_] / $total;
			if ($flag){
				$chi += ( (abs($data->[$line][$_] - $exp ) - 0.5) ** 2 / $exp);
			}else{
				$chi += (abs($data->[$line][$_] - $exp) ** 2 / $exp);
			}
		} (0 .. $col - 1);
	} (0 .. $row - 1);

	my $df = ($row - 1) * ($col - 1);
	my $p_value = 'NA';
	eval {$p_value=Math::Cephes::chdtrc($df, $chi)};
	warn $@ if $@;
	return $chi, $p_value;
}



sub normal_stat{
	my $ref = shift;
	my $type = ref $ref;
	
	return () unless $type;

	my ($total, $mean, $median, $sd, $count);

	if ($type eq 'ARRAY'){
		my %len;
		$count = scalar @$ref;
		map {$len{$_}++; $total += $_;} @$ref;
		$ref = \%len;
	}elsif($type eq 'HASH'){
		map {$total += ($_ * $ref->{$_}); $count += $ref->{$_}} keys %$ref;
	}else{
		warn 'unknown type';
		return ();
	}

	$mean = $total / $count;
	##
	{	
		my $var_sum;
		map{$var_sum += ((($_ - $mean)  ** 2 ) * $ref->{$_})} keys %$ref;
		$sd = sqrt($var_sum / ($count - 1));
	}
	##
	$median = Distribute::median($ref, $count / 2);
	return ($count, $total, $mean, $median, $sd);
}


sub median{
	my ($ref, $count) = @_;
	my ($median, $total);
	for my $key (sort {$b<=>$a} keys %$ref){
		$total += $ref->{$key};
		if ($total > $count){
			$median = $key;
			last;
		}
	}
	return $median;
}

