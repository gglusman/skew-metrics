#!/usr/bin/env perl
$|=1;
use strict;
my $version = "171120";

####
#
# This software computes skew metrics for a bacterial genome.
# The method is described in:
#    Lena M. Joesch-Cohen, Max Robinson, Neda Jabbari, Christopher Lausted, Gustavo Glusman.
#    Novel metrics for quantifying bacterial genome composition skews. bioRxiv. 2017
#    doi: https://doi.org/10.1101/176370
#
# Copyright 2017 by Gustavo Glusman, Institute for Systems Biology, Seattle, WA, USA.
# It is provided by the Institute for Systems Biology as open source software.
#
####
#
# The first parameter is the file holding the genome sequence in FASTA format.
# The second parameter is the file holding the genome annotation in GFF format. The sequence/contig identifiers in the GFF need to match those used in the FASTA file.
# The third (optional) parameter is the minimal length of a sequence/contig to be included. Defaults to 100 kb. Set it to 1 if you want to include all sequences.
# The fourth (optional) parameter is the fraction of sequence length to consider "too close to the end"; Ori and Ter sites closer to either end of the sequence than this cutoff are corrected to the closest sequence end. Defaults to 1%. Set it to a smaller number (e.g., 0.000000001) if you want no correction.
#
# The output gives information about the analysis run (these lines start with the # symbol), and finally gives three lines stating the computed value for the cross_product, the dot_product and the RMSD.
#
####
#
# Examples of usage:
#   skew_metrics.pl genome.fna genome.gff
#   skew_metrics.pl genome.fna genome.gff 1   (include all sequences in the FASTA file)
#
####

my($fna, $gff, $minlen, $fixends) = @ARGV;
$minlen ||= 100000;  #ignore sequences shorter than this
$fixends ||= 0.01; #set to fraction of sequence to fix, e.g. 0.01
my $minLenToFix = 1e5; #minimal length of sequence to consider fixing Ori and Ter

print "# skew_metrics.pl version $version\n";
my $lqs = getLqsParams("lqs_parameters.txt"); #file is optional

# Read sequence, compute basic statistics, identify origin of replication and terminator sites
my($seq, $desc) = readSeq($fna);
my($len, $gc, $ori, $ter, $totalLen, $globalGC, $main) = minMaxGCskew($seq);
print "# Minimal length cutoff = $minlen\n";
print "# Read $fna\n";
print "# Number of suitable contigs = ", scalar keys %$len, ", total length = $totalLen bp, GC = ", sprintf("%.1f", 100*$globalGC), "%\n";
print "# Longest sequence = $main, length = ", $len->{$main}, " bp, GC = ", sprintf("%.1f", 100*$gc->{$main}), "%\n";
print "# ", $desc->{$main}, "\n";
if ($fixends) { fixends($ori); fixends($ter) }
print "# Ori = $ori->{$main}, Ter = $ter->{$main}\n";

# Read features and segment sequences
my($feat, $totalFeat) = readGFF($gff);
print "# Read $gff\n";
print "# Number of features = $totalFeat\n";
print "# Number of features on $main = ", scalar @{$feat->{$main}}, "\n";
my $seg = segmentByFeatures($seq, $feat, $len);
print "# Number of segments on $main = ", scalar @{$seg->{$main}}, "\n";

# Compute metrics and report
my($gc, $ta, $angle, $cp, $dp, $rmsd) = compute_metrics($seg, 100*$globalGC);

print join("\t", qw/#type gc ta angle/), "\n";
foreach my $type (qw/lead lag/) {
	print join("\t",
		"#$type",
		sprintf("%.4f", $gc->{$type}),
		sprintf("%.4f", $ta->{$type}),
		sprintf("%.4f", $angle->{$type})
		), "\n";
}

print "cross_product\t", sprintf("%.4f", $cp), "\n";
print "dot_product\t", sprintf("%.4f", $dp), "\n";
print "RMSD\t", sprintf("%.4f", $rmsd), "\n";




#######################
sub getLqsParams {
	my($file) = @_;
	my %lqs;
	
	if (-e $file) {
		print "# Reading skew fit parameters from $file\n";
		open F, $file;
		$_ = <F>;
		while (<F>) {
			chomp;
			my($type, $GCbin, $intercept, $slope) = split /\t/;
			$lqs{$type}[$GCbin eq 'high'] = [$intercept, $slope];
		}
		close F;
	} else {
		## Values hard-coded here to avoid having to locate a file.
		print "# Using hard-coded skew fit parameters\n";
		$lqs{'leadGC'}[0] = [0.27622, -0.00393];
		$lqs{'leadGC'}[1] = [0.26344, -0.00379];
		$lqs{'leadTA'}[0] = [-0.23960, 0.00476];
		$lqs{'leadTA'}[1] = [0.01815, -0.00035];
		$lqs{'lagGC'}[0] = [-0.23410, 0.00501];
		$lqs{'lagGC'}[1] = [-0.04872, 0.00128];
		$lqs{'lagTA'}[0] = [0.19010, -0.00300];
		$lqs{'lagTA'}[1] = [0.07187, -0.00075];
	}
	
	return \%lqs;
}

sub compute_metrics {
	my($seg, $gc) = @_;
	
	my(%gc, %ta, %totalsize);
	foreach my $cacc (keys %$seg) {
		my $ori = $ori->{$cacc};
		my $ter = $ter->{$cacc};
		foreach my $f (@{$seg->{$cacc}}) {
			my($feature, $start, $end, $strand, $size, $gc, $ta) = @$f;
			next if length($strand)>1;
			my $where = ($start+$end)/2;
			my $type = 'lead';
			my $onBottomStrand = ($strand eq '-');
			my $invert =
				($ter>$ori && ($where<$ori || $where>$ter)) ||
				($ter<$ori && ($where<$ori && $where>$ter));
			if ($invert xor $onBottomStrand) {
				$type = 'lag' if $type eq 'lead';
			}
			my $dir = ($invert ? -1 : 1);
			$gc{$type} += $gc*$dir*$size;
			$ta{$type} += $ta*$dir*$size;
			$totalsize{$type} += $size;
		}
	}
	
	my(%angle, %magnitude, %nx, %ny);
	foreach my $type (sort keys %gc) {
		$nx{$type} = $gc{$type}/$totalsize{$type};
		$ny{$type} = $ta{$type}/$totalsize{$type};
		$angle{$type} = atan2($ny{$type}, $nx{$type});
		$magnitude{$type} = sqrt($ny{$type}**2+$nx{$type}**2);
	}
	my $GCbin = ($gc >= 50);
	my $cp = $magnitude{'lead'}*$magnitude{'lag'}*sin($angle{'lag'}-$angle{'lead'});
	my $dp = $magnitude{'lead'}*$magnitude{'lag'}*cos($angle{'lag'}-$angle{'lead'});
	my $leadGClqs = $nx{'lead'} - ($gc*$lqs->{'leadGC'}[$GCbin][1]+$lqs->{'leadGC'}[$GCbin][0]);
	my $leadTAlqs = $ny{'lead'} - ($gc*$lqs->{'leadTA'}[$GCbin][1]+$lqs->{'leadTA'}[$GCbin][0]);
	my $lagGClqs = $nx{'lag'} - ($gc*$lqs->{'lagGC'}[$GCbin][1]+$lqs->{'lagGC'}[$GCbin][0]);
	my $lagTAlqs = $ny{'lag'} - ($gc*$lqs->{'lagTA'}[$GCbin][1]+$lqs->{'lagTA'}[$GCbin][0]);
	my $rmsd = sqrt($leadGClqs**2 + $leadTAlqs**2 + $lagGClqs**2 + $lagTAlqs**2);
	
	return \%nx, \%ny, \%angle, $cp, $dp, $rmsd;
}

sub segmentByFeatures {
	my($seq, $feat, $len) = @_;
	my %seg;
	foreach my $name (sort {$len->{$b} <=> $len->{$a}} keys %$feat) {
		next if $len->{$name}<$minlen;
		my @feat = @{$feat->{$name}};
		@feat = sort {$a->{'start'} <=> $b->{'start'}} @feat;
		
		my $p;
		foreach my $i (0..$#feat) {
			my $f = $feat[$i];
			if (defined $p) {
				if ($f->{'start'} > $p->{'end'}) {
					consider($seq, \%seg, $name, $p->{'feature'}, $p->{'start'}, $p->{'end'}, $p->{'strand'});
					#consider($seq, \%seg, $name, 'gap', $p->{'end'}+1, $f->{'start'}-1, 	join('_', $p->{'strand'}, $f->{'strand'}));
				} else {
					consider($seq, \%seg, $name, $p->{'feature'}, $p->{'start'}, $f->{'start'}-1, $p->{'strand'});
					#consider($seq, \%seg, $name, 'over', $f->{'start'}, $p->{'end'}, join('=', $p->{'strand'}, $f->{'strand'}));
					$f->{'start'} = $p->{'end'}+1;
				}
			}
			$p = $f;
		}
		consider($seq, \%seg, $name, $p->{'feature'}, $p->{'start'}, $p->{'end'}, $p->{'strand'});
	}
	return \%seg;
}

sub consider {
	my($seq, $seg, $name, $feature, $start, $end, $strand) = @_;
	
	my $len = $end-$start+1;
	my $subseq = substr($seq->{$name}, $start-1, $len);
	my($A, $C, $G, $T, $rest) = basecomp($subseq);
	return unless my $binlen = ($A+$C+$G+$T);
	return unless $G+$C && $T+$A;
	
	my $gc = ($G-$C)/($G+$C);
	my $ta = ($T-$A)/($T+$A);
	push @{$seg->{$name}}, [$feature, $start, $end, $strand, $len, $gc, $ta];
}

sub readGFF {
	my($file) = @_;
	my(%feat, $total);
	if ($file =~ /\.gz$/) {
		open F, "gunzip -c $file |";
	} else {
		open F, $file;
	}
	my $name;
	while (<F>) {
		next if /^#/;
		chomp;
		my($acc, undef, $type, $start, $end, undef, $strand, $frame, $info) = split /\t/;
		next unless $type =~ /CDS|rRNA|tRNA/;
		push @{$feat{$acc}}, {'feature', $type, 'start', $start, 'end', $end, 'strand', $strand, 'info', $info};
		$total++;
	}
	close F;
	
	return \%feat, $total;
}

sub readSeq {
	my($file) = @_;
	my(%seq, %desc);
	my($name, $desc);
	if ($file =~ /\.gz$/) {
		open F, "gunzip -c $file |";
	} else {
		open F, $file;
	}
	
	while (<F>) {
		chomp;
		if (/^>(.+)/) {
			($name, $desc) = split /\s/, $1, 2;
			$desc{$name} = $desc;
		} elsif ($name) {
			s/\s+//g;
			$seq{$name} .= $_;
		} else {
			die "No name for sequence $_\n";
		}
	}
	close F;
	
	return \%seq, \%desc;
}

sub minMaxGCskew {
	my($seqs) = @_;
	my(%len, %minloc, %maxloc, %gc, $total, $globalGC);
	
	foreach my $name (keys %$seqs) {
		my $seq = $seqs->{$name};
		my $len = $len{$name} = length($seq);
		if ($len<$minlen) {
			delete $seqs->{$name};
			delete $len{$name};
			next;
		}
		
		$total += $len;
		my($max, $min, $g, $c);
		for (my $i=0;$i<$len;$i++) {
			my $w = uc substr($seq,$i,1);
			if ($w eq 'G') {
				$g++;
			} elsif ($w eq 'C') {
				$c++;
			} else {
				next;
			}
			my $d = $g-$c;
			if ($d>$max) {
				$max = $d;
				$maxloc{$name} = $i;
			} elsif ($d<$min) {
				$min = $d;
				$minloc{$name} = $i;
			}
		}
		$gc{$name} = ($g+$c)/$len{$name};
		$globalGC += ($g+$c);
	}
	
	my @sorted = sort {$len{$b}<=>$len{$a}} keys %len;
	return \%len, \%gc, \%minloc, \%maxloc, $total, $globalGC/$total, $sorted[0];
}

sub fixends {
	my($set) = @_;
	foreach my $name (keys %$seq) {
		next if $len->{$name}<$minLenToFix;
		my $loctest = $set->{$name}/$len->{$name};
		if ($loctest<$fixends) {
			$set->{$name} = 0;
		} elsif ($loctest>1-$fixends) {
			$set->{$name} = $len->{$name};
		}
	}
}

sub basecomp {
	my($seq, $caseSensitive) = @_;
	my($A, $C, $G, $T);

 	$_ = $seq;
	$_ = uc $_ unless $caseSensitive;
	$G = s/G//g;
 	$C = s/C//g;
 	$T = s/T//g;
 	$A = s/A//g;

	return ($A, $C, $G, $T, length($_));
}

