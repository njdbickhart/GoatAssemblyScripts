#!/usr/bin/perl
# This script is designed to process a tab delimited file comparing Lachesis and RH map orientations of contigs
# where the RH map and Lachesis orientations disagree, the script calculates the nearest problem points
# Output is a bed file with approximate problem point locations

use strict;
use Getopt::Std;

my %opts;
my $usage = "perl $0 -i <input Lachesis tab cluster file> -f <input fasta index file .fai> -o <output bed file> -l < length bed>\n";

getopt('ifol', \%opts);
unless(defined($opts{'i'}) && defined($opts{'f'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

my $worker = SiteIdentifier->new();

# Add chr lens from fai file
open(my $IN, "< $opts{f}") || die "Could not open fai file!\n";
while(my $line = <$IN>){
	chomp $line;
	my @segs = split(/\t/, $line);
	$worker->addChrLens($segs[0], $segs[1]);
}
close $IN;

# Add Lachesis order information
$worker->assignData($opts{'i'});

# Generate problem site bed file
$worker->produceLocations($opts{'o'});

# Generate cluster contig start-end coordinate file
$worker->produceContigCoordLens($opts{'l'});

exit;

BEGIN{
package SiteIdentifier;
use Mouse;
use namespace::autoclean;

has 'chrlens' => (traits => ['Hash'], is => 'rw', isa => 'HashRef[Any]', default => sub{{}},
	handles => {
		'putc' => 'set',
		'getc' => 'get',
		'keyc' => 'keys',
	});
	
has ['chrs', 'clusters', 'lo', 'rho'] => (is => 'rw', isa => 'ArrayRef[Any]');

has 'maxlines' => (is => 'rw', isa => 'Int', default => 0);

sub addChrLens{
	my ($self, $chr, $len) = @_;
	$self->putc($chr, $len);
}

sub assignData{
	my ($self, $tabfile) = @_;
	open(my $IN, "< $tabfile") || die "Could not open tab lachesis/rh file!\n";
	my $head = <$IN>;
	my (@chrs, @clusters, @lo, @rho, $num);
	while(my $line = <$IN>){
		chomp $line;
		$line =~ s/\r//g;
		my @segs = split(/\t/, $line);
		push(@chrs, $segs[1]);
		push(@clusters, $segs[0]);
		push(@lo, $segs[5]);
		push(@rho, $segs[6]);
		$num++;
	}
	close $IN;
	
	print STDERR "Finished loading lachesis tab file\n";
	$self->chrs(\@chrs);
	$self->clusters(\@clusters);
	$self->lo(\@lo);
	$self->rho(\@rho);
	$self->maxlines($num);
}

sub produceLocations{
	my ($self, $bedfile) = @_;
	
	open(my $OUT, "> $bedfile");
	my $lastcluster = $self->clusters->[0];
	my $lastchr;
	my $curpos = 0;
	my $start = 0; my $secondstart = 0; my $error = 0;
	for(my $x = 0; $x < $self->maxlines; $x++){
		my $curchr = $self->chrs->[$x];
		my $curchrlen = $self->getc($curchr);
		my $curlo = $self->lo->[$x];
		my $currho = $self->rho->[$x];
		my $end = $curpos + $curchrlen;
		if($lastcluster != $self->clusters->[$x]){
			# The lachesis cluster changed -- print out any errors that extend to the end!
			if($error){
				print {$OUT} "cluster_$lastcluster\t$start\t$curpos\t$secondstart\n";
			}
			$curpos = 0;
			$start = 0;
			$secondstart = 0;
			$error =0;
			$lastcluster = $self->clusters->[$x];
			
		}
		
		my $rhoisproper = ($currho eq "." || $currho eq "?")? 0 : 1;
		my $rholomatch = ($curlo eq $currho)? 1 : 0;
		if($error && $rhoisproper && $rholomatch){
			# The end point of the error region
			print {$OUT} "cluster_$lastcluster\t$start\t$curpos\t$secondstart\n";
			$error = 0; $secondstart = 0;
		}elsif(!$error && $rhoisproper && !$rholomatch){
			# The beginning of the error region
			$start = $curpos;
			$error = 1;
		}elsif($error && !$rhoisproper){
			# The start of an ambiguous region
			if($secondstart == 0){
				$secondstart = $curpos;
			}
		}
		
		# 5 lowercase 'n's' were padded between each contig
		$curpos += $curchrlen + 5;
		$lastchr = $self->chrs->[$x];
	}
	
	if($error){
		print {$OUT} "cluster_$lastcluster\t$start\t$curpos\t$secondstart\n";
	}
	
	close $OUT;
}

sub produceContigCoordLens{
	my ($self, $lenfile) = @_;
	
	open(my $OUT, "> $lenfile");
	my $lastcluster = $self->clusters->[0];
	my $curpos = 0;
	for(my $x = 0; $x < $self->maxlines; $x++){
		my $curchr = $self->chrs->[$x];
		my $curchrlen = $self->getc($curchr);
		
		if($lastcluster != $self->clusters->[$x]){
			# start of a new cluster
			$curpos = 0;
			$lastcluster = $self->clusters->[$x];
		}
		
		my $end = $curpos + $curchrlen;
		print {$OUT} "cluster_$lastcluster\t$curpos\t$end\t$curchr\n";
		$curpos += $curchrlen + 5;
	}
	
	close $OUT;
}
	
	

__PACKAGE__->meta->make_immutable;
}
