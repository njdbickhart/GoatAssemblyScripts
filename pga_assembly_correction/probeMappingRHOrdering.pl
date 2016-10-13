#!/usr/bin/perl
# This script is designed to process a BAM file and the RH order file in order to assign contigs/scaffolds to chromosomes
# Adding feature to translate Brian's RH map text into an order file, directly

use strict;

chomp(@ARGV);
my $usage = "perl $0 <rh order file> <bam files with probe mappings> <output file>
	OR
perl $0 <rh order file> <output file>\n";
unless(scalar(@ARGV) == 3){
	if(scalar(@ARGV) == 2){

	}
	print $usage;
	exit;
}

my $worker = RHMapCounter->new();

# Read RH map data
$worker->getRHProbeOrder($ARGV[0]);
print STDERR "RH map order stored\n";


# Read bam file
$worker->readBAMinfo($ARGV[1]);
print STDERR "Reading bamfile!\n";

# Print output
$worker->printOrderFile($ARGV[2]);
print STDERR "Ordered entries in: $ARGV[2]\n";

exit;

BEGIN{
package RHMapCounter;
use Mouse;
use namespace::autoclean;

has 'probeOrder' => (traits => ['Array'], is => 'rw', isa => 'ArrayRef[Str]', default => sub{[]}, 
	handles => {
		'addprobe' => 'push',
		'getprobe' => 'get',
		'allprobes' => 'elements',
	});

has 'probeChr' => (traits => ['Hash'], is => 'rw', isa => 'HashRef[Any]', default => sub{{}},
	handles => {
		'setprobechr' => 'set',
		'getprobechr' => 'get',
	});
has 'probePacBio' => (traits => ['Hash'], is => 'rw', isa => 'HashRef[Any]', default => sub{{}},
	handles => {
		'setpacbio' => 'set',
		'getpacbio' => 'get',
		'existspacbio' => 'exists',
	});		
	
has 'probeIndex' => (traits => ['Hash'], is => 'rw', isa => 'HashRef[Any]', default => sub{{}},
	handles => {
		'setprobeidx' => 'set',
		'getprobeidx' => 'get',
	});

sub printOrderFile{
	my ($self, $orderfile) = @_;
	
	my @current;
	my $chrcontig = "0"; my $curchr = "0";
	open(my $OUT, "> $orderfile");
	foreach my $probe ($self->allprobes){
		if(!$self->existspacbio($probe)){next;}
		
		my $idx = $self->getprobeidx($probe);
		if($self->getpacbio($probe)->[0]  ne $chrcontig && $chrcontig ne "0" && scalar(@current) > 1){
			# new contig, and we have data to sort through
			my $values = $self->determineOrient(\@current);
			print {$OUT} "$curchr\t$chrcontig\t" . $values->[0] . "\t" . $values->[1] . "\t" . $values->[2] . "\n";
			@current = ();
		}elsif($self->getpacbio($probe)->[0] ne $chrcontig && scalar(@current) == 1){
			# new contig and there's only one entry for it!
			print {$OUT} "$curchr\t$chrcontig\t" . $current[0]->[1] . "\t" . $current[0]->[1] . "\t?\n";
			@current = ();
		}

		push(@current, $self->getpacbio($probe));
		$chrcontig = $self->getpacbio($probe)->[0];
		$curchr = $self->getprobechr($probe);
	}
	
	if(scalar(@current) > 1){
		my $values = $self->determineOrient(\@current);
		print {$OUT} "$curchr\t$chrcontig\t" . $values->[0] . "\t" . $values->[1] . "\t" . $values->[2] . "\n";
		@current = ();
	}elsif(scalar(@current) == 1){
		print {$OUT} "$curchr\t$chrcontig\t" . $current[0]->[1] . "\t" . $current[0]->[1] . "\t?\n";
		@current = ();
	}
	close $OUT;
}

sub determineOrient{
	my ($self, $array) = @_;
	my $min = 2000000000;
	my $max = 0;
	my $orient = ($array->[0]->[1] - $array->[1]->[1] < 0)? "+" : "-";
	# I'm using empty string concatenation and redundant math to ensure that the Perl interpreter
	# Doesn't just store the address of the array in memory
	#my $ctgname = $array->[0]->[5] . "";
	foreach my $v (@{$array}){
		if($v->[1] < $min){
			$min = $v->[1] + 1 - 1;
		}
		if($v->[1] > $max){
			$max = $v->[1] + 1 - 1;
		}
	}
	return [$min, $max, $orient];
}

sub readBAMinfo{
	my ($self, $bamfile) = @_;
	open(my $IN, "samtools view $bamfile |");
	my %pacbio;
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		# check if this is a secondary alignment
		if(($segs[1] & 2048) == 2048){next;}
		
		# Skip unmapped reads
		if($segs[2] eq '*' || ($segs[1] & 4) == 4){next;}
	
		$pacbio{$segs[0]} = [$segs[2], $segs[3]];
		#my $idx = $self->getprobeidx($segs[0]);
		#$self->setpacbio($segs[0] => [$segs[2], $segs[3]]);
		#$self->getprobe($idx)->pacbiopos($segs[3]);
	}
	$self->probePacBio(\%pacbio);
	close $IN;
}

sub translateRHProbeOrder{
	my ($self, $input) = @_;
	open(my $IN, "< $input") || die "Could not find input file!\n";
        my $header = <$IN>;

        local $| = 1;
        my $idx = 0;
        my @probes; my %probechrs; my %probeidx;
        while(my $line = <$IN>){
                chomp $line;
                my @segs = split(/\t/, $line);

                push(@probes, $segs[0]);
                $probechrs{$segs[0]} = $segs[1];
                $probeidx{$segs[0]} = $idx;
                #$self->addprobe($segs[0]);
                #$self->setprobechr($segs[0] => $segs[1]);
                #$self->setprobeidx($segs[0] => $idx);
                $idx++;
                if($idx % 5000 == 0){
                	print "At\t$idx\r";
                }
        }
}
	
                

sub getRHProbeOrder{
	my ($self, $input) = @_;
	open(my $IN, "< $input") || die "Could not find input file!\n";
	my $header = <$IN>;
	
	local $| = 1;	
	my $idx = 0;
	my @probes; my %probechrs; my %probeidx;
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		
		push(@probes, $segs[0]);
		$probechrs{$segs[0]} = $segs[1];
		$probeidx{$segs[0]} = $idx;
		#$self->addprobe($segs[0]);
		#$self->setprobechr($segs[0] => $segs[1]);
		#$self->setprobeidx($segs[0] => $idx);
		$idx++;
		
		if($idx % 5000 == 0){
			print "At\t$idx\r";
		}
	}

	$self->probeOrder(\@probes);
	$self->probeChr(\%probechrs);
	$self->probeIndex(\%probeidx);
	print "\n";
	close $IN;
}

sub readRHMapData{
	my($self, $input) = @_;
	# TODO: Implement this subroutine properly
	open(my $IN, "< $input") || die "Could not find input file!\n";
	my $header = <$IN>;
	my @current;
	my %data; # {chr}->[]->[contigname, min, max, orient]
	my $chrcontig = "0"; my $curchr = "0";
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		if($segs[5] =~ /^\*/){next;}
		if($segs[5] =~ /\#N\/A/){next;}

		# get the simple contig name
		my @ctgsegs = split(/\_/, $segs[5]);
		$segs[5] = $ctgsegs[0];

		if($segs[5] ne $chrcontig && $chrcontig ne "0" && scalar(@current) > 0){
			# new contig, and we have data to sort through
			push(@{$data{$curchr}}, $self->determineOrder(\@current));
			@current = ();
		}elsif($segs[5] ne $chrcontig && scalar(@current) == 0){
			# new contig and there's only one entry for it!
			push(@{$data{$curchr}}, [$chrcontig, $segs[6], $segs[6], "?"]);
		}

		push(@current, \@segs);
		$chrcontig = $segs[5];
		$curchr = $segs[1];
	}
	if(scalar(@current) >= 1){
		push(@{$data{$curchr}}, determineOrder(\@current));
	}
	close $IN;
}



__PACKAGE__->meta->make_immutable;

package RHMapUnit;
use Mouse;
use namespace::autoclean;

has ['paccontig', 'chr', 'pacorient'] => (is => 'ro', isa => 'Str', required => 1); 
has ['pacmin', 'pacmax'] => (is => 'ro', isa => 'Num', required => 1);

__PACKAGE__->meta->make_immutable;

package RHMapProbe;
use Mouse;
use namespace::autoclean;

has ['probename', 'chr'] => (is => 'ro', isa => 'Str', required => 1);
has 'pacbiochr' => (is => 'rw', isa => 'Str', default => "NA");
has 'pacbiopos' => (is => 'rw', isa => 'Num', default => 0);

__PACKAGE__->meta->make_immutable;

}
