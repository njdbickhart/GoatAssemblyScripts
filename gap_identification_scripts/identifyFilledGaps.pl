#!/usr/bin/perl
# This script is designed to take two reference genome assemblies and compare them in terms of gap content
# The first reference genome is used as a comparison point for gap identification, and then reads are generated around gap sequence
# output format:
#	Type	OChr	OGapstart	OGapend	OGaplen	SChr	Salign1coords	Salign2coords	Sfilledbases	SpercFilled	SObsGapSize
# 10/2/2016: updated logic to identify faulty alignments that were not properly mapped

use strict;
use Getopt::Std;
my %opts;

my $usage = "perl $0 -o <original reference> -s <comparative reference> -g <GetMaskBedFasta jar> -j <Java executable path> -d <output>\n";
getopt('osgjd', \%opts);

unless(defined($opts{'o'}) && defined($opts{'s'}) && defined($opts{'g'}) && defined($opts{'j'}) && defined($opts{'d'})){
	print $usage;
	exit;
}

my $worker = GapFinder->new('Oref' => $opts{'o'}, 'Sref' => $opts{'s'}, 'GetMask' => $opts{'g'}, 'Java' => $opts{'j'});

# Prep fasta
$worker->PrepareGapRegions();

# Align data
$worker->ProcessGapFQ();

# Write output
$worker->GenerateOutput($opts{'d'});

exit;

BEGIN{
package GapFinder;
use Mouse;
use MouseX::NativeTraits;
use namespace::autoclean;

my $gapBed = "temp.gap.bed";
my $gapStats = "temp.gap.stats";
my $gapFQ1 = "temp.gap.1.fq";
my $gapFQ2 = "temp.gap.2.fq";
my $pairSam = "temp.gap.sam";

# Storage -> {readname = Ochr_Ostart_Oend} -> [gaps]
has 'Storage' => (is => 'rw', isa => 'HashRef[Any]', default => sub{{}});
has 'Unmapped' => (traits => ['Hash'], is => 'rw', isa => 'HashRef[Str]', default => sub {{}}, handles => {'isunmapped' => 'exists',}); 
has ['Oref', 'Sref', 'GetMask', 'Java'] => (is => 'ro', isa => 'Str', required => 1);

sub GenerateOutput{
	my ($self, $output) = @_;
	my $store = $self->Storage;
	
	# Determine variant type so that I can generate "filled" counts
	open(my $OUT, "> $output");
	foreach my $k (sort {$a cmp $b} keys(%{$store})){
		if(!$store->{$k}->has_Type){
			$store->{$k}->DetermineType;
		}
		if($self->isunmapped($k)){
			# New logic: set to unmapped if it fit prior criteria
			$store->{$k}->Type("Unmapped");
		}else{
			$store->{$k}->DetermineFilled($self->Sref);
		}
		
		my $str = $store->{$k}->GetOutString();
		print {$OUT} $str;
	}
	close $OUT;
}

sub ProcessGapFQ{
	my $self = shift(@_);
	
	# Align gap reads with bwa
	if( -s $pairSam){
		print STDERR "Found file, processing...\n";
	}else{
		system("bwa mem " . $self->Sref . " $gapFQ1 $gapFQ2 > $pairSam");
	}
	
	print STDERR "Finished alignment!\n";
	
	# Read sam file and process
	open(my $SAM, "< $pairSam") || die "Could not open sam file!\n";
	my $store = $self->Storage;
	my $unmap = $self->Unmapped;
	while(my $line = <$SAM>){
		if($line =~ /^@/){next;}
		
		chomp $line;
		my @segs = split(/\t/, $line);
		my $readnum = (($segs[1] & 64) == 64)? 1 : 2;
		my $orient = (($segs[1] & 16) == 16)? "+" : "-";
		my $S1Start = 0; my $S1End = 0; my $S2Start = 0; my $S2End = 0;
		if(! exists($store->{$segs[0]})){
			print STDERR "Error identifying read! $segs[0]!\n";
			next;
		}
		
		#if(($segs[1] & 2048) == 2048){
			# split read alignment
			# updated logic: split read alignments indicate far too much ambiguity
		#	$unmap->{$segs[0]} = 1;
		#	next;
		#}els
		
		if($segs[2] eq "*"){
			# unmapped chr
			# updated logic: add this to the unmapped list as well
			$unmap->{$segs[0]} = 1;
			next;
		}else{
			my $cigarsoft = $self->_determineCigarSoftClip($segs[5]);
			if($cigarsoft > 50){
				# updated logic: > 10% of read length in soft clipping is unmapped
				$unmap->{$segs[0]} = 1;
			}
			my $aend = $self->_determineCigarLen($segs[3], $segs[5]);
			if($store->{$segs[0]}->has_SChr && $store->{$segs[0]}->SChr ne $segs[2]){
				$store->{$segs[0]}->Type("Trans");
				$store->{$segs[0]}->TChr($segs[2]);
				my $TOrient = ($readnum == 1)? "First" : "Second";
				$store->{$segs[0]}->TOrient($TOrient);
			}elsif(!($store->{$segs[0]}->has_SChr)){
				$store->{$segs[0]}->SChr($segs[2]);
			}

			if($readnum == 1){
				$store->{$segs[0]}->O1($orient);
				$store->{$segs[0]}->S1Start($segs[3]);
				$store->{$segs[0]}->S1End($aend);
			}else{
				$store->{$segs[0]}->O2($orient);
				$store->{$segs[0]}->S2Start($segs[3]);
				$store->{$segs[0]}->S2End($aend);
			}
		}
	}
	$self->Storage($store);
	$self->Unmapped($unmap);
	
	print STDERR "Finished sam file processing!\n";
	close $SAM;
}
sub _determineCigarSoftClip{
	my ($self, $cigar) = @_;
	my $soft = 0;
	while($cigar =~ m/(\d+)[SH]/g){
		my $count = $1;
		$soft += $count;
	}
	return $soft;
}

sub _determineCigarLen{
	my ($self, $astart, $cigar) = @_;
	my $end = $astart;
	while($cigar =~ m/(\d+)(\D{1})/g){
		my $count = $1;
		my $base = $2;
		if($base =~ m/[MDNPX\=]/){
			$end += $count;
		}
	}
	return $end;
}
	
sub PrepareGapRegions{
	my $self = shift(@_);
	# Get chr sizes
	if( -s $self->Oref . ".fai"){
		print STDERR "Generating fasta reference file...\n";
		system("samtools faidx " . $self->Oref);
	}
	my %chrsizes;
	open(my $FAI, "< " . $self->Oref . ".fai") || die "Could not find fasta index!\n";
	while(my $line = <$FAI>){
		chomp $line;
		my @segs = split(/\t/, $line);
		$chrsizes{$segs[0]} = $segs[1];
	}
	close $FAI;

	# check to see if files exist and use them
	if( -s $gapBed){
		if( -s $gapFQ1 && -s $gapFQ2){
			print STDERR "Identified all temporary files, skipping...\n";
			my $prevgapS; my $prevgapE; my $prevgapC = "NA";
			my $count = 0; my %store;

			open(my $BED, "< $gapBed") || die "Could not open Gap Bed file!\n";
			while(my $line = <$BED>){
                		chomp $line;
		                my @segs = split(/\t/, $line);
                		my $chrlen = $chrsizes{$segs[0]};
		                my $gaplen = $segs[2] - $segs[1];

                		if($segs[1] < 1000 || $segs[2] > $chrlen - 500){
                        		# Beginning of chromosome or too close to the end
		                }elsif($gaplen < 5){
                	        }elsif($prevgapC eq $segs[0] && ($segs[1] - $prevgapS < 500) && ($segs[1] - $prevgapS >= 100)){
                		        my $e2 = $segs[2] + 500;
		                        my $rn = "$segs[0]_$segs[1]_$segs[2]";
                        
					$store{$rn} = Gap->new('OChr' => $segs[0], 'OGapS' => $segs[1], 'OGapE' => $segs[2]);
					$count++;
				}else{
					my $s1 = $segs[1] - 500;
					my $e2 = $segs[2] + 500;
					my $rn = "$segs[0]_$segs[1]_$segs[2]";
					$store{$rn} = Gap->new('OChr' => $segs[0], 'OGapS' => $segs[1], 'OGapE' => $segs[2]);
					$count++;
				}
		
				$prevgapS = $segs[1];
				$prevgapE = $segs[2];
				$prevgapC = $segs[0];
			}
			
			$self->Storage(\%store);
			close $BED;
			return;
		}
	}
	
	
	# Generate bed of gaps
	system($self->Java . " -jar " . $self->GetMask . " -f " . $self->Oref . " -o $gapBed -s $gapStats");
	
	# Parse gaps and make fastq
	print STDERR "Generated gap regions. Parsing...\n";
	open(my $BED, "< $gapBed") || die "Could not open Gap Bed file!\n";
	open(my $FQ1, "> $gapFQ1");
	open(my $FQ2, "> $gapFQ2");
	my $prevgapS; my $prevgapE; my $prevgapC = "NA";
	my $count = 0; my %store;
	while(my $line = <$BED>){
		chomp $line; 
		my @segs = split(/\t/, $line);
		my $chrlen = $chrsizes{$segs[0]};
		my $gaplen = $segs[2] - $segs[1];
		
		if($segs[1] < 1000 || $segs[2] > $chrlen - 500){
			# Beginning of chromosome or too close to the end
		}elsif($gaplen < 5){
			# small gap -- not worth working on
		}elsif($prevgapC eq $segs[0] && ($segs[1] - $prevgapS < 500) && ($segs[1] - $prevgapS >= 100)){
			# gap is very close to a previous gap region but is larger than 100bp
			#my $rn1 = "$segs[0]_$prevgapS_$segs[1]";
			my $e2 = $segs[2] + 500;
			#my $rn2 = "$segs[0]_$segs[2]_$e2";
			my $rn = "$segs[0]_$segs[1]_$segs[2]";
			my $seq1 = $self->_selectRegion($segs[0], $prevgapS, $segs[1]);
			my $seq2 = $self->_selectRegion($segs[0], $segs[2], $e2);
			print {$FQ1} "\@$rn\n$seq1\n+\n" . ('I' x length($seq1)) . "\n";
			print {$FQ2} "\@$rn\n$seq2\n+\n" . ('I' x length($seq2)) . "\n";
			$store{$rn} = Gap->new('OChr' => $segs[0], 'OGapS' => $segs[1], 'OGapE' => $segs[2]);
			$count++;
		}else{
			my $s1 = $segs[1] - 500;
			#my $rn1 = "$segs[0]_$s1_$segs[1]";
			my $e2 = $segs[2] + 500;
			#my $rn2 = "$segs[0]_$segs[2]_$e2";
			my $rn = "$segs[0]_$segs[1]_$segs[2]";
			my $seq1 = $self->_selectRegion($segs[0], $s1, $segs[1]);
			my $seq2 = $self->_selectRegion($segs[0], $segs[2], $e2);
			print {$FQ1} "\@$rn\n$seq1\n+\n" . ('I' x length($seq1)) . "\n";
			print {$FQ2} "\@$rn\n$seq2\n+\n" . ('I' x length($seq2)) . "\n";
			$store{$rn} = Gap->new('OChr' => $segs[0], 'OGapS' => $segs[1], 'OGapE' => $segs[2]);
			$count++;
		}
		
		$prevgapS = $segs[1];
		$prevgapE = $segs[2];
		$prevgapC = $segs[0];
	}
	
	$self->Storage(\%store);
	close $BED;
	close $FQ1;
	close $FQ2;			
}

sub _selectRegion{
	my ($self, $chr, $start, $end) = @_;
	open(my $FA, "samtools faidx " . $self->Oref . " $chr:$start-$end | ");
	my $head = <$FA>;
	my $seq = "";
	while(my $line = <$FA>){
		chomp $line;
		$seq .= $line;
	}
	close $FA;
	return $seq;
}
		

__PACKAGE__->meta->make_immutable;

package Gap;
use Mouse;
use namespace::autoclean;

foreach my $j ('Type', 'OChr', 'SChr', 'TChr', 'TOrient', 'O1', 'O2'){ has $j => (is => 'rw', isa => 'Str', predicate => "has_$j");}
foreach my $m ('OGapS', 'OGapE', 'S1Start', 'S1End', 'S2Start', 'S2End', 'Sfilled', 'Stotal'){ has $m => (is => 'rw', isa => 'Num', predicate => "has_$m");}

sub DetermineFilled{
	my ($self, $ref) = @_;
	if($self->has_SChr && $self->has_S1Start && $self->has_S1End && $self->has_S2Start && $self->has_S2End && !($self->Type eq "Trans")){
		my $start; my $end;
		my @coords = ($self->S1Start, $self->S1End, $self->S2Start, $self->S2End);
		@coords = sort{$a <=> $b} @coords;
		
		my $len = $coords[2] - $coords[1]; 
		if($len < 5){
			print STDERR "Error assessing gap sequence: $len\n";
			$self->Sfilled(-1);
			$self->Stotal($len);
		}else{
			$self->Stotal($coords[2] - $coords[1]);
			open(my $FA, "samtools faidx $ref " . $self->SChr . ":$coords[1]-$coords[2] |");
			my $h = <$FA>;
			my $ncount = 0;
			# Count number of N's and subtract from original length to get filled percentages
			while(my $line = <$FA>){
				chomp $line;
				my ($t) = $line =~ tr/Nn/Nn/;
				$ncount += $t;
			}
			close $FA;
			$self->Sfilled($self->Stotal - $ncount);
		}		
	
	}else{
		$self->Sfilled(0);
		$self->Stotal(0);
	}
}

sub DetermineType{
	my $self = shift(@_);
	
	if(!$self->has_Type){
		if($self->has_SChr && $self->has_S1Start && $self->has_S1End && $self->has_S2Start && $self->has_S2End){
			$self->Type("Closed");
		}else{
			$self->Type("Open");
		}
	}
}

sub GetOutString{
	my $self = shift(@_);
	my $outStr = "";
	my $OGaplen = 0;
	if(!$self->has_Type){
		$self->DetermineType();
	}
	if($self->has_OGapS && $self->has_OGapE){
		$OGaplen = $self->OGapE - $self->OGapS;
	}
	my $SPerc = 0.0;
	if($self->has_Sfilled && $self->has_Stotal){
		$SPerc = $self->has_Sfilled / $self->has_Stotal;
	}
	
	$outStr = $self->Type . "\t" . $self->OChr . "\t" . $self->OGapS . "\t" . $self->OGapE . "\t" . $OGaplen . "\t";
	if($self->has_SChr && (($self->has_S1Start && $self->has_S1End) || ($self->has_S2Start && $self->has_S2End))){
		my $S1coords; my $S2coords;
		if($self->has_TChr){
			$S1coords = ($self->TOrient eq "First")? $self->TChr : $self->SChr;
			$S2coords = ($self->TOrient eq "Second")? $self->TChr : $self->SChr;
			$S1coords .= ":" . $self->S1Start . "-" . $self->S1End;
			$S2coords .= ":" . $self->S2Start . "-" . $self->S2End;
		}else{
			$S1coords = $self->SChr . ":" . $self->S1Start . "-" . $self->S1End;
			$S2coords = $self->SChr . ":" . $self->S2Start . "-" . $self->S2End;
		}
		$outStr .= $self->SChr . "\t" . $S1coords . "\t" . $S2coords . "\t" . $self->Sfilled . "\t$SPerc\t" . $self->Stotal . "\n";
	}else{
		$outStr .= "NA\tNA\tNA\t0\t0.0\t0\n";
	}
	return $outStr;
}

__PACKAGE__->meta->make_immutable;

}
