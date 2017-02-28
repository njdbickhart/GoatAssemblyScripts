#!/usr/bin/perl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=10000
#SBATCH --output=./frc_summary_%.j.out
#SBATCH --error=./frc_summary_%.j.err
# This script is designed to automate the summarization of data from FRC combined analysis
# It is designed to run on Slurm and will process as many files as entered in the options

use strict;
use Getopt::Std;
use columnCounter;

my $usage = "perl $0 -b <comma separated list of base output file names with paths> -n <comma separated list of assembly names> -o <output file>\n";
my %opts;

getopt('bno', \%opts);

unless(defined($opts{'b'}) && defined($opts{'n'}) && defined($opts{'o'})){
	print $usage;
	exit;
}

my @batches = split(/,/, $opts{'b'});
my @names = split(/,/, $opts{'n'});

if(scalar(@batches) != scalar(@names)){
	print "Error! Number of assembly names must match the number of output file paths!\n$usage";
	exit;
}

# File subscripts: _Features.txt .lumpy.vcf .qv

my $manager = SummaryManager->new('mkdwn' => 1);
$manager->output($opts{'o'});

my $filecounts = 1;
for my ($x = 0; $x < scalar(@batches); $x++){
	my $baseName = $batches[$x];
	my $asmName = $names[$x];
	my $frcfile = "$baseName\_Features.txt";
	my $lumpyfile = "$baseName.lumpy.bedpe";
	my $qvfile = "$baseName.qv";
	# Check if the files exist
	unless( -s $frcfile && -s $lumpyfile && -s $qvfile){
		print STDERR "Missing the required files for $baseName output basename!\n";
		print STDERR "Please check to see if the following files are present:\n $frcfile\n$lumpy\n$qvfile\n";
		exit;
	}
	
	my $worker = columnCounter->new('colnum' => 1, 'mkdwn' => 1);
	$worker->ignore('#');
	
	# FRC first (it's easiest)
	$worker->readFile($frcfile);
	
	# Now Lumpy
	my $tempRand = int(rand(12323));
	open(my $OUT, "> lumpy.temp.$tempRand");
	open(my $IN, "< $lumpyfile");
	while(my $line = <$IN>){
		if($line =~/^#/){next;}
		chomp $line;
		my @segs = split(/\t/, $line);
		
		$segs[10] =~ s/TYPE://g; 
		print {$OUT} "lumpy\t$segs[10]\n";
	}
	close $IN; close $OUT;
	
	$worker->readFile("lumpy.temp.$tempRand");
	system("rm lumpy.temp.$tempRand");
	
	# Finally the QV estimate
	open(my $IN, "< $qvfile");
	my $qv = <$IN>;
	chomp $qv;
	close $IN;
	
	$worker->set("QV" => sprintf("%.2f", $qv));
	
	$manager->setTable($asmName, $worker->createSumTable($asmName));
}

$manager->PrintResults();
	

exit;
