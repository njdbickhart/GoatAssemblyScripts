#!/usr/bin/perl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=2000
#SBATCH --output=./frc_analysis_%.j.out
#SBATCH --error=./frc_analysis_%.j.err
# This script is designed to automate the analysis of all FRC tools for assembly validation
# Assumes that the job submission system is Slurm and that the system has module control

use strict;
use Getopt::Std;
use File::Basename;
use lib dirname (__FILE__);
use slurmTools;
my %opts;
my $usage = "perl $0 -f <assembly fasta> -g <genome size (bp)> -o <output basename> -b <bams, comma separated>\n";

getopt('fgob', \%opts);
unless(defined($opts{'f'}) && defined($opts{'g'}) && defined($opts{'o'}) && defined($opts{'b'})){
	print $usage;
	exit;
}

my @bams = split(/,/, $opts{'b'});
my $outPrefix = $opts{'o'};
my $inputBam;
if(scalar(@bams) > 1){
	# merge split bams so that downstream processes can use them
	print "Generating merged bam file...\n";
	$inputBam = "$opts{o}.merged.total.bam";
	my $bamWhiteSpace = join(" ", @bams);
	system("module load samtools; samtools merge -c -p --threads=4 $inputBam $bamWhiteSpace; samtools index $inputBam");
	unless( -s $inputBam){
		print STDERR "Error generating merged bam file! Attempted merge filename: $inputBam. Input bams:\n$bamWhiteSpace\n";
		exit;
	}
}else{
	$inputBam = $bams[0];
}

my @modules = ('lumpy-sv/0.2.12-51-g16b6876', 'FRC_align/1.3.0-5b3f53e', 'freebayes/1.1.0-1-gf15e66e');
my $sideWorker = slurmTools->new("workDir" => "./", "scriptDir" => "./tasks", "outDir" => "./stdout", "errDir" => "./error", "modules" => \@modules,
		"nodes" => 1, "tasks" => 2, "mem" => "20000");
		
# Gather generic BAM stats
	# It's ugly, but unless I load another module, this is the best hack to determine the Lumpy script directory
	open(my $IN, "module load lumpy; which lumpy |");
	my $lumpy = <$IN>;
	my $lumpyDir = dirname($lumpy);
	chomp $lumpyDir;
	close $IN;
	my $lumpyScriptDir = "$lumpyDir/../scripts/";

print "Determining bam stats...\n";
my ($readLen, $mean, $sd) = determineStats($inputBam, $lumpyDir);
print "BAM stats: readlen:$readLen mean:$mean sd:$sd\n";

unless( -e "$outPrefix.lumpy.vcf"){
	print "Queuing Lumpy...\n";

	# Determine sample IDs
	my $sampleId;
	open(my $IN, "module load samtools; samtools view -H | grep '@RG' |");
	while(my $line = <$IN>){
		chomp $line;
		my @segs = split(/\t/, $line);
		if($segs[0] eq "@RG"){
			$line =~ /SM:(\S+)\s+/;
			$sampleId = $1;
			last;
		}
	}
	close $IN;
	
	my $PE = "-pe id:$sampleId,bam_file:$inputBam,histo_file:$inputBam.histo,mean:$mean,stdev:$sd,read_length:$readLen,min_non_overlap:$readLen,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20";
	
	my $cmd = "lumpy -mw 4 -tt 0 $PE > $outPrefix.lumpy.vcf";
	
		
	$sideWorker->createGenericCmd($cmd, "lumpyCmd.sh");
}

unless( -e "$outPrefix\_FRC.txt"){
	print "Queueing FRC...\n";
	my $min = $mean - ($sd * 3);
	if($min < 0){ $min = 0;}
	my $max = $mean + ($sd * 3);
	
	print "FRC min and max values: $min $max\n";
	my $PE = "-pe --pe-sam $inputBam --pe-max-insert $max";
	
	my $cmd = "FRC $PE --genome-size $opts{g} --output $outPrefix";
	
	$sideWorker->createGenericCmd($cmd, "frcCmd.sh");
}

unless( -e "$outPrefix.bayes.vcf"){
	print "Queueing Freebayes...\n";
	my @cmds;
	push(@cmds, "freebayes -C 2 -0 -O -q 20 -z 0.02 -E 0 -X -u -p 1 -F 0.5 -b $inputBam -v $outPrefix.bayes.vcf -f $opts{f}");
	push(@cmds, qq(NUM_SNP=`perl -e '$c = 0; while(<>){chomp; @F = split(/\t/); if($F[0] =~ /^#/){next;} ($ab) = $F[7] =~ /AB=(.{1,10})\;ABP/; if($ab < 0.65){next;}else{ $la = length($F[3]); $lb = length($F[4]); if($la == $lb){$c++;}elsif($la < $lb){$c += $lb - $la;}else{$c += $la - $lb;}}} print "$c\n";' < $inputBam`));
	push(@cmds, qq(NUM_BP=`samtools depth $inputBam | perl -e '$c = 0; while(<>){chomp; @s = split(/\t/); if($s[2] >= 3){$c++;}} print "$c\n";'`));
	push(@cmds, qq(QV=`perl -e 'chomp(@ARGV); $ns = $ARGV[0]; $nb = $ARGV[1]; print (-10 * log($ns/$nb)/log(10)); print "\n";' $NUM_SNP $NUM_BP`));
	push(@cmds, "echo $QV > $outPrefix.qv");
	
	$sideWorker->createArrayCmd(\@cmds, "freebayesCmd.sh");
}

$sideWorker->queueJobs;
	
exit;

sub determineStats{
	my ($inputBam, $lumpyScriptDir) = @_;
	
	open(my $IN, "module load samtools; samtools view -F 2048 $inputBam | head -n 1 |");
	my $samLine = <$IN>;
	close $IN;
	chomp $samLine;
	my @samSegs = split(/\t/, $samLine);
	my $readLen = length($samSegs[9]);
	
	open(my $IN, "module load lumpy; samtools view $inputBam | tail -n+100000 | $lumpyScriptDir/pairend_distro.py -r $readLen -X 4 -N 10000 -o $inputBam.histo |");
	my $pars = <$IN>;
	close $IN;
	my ($mean, $sd) = $pars =~ /mean:(.+)\s+stdev:(.+$)/;
	chomp $mean;
	chomp $sd;
	
	return $readLen, $mean, $sd;
}