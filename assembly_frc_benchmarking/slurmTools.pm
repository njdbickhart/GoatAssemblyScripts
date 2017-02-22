#!/usr/bin/perl
# Tools to initiate Slurm shell script generation and execution

package slurmTools;
use Mouse;
use namespace::autoclean;

has ['workDir', 'scriptDir', 'outDir', 'errDir'] => (is => 'ro', isa => 'Str', required => 1);
has 'modules' => (is => 'rw', isa => 'ArrayRef[Any]', handle => 'has_module');
has ['nodes', 'tasks', 'mem', 'time'] => (is => 'rw', isa => 'Any', default => -1);
has 'jobIds' => (is => 'rw', isa => 'ArrayRef[Any]', handle => 'has_jobs');

sub checkJobs{
	my ($self) = @_;
	
	my @jobIds = @{$self->jobIds};
	my @incomplete;
	foreach my $j (@jobIds){
		my $output = `squeue -j $j`;
		if($output =~ /slurm_load_jobs error:/){
			push(@incomplete, $j);
		}
	}
	$self->jobIds(@incomplete);
	if(scalar(@incomplete) > 0){
		return 0;
	}else{
		return 1;
	}
}
			

sub queueJobs{
	my ($self) = @_;
	
	my @scripts = `ls $self->scriptDir`;
	my @jobIds;
	foreach my $s (@scripts){
		my $jid = `sbatch $s`;
		push(@jobIds, $jid);
	}
	$self->jobIds(\@jobIds);
}

sub createArrayCmd{
	my ($Self, $carrayref, $sbase) = @_;
	# Requires an array ref of premade cmds for a single script
	$self->_generateFolders;
	if(!defined($sbase)){
		$sbase = "script_";
	}
	
	my $head = $self->_generateHeader($sname);
	
	foreach my $cmd (@{$carrayref}){
		$head .= "echo $cmd\ntime $cmd\n\n";
	}
	$head .= "wait\n";
	
	my $sFolder = $self->scriptDir;
	open(my $OUT, "> $sFolder/$sname") || die "Could not create script!\n";
	print {$OUT} $head;
	close $OUT;
}

sub createGenericCmd{
	my ($self, $cmd, $sname) = @_;
	# Requires detailed command statement and [optionally] a script name
	$self->_generateFolders;
	if(!defined($sname)){
		my $hash = $self->_generateSHash($cmd);
		$sname = "script_$hash.sh";
	}
	
	my $head = $self->_generateHeader($sname);
	$head .= "echo $cmd\ntime $cmd\nwait\n";
	
	my $sFolder = $self->scriptDir;
	open(my $OUT, "> $sFolder/$sname") || die "Could not create script!\n";
	print {$OUT} $head;
	close $OUT;
}

sub _generateFolders{
	my ($self) = @_;
	mkdir $self->workDir || print "$!\n";
	mkdir $self->scriptDir || print "$!\n";
	mkdir $self->outDir || print "$!\n";
	mkdir $self->errDir || print "$!\n";
}

sub _generateSHash{
	my ($self, $cmd) = @_;
	# Generates unique hash from command name
	
	my $h = 0;
	
	foreach my $c (unpack("C*", $cmd)){
		my $high = $c & 0xf8000000;
		$h = $h << 5;
		$h = $h ^ ($high >> 27);
		$h = $h ^ $c;
	}
	
	return $h;
}
	

sub _generateHeader{
	my ($self, $sname) = @_;
	my $meta = __PACKAGE__->meta;
	my $tag = "#SBATCH";
	my $str = "#!/bin/bash\n";
	if($self->nodes != -1){
		$str .= "$tag --nodes=" . $self->nodes . "\n";
	}
	if($self->tasks != -1){
		$str .= "$tag --ntasks-per-node=" . $self->tasks . "\n";
	}
	if($self->mem != -1){
		$str .= "$tag --mem=" . $self->mem . "\n";
	}
	if($self->time != -1){
		$str .= "$tag --time=" . $self->time . "\n";
	}
	$str .= "$tag --output=" . $self->outDir . "/$sname\_%j.out\n$tag --error=" . $self->errDir . "/$sname\_%j.err\n$tag --workdir=" . $self->workDir . "\n\n";
	
	$str .= "cd " . $self->workDir . "\n\n";
	
	if($self->has_module){
		foreach my $m (@{$self->modules}){
			$str .= "module load $m\n";
		}
		$str .= "\n";
	}
		
	return $str;
}


__PACKAGE__->meta->make_immutable;

1;