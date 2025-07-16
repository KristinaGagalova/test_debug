#!/usr/bin/env perl

$indir = @ARGV[0];
$outdir = @ARGV[1];

$keysp = "NOT_APPLICABLE";
$prefix_old = "NOT_APPLICABLE";
$prefix_new = @ARGV[2];
$novel_num = 1;
$digits = 6;


##usage cat proteinortho6tsv | script.pl inputdir outputdir PREFIX
#inputdir: directory with protein *.fasta files with same names used in proteinortho column headers
#outputdir: new directory where orthogroup and singleton *.fasta will be written
##note: this script relies on inputdir *.fasta to have sequences on a single line to enable use of grep -A1 to extract sequences efficiently

system("rm -Rf $outdir");
system("mkdir -p $outdir");

open(TAB, '>', $outdir . "/ortho_table_revised.tab");

while(<STDIN>) {
	$line = $_;
	chomp $line;
	@data = split(/\t/, $line);
	if ($line =~ m/\# Species/) {
		@headers = split(/\t/, $line);
		print TAB "orthogroupID\t$line\n";
	}
	else {
		@data = split(/\t/, $line);
		# identify key species member and if so assign orthogroup
		$orthogroup = undef;
		$species = undef;
		for ($i=3; $i <= @data; $i++) {
			$species = @headers[$i];
			@seqs = split(/\,/, @data[$i]);
			if ($species eq $keysp) {
				%isohash = undef;
				foreach $seq (sort {$a <=> $b} @seqs) { unless ($seq eq undef) {
					$seq =~ s/$prefix_old//ig;
					$temp = substr($seq, 0, -1);
					$isohash{$temp} .= substr($seq, -1, 1);
				}}
				foreach $seq (sort {$a <=> $b} keys %isohash) { unless ($seq eq undef) {
					$orthogroup .= $seq . $isohash{$seq};
				}}
			}
		}
		if ($orthogroup eq undef) {
			my $padded_number = sprintf("%0${$digits}d", $novel_num);
			$orthogroup = $prefix_new . "_" . $padded_number;
			$novel_num++;
		}
		else {
			$orthogroup = $prefix_new . $orthogroup; 
		}
		print TAB "$orthogroup\t$line\n";
		
		#go through all orthologs within the orthogroup line and output into multiple (orthogroup.fasta)s or dump singletons in one file (singletons.fasta)
		for ($i=3; $i <= @data; $i++) {
			if (@data[1] > 1) {
				open(OUT, '>>', $outdir . "/" . $orthogroup . ".fasta");
			}
			else {
				open(OUT, '>>', $outdir . "/singletons.fasta");
			}			
			$species = @headers[$i];
			@seqs = split(/\,/, @data[$i]);
			foreach $seq (@seqs) { unless ($seq eq undef || $seq eq '*') {
				$temp = `grep -A1 $seq < $indir/$species`;
				$temp2 = $species;
				$temp2 =~ s/\.fasta//;
				$temp2 = $orthogroup ."_". $temp2 ."_";
				$temp =~ s/\>/\>$temp2/;
				print OUT $temp;
			}}
			close OUT;
		}

	}
}
close TAB;
