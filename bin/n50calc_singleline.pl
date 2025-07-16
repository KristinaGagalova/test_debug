#!/usr/bin/env perl
#usage perl script.pl name_of_fasta_file threshold
#NB includes * character when counting


open(FASTA,'<', @ARGV[0]);

if (@ARGV[1] ne undef) {
    $threshold = @ARGV[1];
}
else {
    $threshold = 0;
}

$fastadata = do {local $/; <FASTA>};
@fastarecords = split(/\>/, $fastadata);


$polyn = 0;
%sizes = undef;
$totallength = 0;
$seqcount = 0;
foreach $record (@fastarecords) {
	unless ($record eq undef) {
		@data = split(/\n/, $record, 2);
		@data[1] =~ s/\n|\s//sgi;
		$polyn += (@data[1] =~ tr/Nn/Nn/);
		$len = length(@data[1]);
		if ($len >= $threshold) {
			$sizes{@data[0]} = length(@data[1]);
			$totallength = $totallength + length(@data[1]);
			$seqcount++;
		}
	}
}
$max = 0;
$min = undef;

$runlen = 0;
$seqno = 0;
$check = 0;

foreach $seq (sort {$sizes{$b} <=> $sizes{$a}} keys %sizes) {
	unless ($seq eq undef) {
		$seqno++;
		if ($sizes{$seq} > $max) {
			$max = $sizes{$seq};
		}
		if ($min == undef) {
			$min = $sizes{$seq};
		}
		else {
			if ($sizes{$seq} < $min && $sizes{$seq} != 0) {
				$min = $sizes{$seq};
			}
		}
		$runlen = $runlen + $sizes{$seq};
		if ($runlen > ($totallength / 2) && $check == 0) {
			$n50 = $seqno;
			$l50 = $sizes{$seq};
			$check = 1;
		}
	}
}
##$seqno = $seqno - 1;
$temp = ($totallength / $seqcount);
print @ARGV[0]."\t" . $threshold . "\t" . $totallength . "\t" . $temp . "\t$max\t$min\t$n50\t$l50\t$seqno\t$polyn\n";
#print @ARGV[0]."\t" . $threshold . "\t" . $totallength . "\t$max\t$min\t$n50\t$l50\t$seqno\t$polyn\n";
exit;
