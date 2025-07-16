#!/usr/bin/env perl

$inprot = @ARGV[0];
$intrans = @ARGV[1];
$outdir = @ARGV[2];

open(REPFASP, '>', "$outdir/representatives.prot.fasta");
open(REPFAST, '>', "$outdir/representatives.cds.fasta");
open(REPTAB, '>', "$outdir/representatives.tab");

%rephash = undef;
%cdshash = undef;

opendir(IN, $intrans);
my @files = readdir(IN);
closedir(IN);

foreach my $file (@files) { unless ($file eq undef || $file =~ m/^\.+/ || $file eq "singletons.fasta") {
	open(FAS, '<', "$intrans/$file");
	my $fasta = do { local $/; <FAS> };
	@fasarray = split(/\>/, $fasta);
	foreach $fasrec (@fasarray) { unless ($fasrec eq undef) {
		@temp = split(/\n/, $fasrec, 2);
		$id = @temp[0];
		$id =~ s/mrna\.//sgi;
		$seq = @temp[1];
		$seq =~ s/\n|\s//sgi;
		$cdshash{$id} = $seq;
	}}
	close FAS;
}}

print "printing cdshash keys";

foreach my $key (keys %cdshash) {
    print "$key\n";
}

opendir(IN, $inprot);
my @files = readdir(IN);
closedir(IN);


foreach my $file (@files) { unless ($file eq undef || $file =~ m/^\.+/ ) {
	open(FAS, '<', "$inprot/$file");
	%isolates = undef;
	%seqstats = undef;
	%fashash = undef;
	my $fasta = do { local $/; <FAS> };
	@fasarray = split(/\>/, $fasta);
	$numseqs = scalar(@fasarray);
	foreach $fasrec (@fasarray) { unless ($fasrec eq undef) {
		@temp = split(/\n/, $fasrec, 2);
		$id = @temp[0];
		$seq = @temp[1];
		$seq =~ s/\n|\s//sgi;
		##exclude sequences containing "X" unknown amino acids, or CQPM predictions containing "*" characters, which are santisied by converting to "X" prior to use as proteinortho input
		unless ($seq =~ /X|x/) {
			$fashash{$id} = $seq;
			$length = length($seq);
			$seqstats{"length_by_id"}->{$id} = $length;
			$seqstats{"id_by_length"}->{$length} = $id;
		}		
	}}
	close FAS;
	
	#identify representative
	$average = 0;
	$representative = undef;
	foreach $length (sort {$b <=> $a} keys %{$seqstats{"id_by_length"}}) {
		$representative = $seqstats{"id_by_length"}->{$length};
		last;
	}					

	@temp = split(/_/, $representative);
	$prefix= @temp[0];
	$groupno= @temp[1];
	$orthogroup = $prefix . "_" . $groupno;
	$isolate = @temp[2];
	$gene = @temp[-3] . "_" . @temp[-2] . "_" . @temp[-1];
	print "gene: $gene\n";
	print "representative: $representative\n";
	#output new multi fasta 
	open(OUTP, '>', "$outdir/$file");
	open(OUTT, '>', "$outdir/$file.cds");
	print OUTP ">$orthogroup\n" . $fashash{$representative} . "\n";
	print OUTT ">$orthogroup\n" . $cdshash{$gene} . "\n";	
	foreach $id (sort {$a cmp $b} keys %fashash) { unless ($id eq undef || $id eq $representative) {
		print OUTP ">$id\n" . $fashash{$id} . "\n";
		@temp = split(/_/, $id);
		$tx = @temp[-2] . "_" . @temp[-1];
		print OUTT ">$id\n" . $cdshash{$tx} . "\n";
	}}
	close OUTP;
	close OUTT;
	$length = length( $fashash{$representative} );
	$rephash{$orthogroup}->{'table'} = "$representative\t$orthogroup\t$isolate\t$gene\t" . $length . "\t" . $fashash{$representative} . "\n";
	$rephash{$orthogroup}->{'fastap'} = ">$orthogroup\n" . $fashash{$representative} . "\n";
	$rephash{$orthogroup}->{'fastat'} = ">$orthogroup\n" . $cdshash{$gene} . "\n";
	
}}

foreach $orthogroup (sort {$a cmp $b} keys %rephash) { unless ($orthogroup eq undef || $rephash{$orthogroup}->{'table'} =~ m/ortho/) {
	print REPTAB $rephash{$orthogroup}->{'table'};
	print REPFASP $rephash{$orthogroup}->{'fastap'};
	print REPFAST $rephash{$orthogroup}->{'fastat'};
}}

close REPFASP;
close REPFAST;
close REPTAB;