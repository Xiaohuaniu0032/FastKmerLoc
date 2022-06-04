use strict;
use warnings;
use FindBin qw/$Bin/;

my ($fq,$db,$bio_name,$outdir) = @ARGV;

# load db
my $kmer = "$Bin/db/2019nCOV/uniqness.25.txt";
my %kmer;
open KMER, "$kmer" or die;
while (<KMER>){
	chomp;
	my @arr = split /\t/, $_;
	if ($arr[-1] eq "UNIQ"){
		my $seq = $arr[2];
		my $pos = $arr[1];
		my $genome = $arr[0];
		$kmer{$seq} = "$genome\t$pos"; # 
	}
}
close KMER;


# read fq
my $kmer_len = 25;
my %hit_info;

my $line = 0;
open FQ, "$fq" or die;
while (<FQ>){
	chomp;
	my $seq = <FQ>;
	chomp $seq;
	
	$line += 1;
	if ($line % 100 == 0){
		# 100/200/300 reads
		print "had processed $line reads...\n";
	}

	<FQ>;
	<FQ>;
	my $len = length($seq);
	my $last_pos = $len - $kmer_len - 1; # 0-based
	for my $i (0..$last_pos){
		my $subseq = substr($seq,$i,$kmer_len);
		my $rc_subseq = reverse($subseq);
		$rc_subseq =~ tr/ATCG/TAGC/;
		if (exists $kmer{$subseq}){
			$hit_info{$subseq} += 1;
		}elsif (exists $kmer{$rc_subseq}){
			$hit_info{$rc_subseq} += 1;
		}else{
			next;
		}
	}
}
close FQ;


my %genome_pos_kmer_info;
foreach my $seq (keys %hit_info){
	my $count = $hit_info{$seq};
	#print "$seq\t$count\n";
	my $seq_info = $kmer{$seq}; # genome/pos
	#print "$seq_info\n";
	my @seq_info = split /\t/, $seq_info;
	my $g = $seq_info[0];
	my $pos = $seq_info[1];
	$genome_pos_kmer_info{$g}{$pos} = "$seq\t$count";
}

foreach my $g (keys %genome_pos_kmer_info){
	my @pos = sort {$a <=> $b} keys %{$genome_pos_kmer_info{$g}};
	for my $pos (@pos){
		my $v = $genome_pos_kmer_info{$g}{$pos};
		print "$g\t$pos\t$v\n"; # genome/pos/seq/count
	}
}