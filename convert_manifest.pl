#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;


my $help = 0;
my $threads = 24;
my $window_size = 6;   
my $cortado_path = `which cortado.py`;
chomp($cortado_path);

# Set default fastq fil path here, or give as command line argument
my $fastq_path = "/data/amplicon_fastqs";

GetOptions ("t=i" => \$threads,    # numeric
	    "h=s" => \$help,
	    "w=i" => \$window_size,
            "c=s"   => \$cortado_path,
            "f=s"    => \$fastq_path,
            "h"    => \$help,
            "help"    => \$help,
)    or die("USAGE: perl convert_manifest.pl [-w number] [-t number_of_threads] [-f /path/to/fastqs]  <manifest text file>\n");

if ($help || !@ARGV) {
    print "USAGE: perl convert_manifest.pl [options] <manifest text file>\n";
    print "\tOPTIONS:\n";
    print "\t-h\tThis message\n";
    print "\t-t int\tNumber of threads to use (runs in batches)\n";
    print "\t-w int\tWindow size around cut site within which to check for indels\n";
    print "\t-f /path/to/fastq/dir\n\t-c /path/to/cortado.py \n";
    exit;
}

my $count = 0;
while (<>) {
    if (/^#/ || /Ref_name/) { next; }

    chomp;
    my ($samp,$ref,$refseq,$donorseq,$main_site,$guideseq) = split(/\t/);
    my $name = $samp."_".$ref;

    if (length($refseq) != length($donorseq)) { die("Sample $name: reference and donor sequences must be same length."); }
    $guideseq = uc($guideseq);
    $refseq = uc($refseq);
    if ($refseq !~ /$guideseq/) {
	my $rc = revcomp($refseq);
	if ($rc !~ /$guideseq/) {
	    die ("Sample $name: reference sequence must contain guide sequnce.");
	}
    }
    
    $samp .= "_";
    if ($donorseq eq "-") { # Run as just NHEJ
        print "/usr/bin/python $cortado_path -r1 $fastq_path/$samp\*R1_001.fastq.gz -r2 $fastq_path/$samp\*R2_001.fastq.gz -o output -n $name -a $refseq -g $guideseq --trim_sequences --trimmomatic_options_string ILLUMINACLIP:adapters.fa:2:30:10 --keep_intermediate  --min_identity_score 58 --window_around_sgrna $window_size  --min_frequency_alleles_around_cut_to_plot 0.1 --max_rows_alleles_around_cut_to_plot 250";
	if ($count < $threads) {
		print " & \n";
		$count++
	} else {
		print "  \n";
		$count = 0;
	}
    } else { # Run as HDR
        print "/usr/bin/python $cortado_path -r1 $fastq_path/$samp\*R1_001.fastq.gz -r2 $fastq_path/$samp\*R2_001.fastq.gz -o output -n $name -a $refseq -e $donorseq -g $guideseq --trim_sequences --trimmomatic_options_string ILLUMINACLIP:adapters.fa:2:30:10 --keep_intermediate --all_edits --main_site $main_site --min_identity_score 58 --window_around_sgrna $window_size --min_frequency_alleles_around_cut_to_plot 0.1 --max_rows_alleles_around_cut_to_plot 250 ";
	if ($count < $threads) {
		print " & \n";
		$count++
	} else {
		print "  \n";
		$count = 0;
	}
    }
}

sub revcomp {
    my ($dna) = @_;
    my $rc = reverse($dna);
    $rc =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
    return $rc;
}
