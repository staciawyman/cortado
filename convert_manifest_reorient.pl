#!/usr/bin/perl
#use strict;
#use warnings;
use Getopt::Long;

# This script converts an analysis manifest into cortado commands using the version of cortado 
# that reorients the reads to all be in the same direction matching the reference sequence.

my $help = 0;
my $threads = 30;
my $window_size = 6;   
my $cortado_path = `which cortado_one_dir.py`;
chomp($cortado_path);

GetOptions ("t=i" => \$threads,    # numeric
	    "h=s" => \$help,
	    "w=i" => \$window_size,
            "c=s" => \$cortado_path,
            "f=s" => \$fastq_path,
            "h"   => \$help,
            "help"  => \$help,
)    or die("USAGE: perl convert_manifest.pl [-w number] [-t number_of_threads] [-f /path/to/fastqs]  <manifest text file>\n");

if (!$fastq_path) {
    die("-f option is required. Please supply the path to the fastq directory");
}

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
$orig_window_size = $window_size;
while (<>) {
    if (/^#/ || /Ref_name/) { next; }
    chomp;
    my (@cols) = split(/\t/);
    if (scalar(@cols) < 6 || scalar(@cols) > 7)  { die("Manifest must have 6 or 7 columns"); }
    if ($cols[6] > 0) {
        ($samp,$ref,$refseq,$donorseq,$main_site,$guideseq,$window) = split(/\t/);
	$window_size = $window;
    } else {
        ($samp,$ref,$refseq,$donorseq,$main_site,$guideseq) = split(/\t/);
	$window_size = $orig_window_size;
    }
    my $name = $samp."_".$ref;
    if ($ref =~ /\s/) { die("No spaces allowed in reference names."); }
    $guideseq = uc($guideseq);
    $refseq = uc($refseq);
    if ($refseq !~ /$guideseq/) {
	my $rc = revcomp($refseq);
	if ($rc !~ /$guideseq/) {
	    #die ("ERROR: Sample $name: reference sequence must contain guide sequence.\n");
	    print ("ERROR: Sample $name: reference sequence must contain guide sequence.\n");
            next;
	}
    }
    
    $samp .= "_";
    if ($donorseq eq "-") { # Run as just NHEJ
        print "/usr/bin/python $cortado_path -r1 $fastq_path/$samp\*R1\*.fastq.gz -r2 $fastq_path/$samp\*R2\*.fastq.gz -o output -n $name -a $refseq -g $guideseq --trim_sequences  --keep_intermediate  --min_identity_score 58 --window_around_sgrna $window_size  --min_frequency_alleles_around_cut_to_plot 0.1 --max_rows_alleles_around_cut_to_plot 250";
	if ($count < $threads) {
		print " & \n";
		$count++
	} else {
		print "  \n";
		$count = 0;
	}
    } else { # Run as HDR
        if (length($refseq) != length($donorseq)) { 
		#die("Sample $name: reference and donor sequences must be same length.\n$refseq\n$donorseq\n"); 
		print("ERROR:Sample $name: reference and donor sequences must be same length.\nERROR:$refseq\nERROR:$donorseq\n"); next;
	}
        print "/usr/bin/python $cortado_path -r1 $fastq_path/$samp\*R1\*.fastq.gz -r2 $fastq_path/$samp\*R2\*.fastq.gz -o output -n $name -a $refseq -e $donorseq -g $guideseq --trim_sequences --keep_intermediate --all_edits --main_site $main_site --min_identity_score 58 --window_around_sgrna $window_size --min_frequency_alleles_around_cut_to_plot 0.1 --max_rows_alleles_around_cut_to_plot 250 ";
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
