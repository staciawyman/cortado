#!/usr/bin/perl
#use strict;
#use warnings;
use Getopt::Long;


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
while (<>) {
    if (/^#/ || /Ref_name/) { next; }
    chomp;
    my (@cols) = split(/\t/);
    if (scalar(@cols) != 6)  { die("Manifest must have 6 columns"); }

    my ($samp,$ref,$refseq,$donorseq,$main_site,$guideseq) = split(/\t/);
    my $name = $samp."_".$ref;
    if ($ref =~ /\s/) { die("No spaces allowed in reference names."); }
    if (length($refseq) > 600) {
	$r = substr($refseq,0,300);
	$refseq = $r;
	$name .= "_long_target";
    }

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
        print "/usr/bin/python $cortado_path -r1 $fastq_path/$samp\*R1\*.fastq.gz  -o output -n $name -a $refseq -g $guideseq --trim_sequences  --keep_intermediate  --min_identity_score 58 --window_around_sgrna $window_size  --min_frequency_alleles_around_cut_to_plot 0.1 --max_rows_alleles_around_cut_to_plot 250";
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
