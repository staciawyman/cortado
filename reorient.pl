#!/usr/bin/perl

# given a reference sequence and input SE fastq file, orient the reads all in the same direction
#$ref = "TAGGGTTGGCCAATCTACTCCCAGGAGCAGGGAGGGCAGGAGCCAGGGCTGGGCATAAAAGTCAGGGCAGAGCCATCTATTGCTTACATTTGCTTCTGACACAACTGTGTTCACTAGCAACCTCAAACAGACACCATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACTCTCTCTGCCTATTGGTCTATTTTCCCA";

open ($fastq,  " gunzip -c  $ARGV[0] | ") or die "error opening fastq $ARGV[0]: $!";
$refstr = $ARGV[1];

#open my $outfile, '>:gzip', "output"  or die "Could not write to outut: $!";
#open my $outfile, '>:gzip', $ARGV[2] or die "Could not write to $ARGV[2]: $!";

open (my $gzip_outfile, "| gzip -c > $ARGV[2]") or die "error starting gzip $!";
$start = 15;
$end = 25;
$refstr_f = substr($refstr,$start,$end);
$refstr_r = substr(revcomp($refstr),$start,$end);
$skip = 0;
$rc = 0;

#print "Search for F: $refstr_f, $refstr_r\n";

while (<$fastq>) {
    chomp;
    if ($. % 4 == 1) { # Line 1, description
        $l1 = "$_";
    }
    if ($. % 4 == 2) { # Line 2, DNA sequence
        if (/$refstr_f/) { $l2 = "$_"; }
        elsif (/$refstr_r/) { $l2 = revcomp($_); $rc = 1;}
        else { $skip = 1; }
    }
    if ($. % 4 == 3) { # Line 3, +
	$l3 = "$_";
    }
    if ($. % 4 == 0) {
        if (!$skip) {
            # print first three lines
            print $gzip_outfile "$l1\n";
	    print $gzip_outfile "$l2\n";
	    print $gzip_outfile "$l3\n";

            if ($rc) {
	        $rc = 0;
                # Quality scores have to be reversed, but not reverse complemented
                print $gzip_outfile reverse($_),"\n";
            } else {
                print $gzip_outfile "$_\n";
            }
        }
	$rc = 0;
    	$l1 = "";
    	$l2 = "";
     	$l3 = "";
        $skip = 0;
    }
}

sub revcomp {
    my ($dna) = @_;
    my $rc = reverse($dna);
    $rc =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
    return $rc;
}
