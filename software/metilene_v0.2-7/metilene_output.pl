#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use File::Spec;
use File::Path qw(make_path remove_tree);
use File::Basename;

# -----------------------------------------------------------------------------
# Variables

my ($USAGE, $help, $ret, $call, $in, $out, $path);
my ($threads, $min_cpgs, $min_diff, $id1, $id2, $p_val, $length);
my ($plot_script, $pval_out, $pval_bed, $pdf);
my $SCRIPTNAME = basename($0);

# -----------------------------------------------------------------------------
# OPTIONS

$USAGE = << "USE";

    usage:  perl $SCRIPTNAME  -q <query_file> [-o <path_prefix>] [-p  <number>] [-c <number>] [-d <number>] [-l <number>] [-a <string>] [-b <string>]
    
    [INPUT]     -q          path/filename of metilene DMRs
                -o          path/prefix of output files (default: input_path/)
                -p          maximum (<) adj. p-value (q-value) for output of significant DMRs (default: 0.05)
                -c          minimum (>=) cpgs (default:10)
                -d          minimum mean methylation difference (>=) (default:0.1)
                -l          minimum length of DMR [nt] (>=) (post-processing, default: 0)
                -a          name of group A (default:"g1")
                -b          name of group B (default:"g2")
USE

if (!@ARGV) {
    printf STDERR $USAGE;
    exit -1;
}

unless (GetOptions(
    "q=s"       => \$in,
    "o=s"       => \$out,
    "p=s"       => \$p_val,
    "c=s"       => \$min_cpgs,
    "d=s"       => \$min_diff,
    "l=s"       => \$length,    
    "a=s"       => \$id1,
    "b=s"       => \$id2,
    "h|help"    => \$help
)){
    printf STDERR $USAGE;
    exit -1;
}
if (defined $help){
    printf STDERR $USAGE;
    exit -1;
}

# -----------------------------------------------------------------------------
# MAIN

#################
## check flags ##
#################
print STDERR ("[INFO]" . prettyTime() . "Checking flags\n");

## relative to absolut path ##
if (defined $in){$in = File::Spec->rel2abs($in)};
if (defined $out){$out = File::Spec->rel2abs($out)};
$plot_script = File::Spec->rel2abs($0);

## out dir
my ($vol, $dir, $file) = File::Spec->splitpath($in);
if (!defined $out){
    $path = $dir;
    $file = "metilene";
}
else {
    ($vol, $dir, $file) = File::Spec->splitpath($out);
    $path = $dir;
}
unless (-e $path){
    $ret = make_path($path);
    if ($ret != 1){
        die "##### AN ERROR has occurred: Could not make out directory $path\n";
    }
    else{
        print STDERR  ("[INFO]" . prettyTime() . "Created result-directory $path\n");
    }
}
elsif (-d $path){
    unless (-w $path && -r $path){
        die "##### AN ERROR has occurred: out directory $path not readable or writable\n";
    }
}

## id group 1 ##
if (!defined $id1){
    $id1 = "g1";
}
    
## id group 2 ##
if (!defined $id2){
    $id2 = "g2";
}

## min cpgs ##
if (!defined $min_cpgs){
    $min_cpgs = 10;
}

## min diff ##
if (!defined $min_diff){
    $min_diff = 0.1;
}

## pval filter ##
if (!defined $p_val){
    $p_val = 0.05;
}

## length ##
if (!defined $length){
    $length = 0;
}

## out file ##
$pval_out = File::Spec->catpath($vol, $path, $file);
$pval_out = ($pval_out . "_qval." . $p_val . ".out");
$pval_bed = File::Spec->catpath($vol, $path, $file);
$pval_bed = ($pval_bed . "_qval." . $p_val . ".bedgraph");
$pdf      = File::Spec->catpath($vol, $path, $file);
$pdf      = ($pdf . "_qval." . $p_val . ".pdf");

$plot_script =~ s/.pl$/.R/;

##################
### filter DMRs ##
##################
print STDERR ("[INFO]" . prettyTime() . "Filter DMRs.\n");

my $count = 0;
open(OUT, ">$pval_out") or die "##### AN ERROR has occurred: could not write to $pval_out\n";
open(OUT_B, ">$pval_bed") or die "##### AN ERROR has occurred: could not write to $pval_out\n";

open (IN, "$in") || die "cannot open $in\n";
    while (<IN>) {
    chomp;
    
    my ($chr, $start, $end, $pval, $diff_methyl, $nr_cpgs, $MWU, $KS2D, $mean_g1, $mean_g2) = split(/\t/,$_);
    if (($pval < $p_val) && (abs($diff_methyl) >= $min_diff) && ($nr_cpgs >= $min_cpgs) && (($end-$start) >= $length)) {
        print OUT "$chr\t$start\t$end\t$pval\t$diff_methyl\t$nr_cpgs\t$mean_g1\t$mean_g2\n";
        print OUT_B "$chr\t$start\t$end\t$diff_methyl\n";
        $count++;
    }
}
close (IN);
close (OUT);
close (OUT_B);

print STDERR ("[INFO]" . prettyTime() . "Wrote $count DMRs with adj. p-value<$p_val, a minimum absolute difference>=$min_diff, a minimum length [CpG]>=$min_cpgs and a minimum length [nt]>=$length\n");
print STDERR ("[INFO]" . prettyTime() . "Bedgraph file containing DMR difference: $pval_bed\n");


#################
### plot stats ##
#################
#print "$plot_script $pval_out $pdf\n";
print STDERR ("[INFO]" . prettyTime() . "Plot DMR statistics.\n");
$call = ("Rscript $plot_script $pval_out $pdf");
call($call);

# -----------------------------------------------------------------------------ls
# FUNCTIONS

sub call{
    my ($sub_call) = @_;
    $ret = system ($sub_call);

    if ($ret != 0){
        die "##### AN ERROR has occurred. Occurred at: $sub_call\n";
    }
}

sub prettyTime{
    my @months      = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my @weekDays    = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    my $year        = 1900 + $yearOffset;
    return "\t$weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $hour:$minute:$second, $year\t";
}



