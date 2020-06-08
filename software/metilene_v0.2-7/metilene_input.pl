#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use File::Spec;
use List::Util qw/sum min max/;
use File::Temp qw/tempfile tempdir/;
use File::Basename;

# -----------------------------------------------------------------------------
# Variables

my ($USAGE, $call, $ret, $help);
my ($in1, $in2, $out, $h1, $h2, $l, $bedtools);
my ($NA, $out_metilene, $g1_h, $g2_h, $g1_n, $g2_n);
my (@g1, @g2);
my $order = "";
my $SCRIPTNAME = basename($0);

# -----------------------------------------------------------------------------
# OPTIONS


## R executable! checking of ALL flags

$USAGE = << "USE";

    usage:  perl $SCRIPTNAME --in1 <list> --in2 <list> [--out <string>] [--h1 <string>] [--h2 <string>] [-b <path/prefix>]

    [INPUT]     --in1       comma-seperated list of sorted (!) bedgraph input files of group 1
                --in2       comma-seperated list of sorted (!) bedgraph input files of group 2
                --out       path/file of out file (metilene input) (default: metilene_g1_g2.input, g1 set by -h1 option, g2 set by -h2 option)
                --h1        identifier of group 1 (default: g1)
                --h2        identifier of group 2 (default: g2)
                -b          path/executable of bedtools executable (default: in PATH)
USE

if (!@ARGV) {
    printf STDERR $USAGE;
    exit -1;
}

unless (GetOptions(
    "in1=s" => \$in1,
    "in2=s" => \$in2,
    "out=s" => \$out,
    "h1=s"  => \$h1,
    "h2=s"  => \$h2,
    "b=s"   => \$bedtools,
    "h|help"=> \$help
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

my $cur_dir = File::Spec->curdir();

############
## checks ##
############
print STDERR ("[WARNING]" . prettyTime() . "Input files need to be SORTED, i.e., use \"bedtools sort -i file >file.sorted\"\n");
print STDERR ("[INFO]" . prettyTime() . "Checking flags\n");

if (defined $in1){
    @g1 = split(/,/,$in1);
    
    for (my $i=0; $i<scalar(@g1); $i++){
        $g1[$i] = File::Spec->rel2abs($g1[$i]);
        
        if ((defined $g1[$i]) && (-e $g1[$i])){
            unless (-r $g1[$i]){
                die "##### AN ERROR has occurred: $g1[$i] (option --g1) not readable\n";
            }
        }
        else{
            die "##### AN ERROR has occurred: required option --g1 ($g1[$i]) missing or nonexistent\n";
        }
    }
}
else {
    die "##### AN ERROR has occurred: required option --$in1 missing\n";
}

if (defined $in2){
    @g2 = split(/,/,$in2);
    
    for (my $i=0; $i<scalar(@g2); $i++){
        $g2[$i] = File::Spec->rel2abs($g2[$i]);
        
        if ((defined $g2[$i]) && (-e $g2[$i])){
            unless (-r $g2[$i]){
                die "##### AN ERROR has occurred: $g2[$i] (option --g2) not readable\n";
            }
        }
        else{
            die "##### AN ERROR has occurred: required option --g2 ($g2[$i]) missing or nonexistent\n";
        }
    }
}

$NA = "NA";

if (!defined $h1) {
    $h1 = "g1";
}
if (!defined $h2) {
    $h2 = "g2";
}
$g1_h = $h1;
$g2_h = $h2;
for (my $i = 1; $i < scalar(@g1); $i++){
    $g1_h = ($g1_h . " " . $h1);
}
for (my $i = 1; $i < scalar(@g2); $i++){
    $g2_h = ($g2_h . " " . $h2);
}

if (defined $out){
    $out            = File::Spec->rel2abs($out);
    $out_metilene   = $out;
}
else {
    $out_metilene   = ("metilene_" . $h1 . "_" . $h2 . ".input");
}

## bedtools ##
if (defined $bedtools){
    $bedtools = File::Spec->rel2abs($bedtools);
    if (-e $bedtools){
        unless (-d $bedtools){
            unless (-x $bedtools){
                die "##### AN ERROR has occurred: --bedtools option executable is not executable\n";
            }
        }
        else{
            die "##### AN ERROR has occurred: --bedtools option executable is directory\n";
        }
    }
    else{
        die "##### AN ERROR has occurred: --bedtools option executable nonexistent\n";
    }
}
else{
    $bedtools = "bedtools";
}
$call = "command -v $bedtools &> /dev/null";
$ret = system ($call);
if ($ret != 0){
    die "##### AN ERROR has occurred: No bedtools executable found. Please provide path/filename of bedtools executable with --bedtools option\n";
}

#####################
## union bedgraph ##
####################
print STDERR ("[INFO]" . prettyTime() . "Write metilene input to $out_metilene\n");

$g1_n = join(" ", @g1);
$g2_n = join(" ", @g2);

$call = "$bedtools unionbedg -header -names $g1_h $g2_h -filler $NA -i $g1_n $g2_n | cut -f1,3- | sed 's/end/pos/' >$out_metilene";
call($call);

print STDERR ("\n*****\n\n[BASIC CALL:] " . "metilene_executable -t 4 -a $h1 -b $h2 $out_metilene >out.file\nPlease adjust executable of metilene, number of threads (-t) and names of groups (-a, -b). For further parameters: metilene_executable --help\n");


# -----------------------------------------------------------------------------
# FUNCTIONS

sub call{
    my ($sub_call) = @_;
        
    $ret = system ($sub_call);
    
    if ($ret != 0){
        die "##### AN ERROR has occurred\n";
    }
}

sub prettyTime{
    my @months      = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
    my @weekDays    = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
    my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
    my $year        = 1900 + $yearOffset;
    return "\t$weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $hour:$minute:$second, $year\t";
}
