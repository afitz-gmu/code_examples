#!/usr/bin/perl
use strict;
use warnings;

# File: fasta.pl
# Author: Alex Fitz
# Date: 02 September 2020
#
# Purpose: Read sequences from a FASTA format file
# Based in part on fasta.pl by Jeff Solka

# the argument list should contain the file name
#die "usage: fasta.pl filename\n" if scalar @ARGV < 1;

# get the filename from the argument list

#my ($filename) = @ARGV;
my ($filename) = "test1.fsa";

# Open the file given as the first argument on the command line
open(INFILE, $filename) or die "Can't open $filename\n";

print "Report for file $filename\n";

# variable declarations:
my @header = ();		    # array of headers
my @sequence = ();		    # array of sequences
my $count = 0;	           # number of sequences
my $lengthSequence = 0;    # used to hold the value of the length of each sequence before added to array
my @lengthArray = ();      #array that holds sequence lengths
my $total = 0;             # total number of sequences in file
my $average = 0;           # average length of sequence being analyzed
my $sub = 0;              # specific based being analyzed in for loop


# read FASTA file
my $n = -1;			    # index of current sequence
while (my $line = <INFILE>) {
    chomp $line;		    # remove training \n from line
    if ($line =~ /^>/) { 	    # line starts with a ">"
	$n++;			    # this starts a new header
	$header[$n] = $line;	    # save header line
	$sequence[$n] = "";	    # start a new (empty) sequence
    }
    else {
	next if not @header;	    # ignore data before first header
	$sequence[$n] .= $line     # append to end of current sequence
    }
}
$count = $n+1;			  # set count to the number of sequences
close INFILE;



# remove white space from all sequences
for (my $i = 0; $i < $count; $i++) {
    $sequence[$i] =~ s/\s//g;
	$sequence[$i] = uc $sequence[$i];
	$lengthSequence = length $sequence[$i];  # holds length of sequence
	push(@lengthArray, $lengthSequence);     # length added to array holding all lengths
    $total += $lengthSequence;               # cumulative sum of all lengths
}

########## Sequence processing starts here:

# Number of sequences
print "There are $count sequences in the file \n";
print "Total sequence length = $total\n";


# Sorting of array
my @nums = sort { $a <=> $b } @lengthArray;              # sorts the sequence lengths from smallest to largest
print "Maximum sequence length = $nums[$count-1]\n";     # prints the largest by pulling the last value in the array
print "Minimum sequence length = $nums[0]\n";            # prints the smallest by pulling the first value in the array


# Average calculation
$average = $total / $count;                      # calculates the average sequence length       
print("Ave sequence length = $average\n");		 # prints that length

my $currentSequence = 0;

# Individual Sequence Processing
# The for loop iterates through all of the sequences present in the fasta file
for (my $x = 0; $x < $count; $x++) {
	
	# I initialize the count and percent variables at a value of 0 so they clear after each iteration
	my $countA = 0;
	my $countT = 0;
	my $countG = 0;
	my $countC = 0;
	my $cpg = 0;
	my $percentA = 0;
	my $percentT = 0;
	my $percentG = 0;
	my $percentC = 0;

	# The program prints out the header of the sequence being analyzed
    print "$header[$x]\n";
	
	# The length of the current sequence is calculated and printed
	$lengthSequence = length $sequence[$x];
	print "Length:$lengthSequence\n";
	
	# The current sequence is saved in a variable so it can be used in the following for loop for calculations
	$currentSequence = $sequence[$x];
	
	
	my @cpg = ($currentSequence =~ /CG/g);    # Based on the current sequence the number of CpGs are calculated and proportion is calculated
	my $d = @cpg;							  # saved as a scalar variable to print out later
	my $percent = @cpg / $lengthSequence;     # proportion of CpG calculated
	
	# This for loop iterations through every base in the sequence counting how many of each base exist and adding them up in their individual variables
	for (my $b = 0; $b < $lengthSequence; $b++){
		$sub = substr($currentSequence, $b, 1);     # varaible holds the value of the individual sequence being analyzed
		 
		if ($sub eq 'G'){
			$countG += 1;
		}
		
		elsif ($sub eq 'C'){
			$countC += 1;
		}
		
		elsif ($sub eq 'A'){
			$countA += 1;
		}		
		
		elsif ($sub eq 'T'){
			$countT += 1;
		}
		
		# The proportion of each base is calculated and saved for G / C / A / T
		
		$percentG = $countG / $lengthSequence;
		$percentC = $countC / $lengthSequence;
		$percentA = $countA / $lengthSequence;
		$percentT = $countT / $lengthSequence;
	}
	
	#Using the correct formatting printf command the number of each base and the proportion are printed sequentially
	printf "A:$countA %0.2f\n", $percentA;
	printf "C:$countC %0.2f\n", $percentC;
	printf "G:$countG %0.2f\n", $percentG;
	printf "T:$countT %0.2f\n", $percentT;
	printf "CpG:$d %0.2f \n", $percent;
}

exit;

