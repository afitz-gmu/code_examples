#!/usr/bin/perl -w
# File: project4.pl
# use strict;
# use warnings;
use CGI qw(:standard);  # module to output HTML commands

# set fatal error message to the browser window
use CGI::Carp qw/fatalsToBrowser/;

# variable declarations:
my @header = ();		    # array of headers
my @sequence = ();		    # array of sequences
my $count = 0;	           # number of sequences
my $lengthSequence = 0;    # used to hold the value of the length of each sequence before added to array
my @lengthArray = ();      #array that holds sequence lengths
my $total = 0;             # total number of sequences in file
my $average = 0;           # average length of sequence being analyzed
my $sub = 0;              # specific based being analyzed in for loop


# based in part on code from lecture 10 by Dr. Solka 
my $url = "/afitz/cgi-bin/project4.pl";
print header;
print start_html('Project 4'),
    h3('Project 4: Nucleotide Counter'),
    start_multipart_form, p,
    "Click the button to choose a FASTA file:",
    br, filefield(-name=>'filename'), p,
    reset, submit('submit','Submit File'), hr;


if (param()) {
    my $filehandle = upload('filename');
    # read whole file into @data array
    # my @data = <$filehandle>;
    
    # read FASTA file
    my $n = -1;			    # index of current sequence
    while (my $line = <$filehandle>) {
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
    close $filehandle;


    # remove white space from all sequences
    for (my $i = 0; $i < $count; $i++) {
        $sequence[$i] =~ s/\s//g;
        $sequence[$i] = uc $sequence[$i];
        $lengthSequence = length $sequence[$i];  # holds length of sequence
        push(@lengthArray, $lengthSequence);     # length added to array holding all lengths
        $total += $lengthSequence;               # cumulative sum of all lengths

    }

 ########## Sequence processing starts here:

print p, "File Sequence Statistics", p;
# Number of sequences
print "There are $count sequences in the file \n", br;
print "Total sequence length = $total\n", br;


# Sorting of array
my @nums = sort { $a <=> $b } @lengthArray;              # sorts the sequence lengths from smallest to largest
print "Maximum sequence length = $nums[$count-1]\n", br;     # prints the largest by pulling the last value in the array
print "Minimum sequence length = $nums[0]\n", br;            # prints the smallest by pulling the first value in the array


# Average calculation
$average = $total / $count;                      # calculates the average sequence length       
print("Ave sequence length = $average\n"), br, p;		 # prints that length

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
    print hr, p, "$header[$x]\n", p;
	
	# The length of the current sequence is calculated and printed
	$lengthSequence = length $sequence[$x];
	print "Length:$lengthSequence\n", br;
	
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
	printf "CpG:$d %0.2f \n", $percent, br;

}   
    
    print address( a({href=>$url},"Click here to submit another file."));

}

print end_html;
exit;



