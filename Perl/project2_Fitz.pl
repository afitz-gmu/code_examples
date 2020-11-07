#!/usr/bin/perl
use strict;
use warnings;
use List::Util 'shuffle';

# File: program2_Fitz.pl 
# Author: Alex Fitz
# Date: 09 October 2020
#
# Purpose: Read sequences from a FASTA format file
# Based in part on fasta.pl by Jeff Solka
# shuffle list used to randomize a numerical sequence in to a randomlized numerical sequence

# variable declarations:
my @header = ();		    # array of headers
my @sequence = ();		    # array of sequences
my $count = 0;	           # number of sequences
my $average = 0;           # average length of the sequences
my $lengthSequence = 0;    # used to hold the value of the length of each sequence before added to array
my @lengthArray = ();      #array that holds sequence lengths
my $total = 0;             # total number of sequences in file
my $sub = 0;              # specific based being analyzed in for loop
my $currentSequence = 0; # holds the current sequence in a scalar variable
my @currentSequence = (); # holds the current sequence being analyzed
my @newRandomSequences = (); # holds the random sequence
my $permuteFile = (); # new permute file name
my @acounts = (); # used to hold the number of a counts
my @ccounts = (); # used to hold the number of c counts
my @gcounts = (); # used to hold the number of g counts
my @tcounts = (); # used to hold the number of t counts
my @cgcounts = (); # used to hold the number of cg counts
my @aprops = (); # used to hold the proportion of a's
my @cprops = (); # used to hold the proportion of c's
my @tprops = (); # used to hold the proportion of t's
my @gprops = (); # used to hold the proportion of g's
my @cgprops = (); # used to hold the proportion of cg's


my ($filename) = @ARGV;    

sub read_fasta {

    # Open the file given as the first argument on the command line
    open(INFILE,  $filename) or die "Can't open $filename\n";

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

    return @header, @sequence;

}

sub stat_fasta {

    my ($acounts, $ccounts, $gcounts, $tcounts, $cgcounts, $aprops, $tprops, $gprops, $cprops, $cgprops, $file) = @_;
    
    #creating array so I can extract the filename without the suffix to create .ot file
    my @newfilename = split('\.', $file, -1);

    # creating new output file name with .ot suffix
    my $outputname = $newfilename[0] . ".ot";

    # write the following to the new filename .ot
    open(NewOT, '>', $outputname) or die $!;

    # prints the filename in the report
    print NewOT "Report for file $file\n";

    # Number of sequences
    print NewOT "There are $count sequences in the file \n";
    print NewOT "Total sequence length = $total\n";

    # Sorting of array
    my @nums = sort { $a <=> $b } @lengthArray;              # sorts the sequence lengths from smallest to largest
    print NewOT "Maximum sequence length = $nums[$count-1]\n";     # prints the largest by pulling the last value in the array
    print NewOT "Minimum sequence length = $nums[0]\n";            # prints the smallest by pulling the first value in the array


    # Average calculation
    $average = $total / $count;                      # calculates the average sequence length       
    print NewOT "Ave sequence length = $average\n";		 # prints that length

    # Individual Sequence Processing
    # The for loop iterates through all of the sequences present in the fasta file
    for (my $x = 0; $x < $count; $x++) {
        
        my($gprops, $cprops, $aprops, $tprops, $cgprops, $acounts, $tcounts, $gcounts, $ccounts, $cgcounts) = @_;

        # creating arrays to hold counts and proportions
        my $ac = 0;
        my $cc = 0;
        my $gc = 0;
        my $tc = 0;
        my $cgc = 0;
        my $ap = 0;
        my $cp = 0;
        my $tp = 0;
        my $gp = 0;
        my $cgp = 0;

        # The program prints out the header of the sequence being analyzed
        print NewOT "$header[$x]\n";
        
        # The length of the current sequence is calculated and printed
        $lengthSequence = length $sequence[$x];
        print NewOT "Length:$lengthSequence\n";
        
        # The current sequence is saved in a variable so it can be used in the following for loop for calculations
        $currentSequence = $sequence[$x];
        
        my @cpg = ($currentSequence =~ /CG/g);    # Based on the current sequence the number of CpGs are calculated and proportion is calculated
        my $d = @cpg;	    					  # saved as a scalar variable to print out later
        my $percent = @cpg / $lengthSequence;     # proportion of CpG calculated        

        
        # This for loop iterations through every base in the sequence counting how many of each base exist and adding them up in their individual variables
        for (my $b = 0; $b < $lengthSequence; $b++){
            
            $sub = substr($currentSequence, $b, 1);     # variable holds the value of the individual sequence being analyzed
            
            if ($sub eq 'G'){
                $gc += 1;
            }
            
            elsif ($sub eq 'C'){
                $cc += 1;
            }
            
            elsif ($sub eq 'A'){
                $ac += 1;
            }		
            
            elsif ($sub eq 'T'){
                $tc += 1;
            }
            
            # The proportion of each base is calculated and saved for G / C / A / T
            
            $gp = $gc / $lengthSequence;
            $cp = $cc / $lengthSequence;
            $ap = $ac / $lengthSequence;
            $tp = $tc / $lengthSequence;
        }
        
        #Using the correct formatting printf command the number of each base and the proportion are printed sequentially
        printf NewOT "A:$ac %0.2f\n", $ap;
        printf NewOT "C:$cc %0.2f\n", $cp;
        printf NewOT "G:$gc %0.2f\n", $gp;
        printf NewOT "T:$tc %0.2f\n", $tp;
        printf NewOT "CpG:$d %0.2f\n", $percent;

        #saving the counts to the arrays so they can be returned to the user
        @$acounts[$x] = $ac;
        @$tcounts[$x] = $tc;
        @$gcounts[$x] = $gc;
        @$ccounts[$x] = $cc;
        @$cgcounts[$x] = $d;

        #saving the proportions to arrays so they can returned to user
        @$aprops[$x] = $ap;
        @$tprops[$x] = $tp;
        @$gprops[$x] = $gp;
        @$cprops[$x] = $cp;
        @$cgprops[$x] = $percent;

    }
    close(NewOT);

    # dereferencing the arrays
    my @acounts = @$acounts;
    my @tcounts = @$tcounts;
    my @gcounts = @$gcounts;
    my @ccounts = @$ccounts;
    my @cgcounts = @$cgcounts;
    my @aprops = @$aprops;
    my @tprops = @$tprops;
    my @gprops = @$gprops;
    my @cprops = @$cprops;
    my @cgprops = @$cgprops;
}

sub permute_fasta {

    for (my $x = 0; $x < $count; $x++){

        # returns the length of the current sequence to use for generating the random sequence
        $lengthSequence = length $sequence[$x];

        # The current sequence is saved in a variable so it can be used in the following for loop for calculations
        my $currentSequence = $sequence[$x];
        my @currentSequence = split(//, $currentSequence);
        
        # saves the new shuffled numerical sequence to an array
        my @randomArray = (shuffle 0 .. $lengthSequence-1)[0..$lengthSequence-1];
        my $lenrandarray = scalar(@randomArray);

        my @randomSequence = (); # will hold the new random sequence in this array
        my $count = 0; # initialize count variable for loop

        # loop that will parse through the current sequence and save the new random one to @randomSequence array
        for(@randomArray){

            my $num = $_;

            @randomSequence[$count] = $currentSequence[$num];

            $count += 1;

        }

        # join the array together in to a scalar variable so it can be added to the array
        my $scal = join("", @randomSequence);       
        
        @newRandomSequences[$x] = $scal;

    }

    @sequence = @newRandomSequences;
}

sub write_fasta {
    
    my(@sequnce) = @_;
    
    #creating array so I can extract the filename without the suffix to create .ot file
    my @newfilename = split('\.', $filename, -1);

     # creating new output file name with .ot suffix
    $permuteFile = $newfilename[0] . "_permute.fsa";   
    
    
    # write the following to the new filename .ot
    open(PermuteFSA, '>', $permuteFile) or die $!;

    #counter for for loop to correctly use the right header for each sequence
    my $counter = 0;

    for(@sequence){
        print PermuteFSA "$header[$counter] \n";
        print PermuteFSA "$_ \n";
        $counter += 1;
    }

    close(PermuteFSA);
    
}

read_fasta($filename);
stat_fasta(\@gprops, \@cprops, \@aprops, \@tprops, \@cgprops, \@acounts, \@tcounts, \@gcounts, \@ccounts, \@cgcounts, $filename);
permute_fasta(@sequence);
write_fasta(\@sequence, $filename);
stat_fasta(\@gprops, \@cprops, \@aprops, \@tprops, \@cgprops, \@acounts, \@tcounts, \@gcounts, \@ccounts, \@cgcounts, $permuteFile);


exit;