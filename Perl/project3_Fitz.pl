use strict;
use warnings;

# print matrix and identify matrix subroutines and code were derived from Dr. Solka's lecture

my @matrix = (); # for holding the matrix
my @matrixB = (); # for holding the matrix
my @current = (); # to convert the string of numbers read in from the file to an array to loop through
my $col = 0; # initializing column variable for loop
my $row = 0;			    # index of current row
my @A = ();
my @B = ();
my $matlen = 0;
my @multiplied = (); # multiplied matrix

print("Please enter the filename for Matrix A: ");
my $filename = <STDIN>;

read_matrix($filename);           # will read in the matrix with the user entered filename A
mult_matrix(\@matrix, \@matrix); # will multiply both matrix's created together
print("\n Prodcut of Matrix A and Matrix B that were entered by user. \n");
print_2Darray(@multiplied);       # prints out the multiplied matrix
print("Please enter n (the power to multiply matrix A to): ");
my $n = <STDIN>;
power_matrix(\@matrix, $n);


sub read_matrix{

    # Open the file given as the first argument on the command line
    open(INFILE, $filename) or die "Can't open $filename\n";
    
    # read text file that contains matrix
    while (my $line = <INFILE>) {
        chomp $line;		    # remove training \n from line

        @current = split / /, $line; # splits the current line from the matrix file by the spaces
        my $length = scalar(@current); #returns the length of the row to use in the for loop

        # loops through all of the columns in the current row being analyzed and adds those values
        # to the @matrix() with the while loop keeping track of what row the program is on and the 
        # for loop advancing the columns
        for($col = 0; $col < $length; $col++){
            $matrix[$row][$col] = $current[$col];
        }

        $row++; # moves to the next row
    }
    close INFILE;
    return @matrix, $row, $col;

}


sub mult_matrix{

    my $r = 0; # variable for current row
    my $c = 0; # variable for current loop
    my $total = 0; # running total for adding multiple multiplications
    my $mult = 0; # variable to use inside matrix multiplication
    my ($A, $B) = @_;
    my @D = @{$A}; # reference by pointer
    my @F = @{$B}; # reference by pointer

    my $var = 5;
    my $t = 0;
    my $lenghtest = 0;
    my $col = 0; # initializing column variable for loop
    my $row = 0;			    # index of current row
    my $rowconvert = 0; #use this to count columns in @B to calculate the number of rows
    my @T = ();
    my @bcollength = ();
    my $total_count = 0;
    my $bcolLength = 0;

    # the following two for loops are used to compare the rows and columns from both input matrix's are compatible for multiplication
    for(@D){
        $rowconvert++;
        @T = split /,/, @$_;
        $total_count = $T[0] + $total_count;
    }

    for(@F){
        $row++;
        @bcollength = split /,/, @$_;
        $bcolLength = $bcollength[0] + $bcolLength;

    }

    $col = $total_count / $rowconvert;
    my $colB = $bcolLength / $row;

    if($col == $row){

        # loop will iterate through all of the rows in the matrix
        for($r = 0; $r < $rowconvert; $r++){
            
            # loops will iterate through all of the columns in the matrix
            for($c = 0; $c < $colB; $c++){

                # for multiplication it will make sure to include all rows
                while($mult < $row){

                    # holds the running total
                    $total = $total + $D[$r][$mult] * $F[$mult][$c];

                    # counter adds one each time it iterates through
                    $mult++;

                }
                
                # new value is appended to the multiplied matrix
                $multiplied[$r][$c] = $total;
                
                $total = 0; # total is reset to 0 for next loop
                $mult = 0; # mult is reset to 0 for next loop

            }

        }

    }
    else{
        print("The number of columns in the first matrix should be equal to the number of rows in the second.");
    }

    return @multiplied;

}


sub power_matrix {

    my ($input, $n) = @_; # will accept inputs from the user entered power and the matrix we are finding the power of (@A) by default
    my @input = @{$input}; # reference by pointer
    my @C = ();
    print("Matrix to the power of: ", $n, "\n");

    if($n == 0){

        # below code was created by Dr. Solka to create the inverse matrix
        my @inv = ();
		my $inrows = scalar(@input);
		my $incols = scalar(@input);
		#create a matrix with 1's on diagonal
		for (my $w = 0; $w < $inrows; $w++) {
			for (my $j = 0; $j < $incols; $j++) {
				$inv[$w][$j] = 0;
			}
				$inv[$w][$w] = 1;
		}

    print_2Darray(@inv);
    }

    # if n = 1 it will return the same matrix
    elsif($n == 1){
        print_2Darray(@input);
    }

    else{
        my $count = 2;
        @C = mult_matrix(\@input, \@input); # at a minimum if it makes it to this part of the loop it will be squaring the matrix's and if it is higher than squared it will
                                    # go in to the while loop
        
        # continues to iterate through until the power equals the user entered power
        while($count < $n){

            @C = mult_matrix(\@C, \@input);
            $count++;
        }
    print_2Darray(@C);
    }
}

# taken from Dr. Solka's lecture
# I added **** for formatting around the matrix's for better viewing in console
sub print_2Darray {
    my (@a) = @_;
    my $rows = scalar @a;
    my $cols = scalar @{$a[0]};
    my $test = ("**********");
    my $g = 0;
    while($g < $cols){
        my $t = ("********");
        $test = $test . $t;
        $g++;
    }
    print($test, "\n");
    for (my $i=0; $i < $rows; $i++) {
        print("* ");
	   for (my $j=0; $j < $cols; $j++) {
	     printf "%8.2f ", $a[$i][$j];
	   }
	   print "* \n"; # newline after each row
    }
    print($test, "\n");
}

