# Various perl routines
# Rok Zitko, zitko@theorie.physik.uni-goettingen.de, Oct 2008

use strict;
use warnings;

# Read a column from a file
sub readcol {
    my ($filename, $column) = @_;
    
    my @l;
    open(F, "<$filename") or die "Can't open $filename for reading, stopped";
    while (<F>) {
	chomp;
	my $x = (split)[$column-1];
	push(@l, $x);
    }
    close(F);
    
    @l
};

1;
