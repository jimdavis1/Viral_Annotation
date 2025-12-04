#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;

# Parse command line options
my %opts;
getopts('hb', \%opts);

# Show help if -h flag is used or no arguments
if ($opts{h} || @ARGV != 1) {
    print_usage();
    exit(0);
}

my $filename = $ARGV[0];

# Check if file exists
unless (-f $filename) {
    die "Error: File '$filename' does not exist.\n";
}

# Create backup if -b flag is used
if ($opts{b}) 
{
    my $backup = $filename . ".bak";
    rename($filename, $backup) or die "Cannot create backup '$backup': $!\n";
    
    # Read the backup file
    open(my $fh, '<', $backup) or die "Cannot open file '$backup': $!\n";
    my $content = do { local $/; <$fh> };
    close($fh);
    
    # Process content
    $content = fix_pssm_content($content);
    
    # Write to original filename
    open(my $out_fh, '>', $filename) or die "Cannot write to '$filename': $!\n";
    print $out_fh $content;
    close($out_fh);
    
} 
else 
{
    # Read the file
    open(my $fh, '<', $filename) or die "Cannot open file '$filename': $!\n";
    my $content = do { local $/; <$fh> };
    close($fh);
    
    # Process content
    $content = fix_pssm_content($content);
    
    # Write back to the same file
    open(my $out_fh, '>', $filename) or die "Cannot write to '$filename': $!\n";
    print $out_fh $content;
    close($out_fh);
}

# Subroutine to fix PSSM content
sub fix_pssm_content {
    my ($content) = @_;
    
    # Pattern to match the two consecutive descr blocks
    my $pattern = qr/
        (descr\s*\{\s*title\s*")      # Start of first descr (capture group 1)
        ([^"]*")                       # Content of first title (capture group 2)
        (\s*\}\s*,\s*)                # End of first descr (capture group 3)
        (descr\s*\{\s*title\s*")      # Start of second descr (capture group 4)
        (1\s+[^"]*")                   # Content of second title starting with "1 " (capture group 5)
        (\s*\}\s*,?)                   # End of second descr (capture group 6)
    /sx;
    
    # Replace the matched pattern
    $content =~ s/$pattern/$4$5$6/g;
    
    # Now remove the "1 " prefix from the remaining title
    $content =~ s/(descr\s*\{\s*title\s*")1\s+/$1/g;
    
    return $content;
}



# Subroutine to print usage information
sub print_usage {
    print <<'USAGE';
Usage: clean_pssm.pl [OPTIONS] <pssm_file>

Description:
    Fixes PSSM files by removing duplicate/incorrect description (descr) fields.
    The script removes the first descr block (typically not starting with "1")
    and keeps the second one (typically starting with "1"), then removes the
    "1 " prefix from the remaining description.

Options:
    -h    Show this help message
    -b    Create a backup file (.bak) before modifying the original

Examples:
    clean_pssm.pl myfile.pssm              # Fix file in place (no backup)
    clean_pssm.pl -b myfile.pssm           # Fix file and create backup
    clean_pssm.pl -h                       # Show this help

USAGE
}