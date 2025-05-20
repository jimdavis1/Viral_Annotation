package FCP_PSSM_Utils;

use strict;
use warnings;
use Exporter qw(import);
use gjoseqlib;
use BlastInterface;

our @EXPORT_OK = qw(create_pssm process_pssm_file);

# Create a PSSM from an alignment
sub create_pssm 
{
    my ($alignment, $options) = @_;
    
    my $out_file = $options->{outfile} || "pssm.out";
    my $pssm_prefix = $options->{pssm_prefix} || "";
    my $tmp = $options->{tmp} || "tmp";
    
    # Call BlastInterface to create the PSSM
    my %pssm_opts = (outPSSM => $out_file);
    my $pssm = BlastInterface::alignment_to_pssm($alignment, \%pssm_opts);
    
    # Process the PSSM file to fix titles if necessary
    if ($options->{process_titles}) 
    {
        my $final_file = $options->{final_file} || "final_pssm.out";
        process_pssm_file($out_file, $final_file);
        
        # Clean up temporary file if requested
        if ($options->{cleanup}) 
        {
            unlink($out_file);
        }
    }
    
    return $pssm;
}

# Process a PSSM file to fix titles
sub process_pssm_file 
{
    my ($input_file, $output_file) = @_;
    
    open(my $in_fh, "<$input_file") or die "Cannot open PSSM file $input_file: $!\n";
    open(my $out_fh, ">$output_file") or die "Cannot open output PSSM file $output_file: $!\n";
    
    my $des = 0;
    
    # BlastInterface is useful for making the pssm, but it doesn't format the title
    # nicely, so below we loop over that, and fix it. Basically, it enumerates 
    # each title, which is ugly, and if you don't give it a title, it will give you an 
    # undefined title plus your annotation string, so we just pick it up from the 
    # annotation string, and delete the undefined one that is printed compulsively.
    while (my $line = <$in_fh>) 
    {
        chomp $line;
        if (($des > 1) && ($line =~ /title/)) 
        {
            next;
        } 
        elsif (($des > 1) && ($line =~ /\}/)) 
        {
            $des = 0;
            next;
        } 
        elsif ($line =~ /descr \{/) 
        {
            $des++;
            if ($des <= 1) 
            {
                print $out_fh "$line\n";
            }
        } 
        elsif ($des <= 1) 
        {
            print $out_fh "$line\n";
        }
    }
    
    close($in_fh);
    close($out_fh);
    return 1;
}

1; # End of module