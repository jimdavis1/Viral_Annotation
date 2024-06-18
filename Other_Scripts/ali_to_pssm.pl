#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;
use gjoseqlib;
use BlastInterface;

my $usage = 'ali_to_pssm.pl [opts] < alignment > pssm
              -h help
              -t title for pssm (optional)
';


my ($help, $title);
my $opts = GetOptions( 'h'         => \$help,
                       't=s'       => \$title);

if ($help){die "$usage\n";}

my @ali = &gjoseqlib::read_fasta();

if ($title)
{
	my %opts = (pssm_title => $title); 
	my $pssm = BlastInterface::alignment_to_pssm(\@ali, \%opts);
}

else
{
	my $pssm = BlastInterface::alignment_to_pssm(\@ali);
}

