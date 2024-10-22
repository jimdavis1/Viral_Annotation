#! /usr/bin/env perl
use strict;
use Getopt::Long;

my $usage = 'list_pssms.pl -b base_directory
		
        This program reads the directory of pssms and provides the complete
        list of pssms for every taxonomic unit.  Note that the annotation field  
        that is displayed from the pssm is deprecated (but still generally informative). 
        
        In earlier versions of the annotation script, the annotation was read directly
        from the pssm. I have changed this to be read directly from the json file to 
        simplify the  management of the annotations. 

        Usage:
		-b base directory (path preceeding /Viral_Annotation)
		   in my case its -b "/home/jjdavis/bin"
		-h help';

my ($help, $base); 		
my $opts = GetOptions( 'h'         => \$help,
                       'b=s'       => \$base);

unless ($base){die "must declare -b base directory preceeding Viral_Annotaton\n$usage\n";}
if ($help){die "$usage\n"; }
		
opendir (DIR, "$base/Viral_Annotation/Viral-PSSMs/"); 
my @fams = grep{$_ !~ /^\./}readdir(DIR);  
close DIR;		

foreach (@fams)
{
	my $fam = $_;  
	opendir (DIR, "$base/Viral_Annotation/Viral-PSSMs/$fam/"); 
	my @prots = grep{$_ !~ /^\./}readdir(DIR);  
	close DIR; 
	
	foreach (@prots)
	{
		my $prot = $_; 
		opendir (DIR, "/home/jjdavis/bin/Viral_Annotation/Viral-PSSMs/$fam/$prot/"); 
		my @pssms = grep{$_ !~ /^\./}readdir(DIR); 
		close DIR; 
		
		foreach (@pssms)
		{ 
			my $pssm = $_; 
			open (IN, "grep title /home/jjdavis/bin/Viral_Annotation/Viral-PSSMs/$fam/$prot/$pssm |"); 
			while (<IN>)
			{
				chomp; 
				s/.+title \"//g; s/\"$//g; 
				print "$fam\t$prot\t$pssm\t$_\n";
			}
		}
	}
}
