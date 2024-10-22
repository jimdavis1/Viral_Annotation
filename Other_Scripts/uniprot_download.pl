#!/usr/bin/env perl
use strict;
use warnings;
use LWP::UserAgent;
use Getopt::Long;

my $usage = 'uniprot_download.pl -q [query] >Output.fasta
        -h   help
        -q   uniprot query eg "Arenaviridae"
        -r   "reviewed entries only (this downloads only SwissProt)
        -j   return json (default = fasta)
';

my ($help, $in_query, $reviewed, $json);
my $format = "fasta"; 
GetOptions(
    'h'   => \$help,
    'q=s' => \$in_query,
    'r'   => \$reviewed,
    'j'   => \$json,
) or die "$usage\n";

if ($help) { die "$usage\n"; }
if (!$in_query) { die "Must declare a query\n$usage\n"; }
if ($json){$format = "json";}



my $ua = LWP::UserAgent->new;

my $query = "($in_query)";
if ($reviewed) { $query = "($in_query) AND (reviewed:true)"; }

my $base_url = "https://rest.uniprot.org/uniprotkb/stream";

#curl "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28Arenaviridae%29" >Lassa.all.prots

my $url = "$base_url?compressed=false&format=$format&query=$query";
    
my $response = $ua->get($url);
my $data;

if ($response->is_success) 
{
	$data = $response->decoded_content;
}
else 
{
	die "Error fetching data: ", $response->status_line;
}

print $data;
