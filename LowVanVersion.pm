package LowVanVersion;
use strict;
use warnings;

# Try to load the generated version file; fall back to "dev" if not present
our $VERSION = "dev";

eval {
    require LowVanVersionData;
    $VERSION = $LowVanVersionData::VERSION if defined $LowVanVersionData::VERSION;
};

sub get_version {
    return $VERSION;
}

1;
