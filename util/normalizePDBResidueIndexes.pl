#!/usr/bin/perl -w
use strict;
my ($utilDirectory) = $0 =~ m/(.*\/)/g;
require $utilDirectory."PDBNormalizer.pm";
my @pdb = readPdb($ARGV[0]);
writeNormPDB(\@pdb);
print "END\n";


