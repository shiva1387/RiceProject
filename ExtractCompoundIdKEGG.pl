#!/usr/bin/env perl

use strict;
die "USAGE: $0 KEGGftp/ligand/compound/compound" unless $#ARGV == 0;

$/ = '///';

print "Extracting compound ID, name, exact MASS, molecular weight\n";
open COMPOUND, $ARGV[0] || die $!;
while(<COMPOUND>) {
    /^ENTRY\s+      (?<CPD>C\d{5}).*
    ^NAME\s+        (?<NAME> .*?)\n.*
    ^EXACT_MASS\s+  (?<EXACT>\S+).*
    ^MOL_WEIGHT\s+  (?<MOL>\S+)\n
    /xsm;
    print join("\t", "cpd:$+{CPD}", $+{NAME}, $+{EXACT}, $+{MOL}),"\n" unless (length($+{CPD}) == 0);
}
close COMPOUND;

#open GLYCAN, $ARGV[1] || die $!;
#while(<GLYCAN>) {
    #my ($cpd) = $_ =~ /ENTRY\s+(G\d{5})/xsm;
    #/NAME\s+(?<NAME> .*?)\n  	#I'm just taking the first symbol name associated with this
    #/xsm;
    #print "gl:$cpd\t$+{NAME}\n" unless (length($cpd) == 0)
#}
#close GLYCAN;


__DATA__

ENTRY       C00001                      Compound
NAME        H2O;
            Water
FORMULA     H2O
EXACT_MASS  18.0106
MOL_WEIGHT  18.0153
REMARK      Same as: D00001
REACTION    R00001 R00002 R00004 R00005 R00009 R00010 R00011 R00017
            R00022 R00024 R00025 R00026 R00028 R00036 R00041 R00044
            R00045 R00047 R00048 R00052 R00053 R00054 R00055 R00056
