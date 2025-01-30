#!/usr/bin/perl -w

use strict;

my $rootdir = "..";




#
# Create C-code.
#

open CODE, "> insr_code.h";

print CODE "if (c == 'Y') {\n";
print CODE "  o = 0;\n";
print CODE "  o += feed_forward(s, netphosk_InsR_1, 17, 8);\n";
print CODE "  o += feed_forward(s, netphosk_InsR_2, 9, 8);\n";
print CODE "  o += feed_forward(s, netphosk_InsR_3, 7, 4);\n";
print CODE "  o /= 3;\n";
if (-e "$rootdir/evaluation_run.tab") {
  open EVALUATION, "< $rootdir/evaluation_run.tab";
  while (<EVALUATION>) {
    s/\r?\n//;
    next if /^#/;
    my ($name1, $name2, $source, $type, $res, undef, $a0, $a1, $a2, $a3, undef, undef, $prior, undef) = split /\t/;
    next unless $source eq "netphosk" and $name1 eq "InsR_group";
    print CODE "  o = $a2+($a3-$a2)/(1+exp($a1*($a0-o)));\n";
    #print CODE "  if (o > ", $prior, ") {\n";
    print CODE "    printf(\"%s\\t%d\\t%c\\t\", name, *n, c);\n";
    print CODE "    print_peptide(s);\n";
    print CODE "    printf(\"\\tnetphosk\\tKIN\\tInsR_group\\t%.6f\\t%.6f\\n\", o, ", $prior, ");\n";
  }
  close EVALUATION;
} else {
  #print CODE "  if (o > 0) {\n";
  print CODE "    printf(\"%s\\t%d\\t%c\\t\", name, *n, c);\n";
  print CODE "    print_peptide(s);\n";
  print CODE "    printf(\"\\tnetphosk\\tKIN\\tInsR_group\\t%.6f\\t\\n\", o);\n";
}
#print CODE "  }\n";
print CODE "}\n";
close CODE;
