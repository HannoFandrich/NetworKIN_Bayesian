#!/usr/bin/perl -w

use strict;

my $rootdir = "..";


open CODE, "> nn_code.h";
open DATA, "> nn_data.h";
if (-e "$rootdir/evaluation_run.tab") {
  open EVAL, "< $rootdir/evaluation_run.tab";
  while (<EVAL>) {
    s/\r?\n//;
    next if /^#/;
    my (undef, $name, $source, $org, $type, $res, undef, $a0, $a1, $a2, $a3, undef, undef, $prior, undef) = split /\t/;
    next unless $source eq "nn";
    if ($res eq "STY") {
      print CODE "if (c == 'S' || c == 'T' || c == 'Y') {\n";
    } elsif ($res eq "ST") {
      print CODE "if (c == 'S' || c == 'T') {\n";
    } elsif ($res eq "SY") {
      print CODE "if (c == 'S' || c == 'Y') {\n";
    } elsif ($res eq "TY") {
      print CODE "if (c == 'T' || c == 'Y') {\n";
    } elsif ($res eq "S") {
      print CODE "if (c == 'S') {\n";
    } elsif ($res eq "T") {
      print CODE "if (c == 'T') {\n";
    } elsif ($res eq "Y") {
      print CODE "if (c == 'Y') {\n";
    } else {
      next;
    }
    my $enum = 0;
    print CODE "  o = 0;\n";
    open LIST, "< ".$rootdir."/nn_".$org."_".$type."_".$name."/syn.lst";
    while (<LIST>) {
      s/\r?\n//;
      next unless /([0-9]+)_[0-9.]+_([0-9]+)\.syn$/;
      $enum++;
      print CODE "  o += feed_forward(s, nn_", $org, "_", $type, "_", $name, "_", $enum, ", ", $1, ", ", $2, ");\n";
      open SYN, "< $_";
      for (my $i = 0; $i < 9; $i++) {
        $_ = <SYN>;
      }
      print DATA "float const nn_", $org, "_", $type, "_", $name, "_", $enum, "[] = {\n";
      while (<SYN>) {
        s/\r?\n//;
        s/^     /  /;
        s/     /, /g;
        print DATA $_, ",\n";
      }
      print DATA "};\n\n";
      close SYN;
    }
    close LIST;
    print CODE "  o /= ", $enum, ";\n";
    print CODE "  o = $a2+($a3-$a2)/(1+exp($a1*($a0-o)));\n";
    #print CODE "  if (o > ", $prior, ") {\n";
    print CODE "    printf(\"%s\\t%d\\t%c\\t\", name, *n, c);\n";
    print CODE "    print_peptide(s);\n";
    print CODE "    printf(\"\\tnn\\t", $org, "\\t", $type, "\\t", $name,"\\t%.6f\\t%.6f\\n\", o, ", $prior, ");\n";
    #print CODE "  }\n";
    print CODE "}\n\n";
  }
  close EVAL;
}
close CODE;
close DATA;
