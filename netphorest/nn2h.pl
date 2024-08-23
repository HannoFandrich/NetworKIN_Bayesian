#!/usr/bin/perl -w

use strict;

my $rootdir = "/mnt/data/NetPhorest/v2.1";


open CODE, "> nn_code.h";
open DATA, "> nn_data.h";
if (-e "$rootdir/evaluation_run.tab") {
  open EVAL, "< $rootdir/evaluation_run.tab" or die "Cannot open $rootdir/evaluation_run.tab: $!"; 
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
    print CODE "  gbHas_unsupported_amino_acid = 0;\n";

    my $synlst = $rootdir."/nn_".$org."_".$type."_".$name."/syn.lst";
    open LIST, "<", $synlst or die "Cannot open syn.lst file $synlst: $!";
    while (<LIST>) {
      s/\r?\n//;
      next unless /([0-9]+)_[0-9.]+_([0-9]+)\.syn$/;
      my $synapse_file = $_;
      $enum++;
      my $array_var_name = "nn_" . $org . "_" . $type . "_" . $name . "_" . $enum;
      my ($nw, $nh) = ($1, $2);
      print CODE "  o += feed_forward(s, $array_var_name", ", ", $nw, ", ", $nh, ");\n";
      open SYN, "< $_" or die "Cannot open synapse file $_: $!";
      print DATA "float const ${array_var_name}[] = {\n";
      my $nskip = 0;
      while (<SYN>) {
        $nskip++, next if /[A-Z]/;
        s/\r?\n//;
        s/^     /  /;
        s/     /, /g;
        print DATA $_, ",\n";
      }
      print DATA "};\n";
      if ($nh == 0) {
	      print DATA "STATIC_ASSERT(sizeof($array_var_name)/sizeof(${array_var_name}[0]) == 21 * $nw * 2 + 2, \"$array_var_name has incorrect size\");";
      }
      else {
	      print DATA "STATIC_ASSERT(sizeof($array_var_name)/sizeof(${array_var_name}[0]) == 21 * $nw * $nh + $nh * 2 + $nh + 2, \"$array_var_name has incorrect size\");";
      }
      print DATA "\n\n";
      close SYN;
      unless (($2 == 0 && $nskip == 8) || ($2 > 0 && $nskip == 9)) {
	      die sprintf("Skipped an unexpected number of lines (%d) in %s!", $nskip, $synapse_file);
      }
    }
    close LIST;
    print CODE "  o /= ", $enum, ";\n";
    print CODE "  o = $a2+($a3-$a2)/(1+exp($a1*($a0-o)));\n";
    #print CODE "  if (o > ", $prior, ") {\n";
    print CODE "  if (gbHas_unsupported_amino_acid == 0){\n";
    print CODE "    printf(\"%s\\t%d\\t%c\\t\", name, *n, c);\n";
    print CODE "    print_peptide(s);\n";
    print CODE "    printf(\"\\tnn\\t", $org, "\\t", $type, "\\t", $name,"\\t%.6f\\t%.6f\\n\", o, ", $prior, ");\n";
    print CODE "  }\n";
    #print CODE "  }\n";
    print CODE "}\n\n";
  }
  close EVAL;
}
close CODE;
close DATA;
