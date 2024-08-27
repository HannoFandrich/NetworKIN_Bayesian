#!/usr/bin/perl -w

use strict;

my $rootdir = "/mnt/data/NetPhorest/v2.1";


my $alphabet = "FIVWMLCHYAGNRTPDEQSK";
my @alphabet = split //, $alphabet;
my %alphabet = (
  "PHE" => 0,
  "ILE" => 1,
  "VAL" => 2,
  "TRP" => 3,
  "MET" => 4,
  "LEU" => 5,
  "CYS" => 6,
  "HIS" => 7,
  "TYR" => 8,
  "PTR" => 8,
  "ALA" => 9,
  "GLY" => 10,
  "ASN" => 11,
  "ARG" => 12,
  "THR" => 13,
  "TPO" => 13,
  "PRO" => 14,
  "ASP" => 15,
  "GLU" => 16,
  "GLN" => 17,
  "SER" => 18,
  "SEP" => 18,
  "LYS" => 19
);




#
# Read command line options.
#

my $private;
if (defined $ARGV[0] and $ARGV[0] eq "private" and defined $ARGV[1]) {
  $private = 1;
} elsif (defined $ARGV[0] and $ARGV[0] eq "public" and defined $ARGV[1]) {
  $private = 0;  
}
else {
  die "Syntax: pssm2h.pl [private|public] [datadir]\n";
}

my $datadir = $ARGV[1];


#
# Read list of private matrices.
#

my %source_org_type_name_private = ();
open PRIVATE, "< $datadir/PRIVATE.tsv" or die "Cannot open $datadir/PRIVATE.tsv: $!";
while (<PRIVATE>) {
  s/\r?\n//;
  my ($source, $org, $type, $name) = split /\t/;
  $source_org_type_name_private{$source}{$org}{$type}{$name} = 1;
}
close PRIVATE;




#
# Read sigmoids.
#

my %source_org_type_name_name = ();
my %source_org_type_name_prior = ();
my %source_org_type_name_sigmoid = ();
if (-e "$rootdir/evaluation_run.tab") {
  open EVALUATION, "< $rootdir/evaluation_run.tab" or die "Cannot open $datadir/evaluation_run.tab: $!";
  while (<EVALUATION>) {
    s/\r?\n//;
    next if /^#/;
    my ($name1, $name2, $source, $org, $type, $res, undef, $a0, $a1, $a2, $a3, undef, undef, $prior, undef) = split /\t/;
    $source_org_type_name_name{$source}{$org}{$type}{$name2} = $name1;
    $source_org_type_name_prior{$source}{$org}{$type}{$name2} = $prior;
    $source_org_type_name_sigmoid{$source}{$org}{$type}{$name2} = "o = $a2+($a3-$a2)/(1+exp($a1*($a0-o)));\n";
  }
  close EVALUATION;
}




#
# Read pELM name mappings.
#

my %hugo_pelm = ();
open PELM, "< hugo_pelm.tsv" or die "Cannot open hugo_pelm.tsv: $!";
while (<PELM>) {
  s/\r?\n//;
  my ($hugo, $pelm) = split /\t/;
  $hugo_pelm{$hugo} = $pelm;
}
close PELM;




#
# Read matrices from the Scansite database.
#

my %source_org_type_name_pos_res_score = ();
{
  open SCANSITE, "< $datadir/scansite_matrix.tsv" or die "Cannot open $datadir/scansite_matrix.tsv: $!";
  $_ = <SCANSITE>;
  while (<SCANSITE>) {
    s/\r?\n//;
    my ($type, $name, $pos, @scores) = split /\t/;
    $name =~ s/[ -]//g;
    foreach my $res (sort @alphabet) {
      $source_org_type_name_pos_res_score{"scansite"}{"human"}{$type}{$name}{$pos}{$res} = shift(@scores);
    }
  }
  close SCANSITE;
}




#
# Read additional matrices in Scansite format.
#

foreach my $file (split /\n/, `ls -1 $datadir`) {

  next unless $file =~ /^([^.]+)\.([^.]+)\.([^.]+)\.(scansite(_[0-9A-Za-z]+)?)$/;
  my $org = $1;
  my $type = $2;
  my $name = $3;
  my $source = $4;
  $name =~ s/[ -]//g;

  open SCANSITE, "< $datadir/$file" or die "Cannot open scansite matrix file $datadir/$file: $!";

  # Read header line.
  $_ = <SCANSITE>;
  s/\r?\n//;
  my @head = split /[\t ]/;

  # Read data lines.
  {
    my $warn = 0;
    my $i = 0;
    while (<SCANSITE>) {
      s/\r?\n//;
      my $count = 0;
      my @data = split /[\t ]/;
      for (my $j = 0; $j <= $#data; $j++) {
        next unless exists $head[$j] and $head[$j] =~ /^ *([$alphabet]) *$/;
        $count++;
        $source_org_type_name_pos_res_score{$source}{$org}{$type}{$name}{$i}{$1} = $data[$j];
      }
      $warn = 1 if $count < 20;
      $i++;
    }
    warn "$datadir/$file does not appear to be in Scansite format!\n" if $warn;
  }

  close SCANSITE;

}




#
# Create C-code.
#

open CODE, "> pssm_code.h";
open DATA, "> pssm_data.h";

srand(10001);
foreach my $source (sort keys %source_org_type_name_pos_res_score) {
  foreach my $org (sort keys %{$source_org_type_name_pos_res_score{$source}}) {
    foreach my $type (sort keys %{$source_org_type_name_pos_res_score{$source}{$org}}) {
      foreach my $name (sort keys %{$source_org_type_name_pos_res_score{$source}{$org}{$type}}) {
        next if not $private and exists $source_org_type_name_private{$source} and exists $source_org_type_name_private{$source}{$org} and exists $source_org_type_name_private{$source}{$org}{$type} and exists $source_org_type_name_private{$source}{$org}{$type}{$name};
        next if -e "$rootdir/evaluation_run.tab" and not (exists $source_org_type_name_sigmoid{$source} and exists $source_org_type_name_sigmoid{$source}{$org} and exists $source_org_type_name_sigmoid{$source}{$org}{$type} and exists $source_org_type_name_sigmoid{$source}{$org}{$type}{$name});

        # Generate random score distribution.
        my @scores = ();
        for (my $i = 0; $i < 100000; $i++) {
          my $score = 1;
          foreach my $pos (keys %{$source_org_type_name_pos_res_score{$source}{$org}{$type}{$name}}) {
            my @res = keys %{$source_org_type_name_pos_res_score{$source}{$org}{$type}{$name}{$pos}};
            $score *= $source_org_type_name_pos_res_score{$source}{$org}{$type}{$name}{$pos}{$res[int rand scalar @res]};
          }
          push @scores, $score;
        }
        @scores = sort {$b <=> $a} @scores;

        # Determine normalization factor.
        my $norm = 0;
        for (my $i = 0; $i < 10 and $scores[$i] > 0; $i++) {
          $norm = $scores[$i];
        }
        $norm = sprintf "%.3e", $norm;

        # Create C-code.
        print DATA "float const ", $source, "_", $org, "_", $type, "_", $name, "_pssm[] = {\n";
        my @pos = sort {$a <=> $b} keys %{$source_org_type_name_pos_res_score{$source}{$org}{$type}{$name}};
        my $size = 0;
        foreach my $pos (@pos) {
          $size++;
          my $min = 1;
          foreach my $res (@alphabet) {
            my $score = $source_org_type_name_pos_res_score{$source}{$org}{$type}{$name}{$pos}{$res};
            $min = $score if $score < $min;
            printf DATA "%.3f, ", $score;
          }
          printf DATA "%.3f", $min;
          print DATA "," unless $pos == scalar(@pos)-1;
          print DATA "\n";
        }
        print DATA "};\n\n";
        print CODE "o = pssm(s, ", $source, "_", $org, "_", $type, "_", $name, "_pssm, ", $size, ")/", $norm, ";\n";
        print CODE "if (o > 0) {\n";
        print CODE "  o = log(o);\n";
        if (-e "$rootdir/evaluation_run.tab") {
          print CODE "  ", $source_org_type_name_sigmoid{$source}{$org}{$type}{$name};
          #print CODE "  if (o > ", $source_org_type_name_prior{$source}{$org}{$type}{$name}, ") {\n";
          print CODE "    printf(\"%s\\t%d\\t%c\\t\", name, *n, c);\n";
          print CODE "    print_peptide(s);\n";
          print CODE "    printf(\"\\t",$source, "\\t", $org, "\\t", $type, "\\t", $source_org_type_name_name{$source}{$org}{$type}{$name}, "\\t%.6f\\t%.6f\\n\", o, ", $source_org_type_name_prior{$source}{$org}{$type}{$name}, ");\n";
          #print CODE "  }\n";
        } else {
          print CODE "  printf(\"%s\\t%d\\t%c\\t\", name, *n, c);\n";
          print CODE "  print_peptide(s);\n";
          print CODE "  printf(\"\\t", $source, "\\t", $org, "\\t", $type, "\\t", $name, "\\t%.6f\\t\\n\", o);\n";
        }
        print CODE "}\n\n";

      }
    }
  }
}

close CODE;
close DATA;
