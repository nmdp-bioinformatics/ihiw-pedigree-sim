#!/usr/bin/env perl
############################################################################
# SCRIPT NAME:  pedigree_sim_children.pl
# DESCRIPTION:  simulate children genotypes from parental genotypes
#
# DATE WRITTEN: July 21, 2019
# WRITTEN BY:   Loren Gragert
#
# REVISION HISTORY:
# REVISION DATE         REVISED BY      DESCRIPTION
# ------- ----------    --------------  -------------------------------------
#
##############################################################################
use strict; # always
use warnings; # always
use ARS;
use Math::Rand48; # 48-bit RNG
use Math::Round;

# WARNING - THIS CODE WILL NOT RUN BECAUSE ARS.PM 
# CONNECTS TO NMDP DATABASE
# CODE WILL BE PORTED TO PYTHON TO USE PYARD MODULE

# initialize ARS
use Connect;
my $dbh = &Connect::Connect;
my $ARS_reduce = 1;
my $ars = new ARS("nmdp","g","3.4.0",1);

# command line parameters
my $exp_name = shift @ARGV; # name of experiment
my $fully_typed = shift @ARGV; # will randomly mask C and DQB1 if set to 1
my $vary_family_size = shift @ARGV; # will vary between 2 and 5 if set to 1
my $resolution = shift @ARGV; # only set to DNA2 - can't use simulated typing files
my $recomb_rate = shift @ARGV; # set to 0.0 to have no recombination

# fix different random seed based on experminent
if ($exp_name eq "EASY") {
  srand("20190721");
}
elsif ($exp_name eq "MM") {
  srand("20200721");
}
elsif ($exp_name eq "LESSEASY") {
  srand("20210721");
}
elsif ($exp_name eq "HARD") {
  srand("20220721");
}
elsif ($exp_name eq "HARDER") {
  srand("20230721");
}
elsif ($exp_name eq "HARDEST") {
  srand("20240721");
}
else {
  srand("20250721");
}

# fixed default parameters
my $typed_c = 0.25;
my $typed_dqb1 = 0.1;
my $nfamilies = 500;
my $nchildren_per_family = 5;
# my $nchildren = $nfamilies * $nchildren_per_family;
my $race = "CAU";

# load high res to DNA 2-digit map
my (%dna2_map);
&loadAlleleFamilyMap(\%dna2_map);

# load parental and maternal haplotypes from family file
my %family_haplo1_paternal;
my %family_haplo2_paternal;
my %family_haplo1_maternal;
my %family_haplo2_maternal;
my %family_paternal_id;
my %family_maternal_id;
my $file_family_master = "FAMILY_HAPLO_MASTER_" . $exp_name . ".csv";
open (FAMILY,"$file_family_master") || die "File missing: $file_family_master\n";
while (<FAMILY>) {
  chomp;
  my ($family_id,$paternal_id,$maternal_id,$haplo1_paternal,$haplo2_paternal,$haplo1_maternal,$haplo2_maternal) = split /,/,$_;
  if ($family_id eq "Family_ID") { next; } # skip header row
  $family_haplo1_paternal{$family_id} = $haplo1_paternal;
  $family_haplo2_paternal{$family_id} = $haplo2_paternal;
  $family_haplo1_maternal{$family_id} = $haplo1_maternal;
  $family_haplo2_maternal{$family_id} = $haplo2_maternal;
  $family_paternal_id{$family_id} = $paternal_id;
  $family_maternal_id{$family_id} = $maternal_id;
}
close FAMILY;

# generate child haplo pairs by randomly selecting haplotype pairs from parents

# open child haplotype master file
my $file_child_master = "CHILD_HAPLO_MASTER_" . $exp_name . ".csv";
open (CHILD,">$file_child_master") || die "Can't open file: $file_child_master\n";

my $file_child_workshop = "CHILD_HAPLO_WORKSHOP_" . $exp_name . ".csv";
open (CHILD_WORKSHOP,">$file_child_workshop") || die "Can't open file: $file_child_workshop\n";

my $file_recomb_master = "RECOMB_HAPLO_MASTER_" . $exp_name . ".csv";
open (RECOMB,">$file_recomb_master") || die "Can't open file: $file_recomb_master\n";

print RECOMB "INDIVIDUAL_ID,RECOMB_PAT,RECOMB_MAT\n";

my $individual_id = "100000";
foreach my $family_id (sort keys %family_paternal_id) {

  # vary number of children in the family between 2 and 5 if selected
  if ($vary_family_size == 1) {
    $nchildren_per_family = 5 * drand48();
    if ($nchildren_per_family < 2) {
      $nchildren_per_family = 2;
    }
  }

  for (my $i=0;$i<$nchildren_per_family;$i++) {


    my $haplo1_paternal = $family_haplo1_paternal{$family_id};
    my $haplo2_paternal = $family_haplo2_paternal{$family_id};
    my $haplo1_maternal = $family_haplo1_maternal{$family_id};
    my $haplo2_maternal = $family_haplo2_maternal{$family_id};
    my $recombined_paternal = 0; # 1 if recombined
    my $recombined_maternal = 0; # 1 if recombined
    my $recombined_loc_pair_pat = "NA"; # NA if not recombined
    my $recombined_loc_pair_mat = "NA"; # NA if not recombined

    # chance at recombination of paternal and maternal haplotypes
    ($haplo1_paternal,$haplo2_paternal,$recombined_paternal,$recombined_loc_pair_pat) = 
      recombine_haplos($haplo1_paternal,$haplo2_paternal,$recomb_rate);
    ($haplo1_maternal,$haplo2_maternal,$recombined_maternal,$recombined_loc_pair_mat) = 
      recombine_haplos($haplo1_maternal,$haplo2_maternal,$recomb_rate);
    
    # print STDERR "LOC PAIR $recombined_loc_pair_pat $recombined_loc_pair_mat\n";

    my $paternal_id = $family_paternal_id{$family_id};
    my $maternal_id = $family_maternal_id{$family_id};
    $individual_id++;

    my $child_paternal_haplo;
    if (drand48() < 0.5) {
      $child_paternal_haplo = $haplo1_paternal;
    }
    else {
      $child_paternal_haplo = $haplo2_paternal;
    }

    my $child_maternal_haplo;
    if (drand48() < 0.5) {
      $child_maternal_haplo = $haplo1_maternal;
    }
    else {
      $child_maternal_haplo = $haplo2_maternal;
    }

    my $sex;
    if (drand48() < 0.5) {
      $sex = 0;
    }
    else {
      $sex = 1;
    }

    # print STDERR "$child_paternal_haplo $child_maternal_haplo\n";

    my ($a1,$c1,$b1,$drb1_1,$dqb1_1) = split /~/,$child_paternal_haplo;
    my ($a2,$c2,$b2,$drb1_2,$dqb1_2) = split /~/,$child_maternal_haplo;

    my $pheno_a = sort_pheno($a1,$a2);
    my $pheno_c = sort_pheno($c1,$c2);
    my $pheno_b = sort_pheno($b1,$b2);
    my $pheno_drb1 = sort_pheno($drb1_1,$drb1_2);
    my $pheno_dqb1 = sort_pheno($dqb1_1,$dqb1_2);

    if ($resolution eq "DNA2") {

      # make expanded GL Strings

      my ($a1_loc,$a1_typ) = split /\*/,$a1;
      my ($a2_loc,$a2_typ) = split /\*/,$a2;
      my ($c1_loc,$c1_typ) = split /\*/,$c1;
      my ($c2_loc,$c2_typ) = split /\*/,$c2;
      my ($b1_loc,$b1_typ) = split /\*/,$b1;
      my ($b2_loc,$b2_typ) = split /\*/,$b2;
      my ($dr1_loc,$dr1_typ) = split /\*/,$drb1_1;
      my ($dr2_loc,$dr2_typ) = split /\*/,$drb1_2;
      my ($dq1_loc,$dq1_typ) = split /\*/,$dqb1_1;
      my ($dq2_loc,$dq2_typ) = split /\*/,$dqb1_2;

      my $a1_2dig = substr($a1_typ,0,2);
      my $a2_2dig = substr($a2_typ,0,2);
      my $c1_2dig = substr($c1_typ,0,2);
      my $c2_2dig = substr($c2_typ,0,2);
      my $b1_2dig = substr($b1_typ,0,2);
      my $b2_2dig = substr($b2_typ,0,2);
      my $dr1_2dig = substr($dr1_typ,0,2);
      my $dr2_2dig = substr($dr2_typ,0,2);
      my $dq1_2dig = substr($dq1_typ,0,2);
      my $dq2_2dig = substr($dq2_typ,0,2);
      
      my $a1_glstring = $dna2_map{$a1_loc}{$a1_2dig};
      my $a2_glstring = $dna2_map{$a2_loc}{$a2_2dig};
      my $c1_glstring = $dna2_map{$c1_loc}{$c1_2dig};
      my $c2_glstring = $dna2_map{$c2_loc}{$c2_2dig};
      my $b1_glstring = $dna2_map{$b1_loc}{$b1_2dig};
      my $b2_glstring = $dna2_map{$b2_loc}{$b2_2dig};
      my $dr1_glstring = $dna2_map{$dr1_loc}{$dr1_2dig};
      my $dr2_glstring = $dna2_map{$dr2_loc}{$dr2_2dig};
      my $dq1_glstring = $dna2_map{$dq1_loc}{$dq1_2dig};
      my $dq2_glstring = $dna2_map{$dq2_loc}{$dq2_2dig};

      $pheno_a = $a1_glstring . "+" . $a2_glstring;
      $pheno_c = $c1_glstring . "+" . $c2_glstring;
      $pheno_b = $b1_glstring . "+" . $b2_glstring;
      $pheno_drb1 = $dr1_glstring . "+" . $dr2_glstring;
      $pheno_dqb1 = $dq1_glstring . "+" .  $dq2_glstring;
    }


    # mask C or DQB1 typing
    if ($fully_typed == 0) {
      my $rand_c = drand48();
      if ($rand_c > $typed_c) {
        $pheno_c = "C*UUUU+C*UUUU";
      }

      my $rand_dqb1 = drand48();
      if ($rand_dqb1 > $typed_dqb1) {
        $pheno_dqb1 = "DQB1*UUUU+DQB1*UUUU";
      }
    }

    my $pheno = join("^",$pheno_a,$pheno_c,$pheno_b,$pheno_drb1,$pheno_dqb1);

    my $child_out = join (",",$family_id,$individual_id,$paternal_id,$maternal_id,$sex,$pheno,$child_paternal_haplo,$race,$child_maternal_haplo,$race);
    my $child_workshop_out = join (",",$family_id,$individual_id,$paternal_id,$maternal_id,$sex,$pheno);
    my $recomb_out = join (",",$individual_id,$recombined_loc_pair_pat,$recombined_loc_pair_mat);
    print CHILD "$child_out\n";
    print CHILD_WORKSHOP "$child_workshop_out\n";
    print RECOMB "$recomb_out\n";
  }
}

close CHILD;
close CHILD_WORKSHOP;
close RECOMB;

##############################################################################
# Function: sort_pheno - sort HLA alleles in phenotype
##############################################################################
sub sort_pheno {
   my ($typ1,$typ2) = @_;
   if ($typ1 ge $typ2) {
      ($typ1, $typ2) = ($typ2, $typ1);
   }
   return "$typ1+$typ2";
}

##############################################################################
# Function: recombine_haplos - sort HLA alleles in phenotype
##############################################################################
sub recombine_haplos {
  my ($hap1,$hap2,$recomb_rate) = @_;

  my ($a1,$c1,$b1,$drb1_1,$dqb1_1) = split /~/,$hap1;
  my ($a2,$c2,$b2,$drb1_2,$dqb1_2) = split /~/,$hap2;

  # pick random number to determine if recombination will occur
  my $recomb_loc_pair = "NA";
  my $rand_recomb = drand48();
  if ($rand_recomb > $recomb_rate) {
    my $recombined = 0;
      # return non-recombined haplotypes
    return ($hap1, $hap2, $recombined, $recomb_loc_pair);
  }

  # random number between 0 and 3 
  my $rand_loc = int(4 * drand48());


  if ($rand_loc == 0) { # between A and C
    ($a1, $a2) = ($a2, $a1);
    $recomb_loc_pair = "A~C";
  }
  elsif ($rand_loc == 1) { # between C and B
    ($a1, $a2) = ($a2, $a1);
    ($c1, $c2) = ($c2, $c1);
    $recomb_loc_pair = "C~B";
  }
  elsif ($rand_loc == 2) { # between B and DRB1
    ($drb1_1, $drb1_2) = ($drb1_2, $drb1_1);
    ($dqb1_1, $dqb1_2) = ($dqb1_2, $dqb1_1);
    $recomb_loc_pair = "B~DRB1";
  }
  elsif ($rand_loc == 3) { # between DRB1 and DQB1
    ($dqb1_1, $dqb1_2) = ($dqb1_2, $dqb1_1);
    $recomb_loc_pair = "DRB1~DQB1";
  }

  my $hap1_recomb = join("~",$a1,$c1,$b1,$drb1_1,$dqb1_1);
  my $hap2_recomb = join("~",$a2,$c2,$b2,$drb1_2,$dqb1_2);
  my $recombined = 1;

  # print STDERR "Recombined $rand_recomb $rand_loc $hap1 $hap2 $hap1_recomb $hap2_recomb $recombined $recomb_loc_pair\n";

  return ($hap1_recomb,$hap2_recomb,$recombined,$recomb_loc_pair);
}



##############################################################################
# Function: loadAlleleFamilyMap - load 2-digit DNA mappings
##############################################################################
sub loadAlleleFamilyMap {

    my ($rdna2_map) = @_;

open (DNA,"rel_dna_ser.txt") || die "rel_dna_ser.txt";

my %family;
while (<DNA>) {
  chomp;
  my ($loc,$allele,$sero,$sero_possible,$sero_assumed,$expert_exceptions) = split /\;/, $_;

  $loc =~ tr/*//d; # trim "*" off locus name
 
  if (($loc ne "A") && ($loc ne "C") && ($loc ne "B") && ($loc ne "DRB1") && ($loc ne "DQB1")) { next; }

  my ($family,$protein,$synony,$noncoding) = split /:/,$allele;
  $allele = "$family:$protein";
  if (defined $noncoding) {
    if ($noncoding =~ m/N/) { $allele .= "N"; }
    if ($noncoding =~ m/Q/) { $allele .= "Q"; }
    if ($noncoding =~ m/L/) { $allele .= "L"; }
  }
  if (defined $synony) {
    if ($synony =~ m/N/) { $allele .= "N"; }
    if ($synony =~ m/Q/) { $allele .= "Q"; }
    if ($synony =~ m/L/) { $allele .= "L"; }
  }
  push @{$family{"$loc*$family"}},"$loc*$allele";

  if (substr ($loc,0,1) eq "#") { next; }

}

close DNA;

foreach my $family (sort keys %family) {

  my %ars_dupe_check;
  my @alleles_family_ars;
  my (@alleles_family) = @{$family{$family}};

  foreach my $allele (@alleles_family) {
    my $allele_ars;
    if ($ARS_reduce == 1) {
      # $allele_ars = Identical::getIdenticalg($allele);
      $allele_ars = $ars->redux($allele);
    }
    else {
      $allele_ars = "$allele";
    }
    if (!exists $ars_dupe_check{$allele_ars}) {
      push @alleles_family_ars,$allele_ars;
      $ars_dupe_check{$allele_ars} = 1;
    }
  }
  my $glstring_ars = join "/",@alleles_family_ars;

  my ($loc,$family_typ) = split /\*/,$family;

  # print STDERR "$loc $family_typ $glstring_ars\n";

  $$rdna2_map{$loc}{$family_typ} = $glstring_ars;
}


} # end loadAlleleFamilyMap



exit 0;