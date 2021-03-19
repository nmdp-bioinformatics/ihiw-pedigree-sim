#!/usr/bin/env perl
############################################################################
# SCRIPT NAME:  pedigree_sim_parents.pl
# DESCRIPTION:  simulate parental genotypes from HR typing file
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
use Math::Rand48; # 48-bit RNG
use Math::Round;


# WARNING - THIS CODE WILL NOT RUN BECAUSE ARS.PM 
# CONNECTS TO NMDP DATABASE
# CODE WILL BE PORTED TO PYTHON TO USE PYARD MODULE

my $impute_sim_path = "REDACTED";

# number of donors to select from simulated typings

# command line parameters
my $exp_name = shift @ARGV; # name of experiment
my $fully_typed = shift @ARGV; # will randomly mask C and DQB1 if set to 1
my $vary_nparents_output = shift @ARGV; # parents will vary between 0 and 2 if set to 1
my $resolution = shift @ARGV; # only impacts parent file - not family file

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
my $nparents = 2;
my $nfamilies = 500;
my $ndonor_hr = 2 * $nfamilies;
my $race = "CAU";

my $file_hr = $impute_sim_path . $resolution . "_FULL_PROTEIN_" . $race . ".out";
my @sampled_hr;
my @sampled;
open (HR,"$file_hr") || die "File missing: $file_hr\n";
while (<HR>) {
  chomp;
  my ($id,$pheno,$haplo1,$haplo2) = split /\%/,$_;
  my ($pheno_nomen,$pheno_gl) = split /\;/,$pheno;
  $pheno_nomen =~ tr/a-z//d;    # remove 'g' from ARS rollup
  my $sample_string = $id . "%" . $pheno_nomen . "%" . $haplo1 . "%" . $haplo2;
  push @sampled_hr,$sample_string;
}
close HR;

my $max_index = $ndonor_hr - 1;
while ($ndonor_hr >= 1) {
  my $index = Math::Round::nearest_floor(1,int (drand48() * $max_index));
  push @sampled,$sampled_hr[$index];
  # print STDERR "$sampled_hr[$index]\n";
  $sampled_hr[$index] = $sampled_hr[$max_index];
  $max_index--;
  $ndonor_hr--;
}

# open family haplotype master file
my $file_family_master = "FAMILY_HAPLO_MASTER_" . $exp_name . ".csv";
open (FAMILY,">$file_family_master") || die "Can't open file: $file_family_master\n";
print FAMILY "Family_ID,Paternal_ID,Maternal_ID,Haplo_P1,Haplo_P2,Haplo_M1,Haplo_M2\n";

# open parent haplotype master file
my $file_parent_master = "PARENT_HAPLO_MASTER_" . $exp_name . ".csv";
open (PARENT,">$file_parent_master") || die "Can't open file: $file_parent_master\n";
print PARENT "Family_ID,Individual_ID,Paternal_ID,Maternal_ID,Sex,Phenotype,Hap1,Race1,Hap2,Race2\n";

# open parent haplotype workshop file
my $file_parent_workshop = "PARENT_HAPLO_WORKSHOP_" . $exp_name . ".csv";
open (PARENT_WORKSHOP,">$file_parent_workshop") || die "Can't open file: $file_parent_workshop\n";
print PARENT_WORKSHOP "Family_ID,Individual_ID,Paternal_ID,Maternal_ID,Sex,Phenotype\n";


my $family_id = "000000";
my $individual_id = "000000";
my $sampled_index = 0;

# generate families one by one
for (my $i=0;$i<$nfamilies;$i++) {

  # randomly decide if parent should be output
  if ($vary_nparents_output == 1) {
    $nparents = int(3 * drand48());
  }

  $family_id++;

  # FATHER
  my $sampled_paternal_line = $sampled[$sampled_index];
  $sampled_index++;
  $individual_id++;
  my $paternal_id = $individual_id;
  my ($old_id_paternal,$pheno_paternal,$haplo1_paternal,$haplo2_paternal) = split /\%/,$sampled_paternal_line;

  # print STDERR "$pheno_paternal\n";

  # use phenotype not true haplotypes because typing can be ambiguous
  my ($pheno_a,$pheno_c,$pheno_b,$pheno_drb1,$pheno_dqb1) = split /\^/,$pheno_paternal;


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

  # print father to parent haplotype master file
  if ($nparents == 2) {
    my $paternal_out = join (",",$family_id,$individual_id,"NULL","NULL",1,$pheno,$haplo1_paternal,$race,$haplo2_paternal,$race);
    print PARENT "$paternal_out\n";
    my $paternal_workshop_out = join (",",$family_id,$individual_id,"NULL","NULL",1,$pheno);
    print PARENT_WORKSHOP "$paternal_workshop_out\n";    
  }

  # MOTHER
  my $sampled_maternal_line = $sampled[$sampled_index];
  $sampled_index++;
  $individual_id++;
  my $maternal_id = $individual_id;

  my ($old_id_maternal,$pheno_maternal,$haplo1_maternal,$haplo2_maternal) = split /\%/,$sampled_maternal_line;

  # use phenotype not true haplotypes because typing can be ambiguous
  ($pheno_a,$pheno_c,$pheno_b,$pheno_drb1,$pheno_dqb1) = split /\^/,$pheno_maternal;

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

  $pheno = join("^",$pheno_a,$pheno_c,$pheno_b,$pheno_drb1,$pheno_dqb1);

  # print mother to parent haplotype master file
  if ($nparents >= 1) {
    my $maternal_out = join (",",$family_id,$individual_id,"NULL","NULL",2,$pheno,$haplo1_maternal,$race,$haplo2_maternal,$race);
    print PARENT "$maternal_out\n";
    my $maternal_workshop_out = join (",",$family_id,$individual_id,"NULL","NULL",2,$pheno);
    print PARENT_WORKSHOP "$maternal_workshop_out\n";
  }

  # FAMILY
  # print parental haplotypes to family haplotype master file
  my $family_out = join (",",$family_id,$paternal_id,$maternal_id,$haplo1_paternal,$haplo2_paternal,$haplo1_maternal,$haplo2_maternal);
  print FAMILY "$family_out\n";

}

close PARENT;
close FAMILY;

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



exit 0;