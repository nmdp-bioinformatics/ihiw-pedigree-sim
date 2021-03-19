#!/bin/bash

# Parents Params
# my $exp_name = shift @ARGV; # name of experiment
# my $fully_typed = shift @ARGV; # will randomly mask C and DQB1 if set to 1
# my $vary_nparents_output = shift @ARGV; # parents will vary between 0 and 2 if set to 1
# my $resolution = shift @ARGV; # only impacts parent file - not family file

# Children Params
# my $exp_name = shift @ARGV; # name of experiment
# my $fully_typed = shift @ARGV; # will randomly mask C and DQB1 if set to 1
# my $vary_family_size = shift @ARGV; # will vary between 2 and 5 if set to 1
# my $resolution = shift @ARGV; # only set to DNA2 - can't use simulated typing files
# my $recomb_rate = shift @ARGV; # set to 0.0 to have no recombination


# EASY - Fully Typed, 2 parents, 5 children, No ambiguity, No recombination

perl pedigree_sim_parents.pl EASY 1 0 HR
perl pedigree_sim_children.pl EASY 1 0 HR 0.0

cat PARENT_HAPLO_MASTER_EASY.csv CHILD_HAPLO_MASTER_EASY.csv > PEDIGREE_HAPLO_MASTER_EASY.csv
cat PARENT_HAPLO_WORKSHOP_EASY.csv CHILD_HAPLO_WORKSHOP_EASY.csv > PEDIGREE_HAPLO_WORKSHOP_EASY.csv

# MM - Fully Typed, 2 parents, 5 children, No ambiguity, 5% recombination

perl pedigree_sim_parents.pl MM 1 0 HR
perl pedigree_sim_children.pl MM 1 0 HR 0.05

cat PARENT_HAPLO_MASTER_MM.csv CHILD_HAPLO_MASTER_MM.csv > PEDIGREE_HAPLO_MASTER_MM.csv
cat PARENT_HAPLO_WORKSHOP_MM.csv CHILD_HAPLO_WORKSHOP_MM.csv > PEDIGREE_HAPLO_WORKSHOP_MM.csv

# LESSEASY - Fully Typed, 2 parents, 5 children, Typing ambiguity, No recombination

perl pedigree_sim_parents.pl LESSEASY 1 0 SSO
perl pedigree_sim_children.pl LESSEASY 1 0 DNA2 0.00

cat PARENT_HAPLO_MASTER_LESSEASY.csv CHILD_HAPLO_MASTER_LESSEASY.csv > PEDIGREE_HAPLO_MASTER_LESSEASY.csv
cat PARENT_HAPLO_WORKSHOP_LESSEASY.csv CHILD_HAPLO_WORKSHOP_LESSEASY.csv > PEDIGREE_HAPLO_WORKSHOP_LESSEASY.csv

# HARD - Fully Typed, No parents, 5 children, Typing ambiguity, No recombination

perl pedigree_sim_parents.pl HARD 1 0 SSO
perl pedigree_sim_children.pl HARD 1 0 DNA2 0.00

cp CHILD_HAPLO_MASTER_HARD.csv PEDIGREE_HAPLO_MASTER_HARD.csv
cp CHILD_HAPLO_WORKSHOP_HARD.csv PEDIGREE_HAPLO_WORKSHOP_HARD.csv

# HARDER - Fully Typed, No parents, 5 children, Typing ambiguity, 5% recombination

perl pedigree_sim_parents.pl HARDER 1 0 SSO
perl pedigree_sim_children.pl HARDER 1 0 DNA2 0.05

cp CHILD_HAPLO_MASTER_HARDER.csv PEDIGREE_HAPLO_MASTER_HARDER.csv
cp CHILD_HAPLO_WORKSHOP_HARDER.csv PEDIGREE_HAPLO_WORKSHOP_HARDER.csv

# HARDEST - Some missing C/DQB1 Typing, 0-2 parents, 2-5 children, Typing ambiguity, 5% recombination

perl pedigree_sim_parents.pl HARDEST 0 1 SSO
perl pedigree_sim_children.pl HARDEST 0 1 DNA2 0.05

cat PARENT_HAPLO_MASTER_HARDEST.csv CHILD_HAPLO_MASTER_HARDEST.csv > PEDIGREE_HAPLO_MASTER_HARDEST.csv
cat PARENT_HAPLO_WORKSHOP_HARDEST.csv CHILD_HAPLO_WORKSHOP_HARDEST.csv > PEDIGREE_HAPLO_WORKSHOP_HARDEST.csv