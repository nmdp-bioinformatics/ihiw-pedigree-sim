# ihiw-pedigree-sim
Simulation of ambiguous HLA typing data inherited within family pedigrees.

This repository includes the dataset, instructions, and code for the 18th International Histocompatibility and Immunogenetics Workshop project on development of haplotype and pedigree analyses software tools.

https://www.ihiw18.org/component-bio-informatics/project-kir-haplotypes/





## Pedigree Simulation Files to Participants to Run Through Pedigree Analysis Tools:



#### Box Link to Simulation Files:

https://tulane.box.com/s/8si3cvuekjkhaja7o2jckgcilpmt8v0x



#### WORKSHOP data file for distribution to participants:

Input file for pedigree tool for cross-validation study is distributed to IHIW pedigree tools project participants

Same as PEDIGREE master table except the last 4 columns with true haplotypes are removed to blind the project participants.

`Family_ID,Individual_ID,Paternal_ID,Maternal_ID,Sex,Phenotype`



#### Modified PED File Format:

File format was inspired by PED format (https://software.broadinstitute.org/gatk/documentation/article?id=11016):

```
Family_ID
Individual_ID
Paternal_ID
Maternal_ID
Sex (1=male; 2=female; other=unknown) - always use other - doesn’t matter - remove
Phenotype - GL String of Ambiguous HLA typing
```



Validation scripts use the `PEDIGREE` and `RECOMB` files to determine if predictions were correct. These have additional columns that are extensions to format to show the true phased haplotypes and population sampled:

```
Haplotype1
Race1
Haplotype2
Race2
```



#### Suggested Format for Pedigree Tool Output File to Submit to IHIW Workshop:

`Individual_ID,Hap1_Inferred,Hap1_Parent,Hap1_Recombined,Hap2_Inferred,Hap2_Parent,Hap2_Recombined`

Assign pedigree-inferred imputed haplotype and parent for each child. 

Parents should have "NA" for `Hap1_Parent` and `Hap2_Parent`

`Hap1_Recombined` should provide pair of loci between which the recombination event occurred (e.g. `A~C` or `B~DRB1`), or `NA` if no recombination event occurred.



## Pedigree Simulation Files to Cross-Validate Performance of Pedigree Tools:



The following files are described for information purposed. These files contain the true sampled haplotypes that pedigree tools aim to recover. These files will be used for cross-validation of pedigree tool performance, but are not provided to project participants.



#### PARENT haplotype master table: 

Contains unphased phenotype for each individual and true sampled haplotypes.

`Individual_ID` are sequentially numbered across parent and children.

`Family_ID` - Links mother and father as well as the children

`Paternal_ID` and `Maternal_ID` are both blank for parents. They are the parents, they don't have parents.

Example data for two pairs of parents from `PARENT_HAPLO_MASTER_EASY.csv`:

```
Family_ID,Individual_ID,Paternal_ID,Maternal_ID,Sex,Phenotype,Hap1,Race1,Hap2,Race2
000001,000001,NULL,NULL,1,A*02:01+A*30:04^C*03:04+C*14:02^B*40:01+B*51:22^DRB1*04:01+DRB1*14:02^DQB1*03:02+DQB1*03:02,A*02:01g~C*03:04g~B*40:01g~DRB1*04:01~DQB1*03:02g,CAU,A*30:04~C*14:02g~B*51:22~DRB1*14:02~DQB1*03:02g,CAU
000001,000002,NULL,NULL,2,A*24:02+A*02:01^C*12:03+C*06:02^B*38:01+B*57:01^DRB1*13:01+DRB1*07:01^DQB1*06:03+DQB1*03:03,A*24:02g~C*12:03g~B*38:01~DRB1*13:01~DQB1*06:03g,CAU,A*02:01g~C*06:02g~B*57:01g~DRB1*07:01g~DQB1*03:03g,CAU
000002,000003,NULL,NULL,1,A*11:01+A*11:01^C*04:01+C*04:01^B*44:03+B*35:01^DRB1*15:01+DRB1*01:01^DQB1*06:02+DQB1*05:01,A*11:01g~C*04:01g~B*44:03g~DRB1*15:01g~DQB1*06:02,CAU,A*11:01g~C*04:01g~B*35:01g~DRB1*01:01~DQB1*05:01g,CAU
000002,000004,NULL,NULL,2,A*03:01+A*02:01^C*07:02+C*06:02^B*07:02+B*13:02^DRB1*15:01+DRB1*07:01^DQB1*06:02+DQB1*02:01,A*03:01g~C*07:02g~B*07:02g~DRB1*15:01g~DQB1*06:02,CAU,A*02:01g~C*06:02g~B*13:02g~DRB1*07:01g~DQB1*02:01g,CAU
```



#### CHILD haplotype master table:

`Individual_ID` for children are sequentially numbered from 100001 to avoid collision with parents.

Example data for one family from `CHILD_HAPLO_MASTER_EASY.csv`:

```
Family_ID,Individual_ID,Paternal_ID,Maternal_ID,Sex,Phenotype,Hap1,Race1,Hap2,Race2
000001,100001,000001,000002,0,A*02:01g+A*24:02g^C*03:04g+C*12:03g^B*38:01+B*40:01g^DRB1*04:01+DRB1*13:01^DQB1*03:02g+DQB1*06:03g,A*02:01g~C*03:04g~B*40:01g~DRB1*04:01~DQB1*03:02g,CAU,A*24:02g~C*12:03g~B*38:01~DRB1*13:01~DQB1*06:03g,CAU
000001,100002,000001,000002,0,A*02:01g+A*02:01g^C*03:04g+C*06:02g^B*40:01g+B*57:01g^DRB1*04:01+DRB1*07:01g^DQB1*03:02g+DQB1*03:03g,A*02:01g~C*03:04g~B*40:01g~DRB1*04:01~DQB1*03:02g,CAU,A*02:01g~C*06:02g~B*57:01g~DRB1*07:01g~DQB1*03:03g,CAU
000001,100003,000001,000002,1,A*02:01g+A*30:04^C*06:02g+C*14:02g^B*51:22+B*57:01g^DRB1*07:01g+DRB1*14:02^DQB1*03:02g+DQB1*03:03g,A*30:04~C*14:02g~B*51:22~DRB1*14:02~DQB1*03:02g,CAU,A*02:01g~C*06:02g~B*57:01g~DRB1*07:01g~DQB1*03:03g,CAU
000001,100004,000001,000002,0,A*02:01g+A*24:02g^C*03:04g+C*12:03g^B*38:01+B*40:01g^DRB1*04:01+DRB1*13:01^DQB1*03:02g+DQB1*06:03g,A*02:01g~C*03:04g~B*40:01g~DRB1*04:01~DQB1*03:02g,CAU,A*24:02g~C*12:03g~B*38:01~DRB1*13:01~DQB1*06:03g,CAU
000001,100005,000001,000002,1,A*02:01g+A*30:04^C*06:02g+C*14:02g^B*51:22+B*57:01g^DRB1*07:01g+DRB1*14:02^DQB1*03:02g+DQB1*03:03g,A*30:04~C*14:02g~B*51:22~DRB1*14:02~DQB1*03:02g,CAU,A*02:01g~C*06:02g~B*57:01g~DRB1*07:01g~DQB1*03:03g,CAU
```

This file is actually missing the header because it is concatenated at end of parent file.

Sex is randomized for children. Sex of children should not have an impact on the pedigree assignment.



#### FAMILY haplotype master table:

Row for each family instead of for each individual.

Shows true haplotype pairs of each father and mother pair.

```
Family_ID,Paternal_ID,Maternal_ID,Haplo_P1,Haplo_P2,Haplo_M1,Haplo_M2
000001,000001,000002,A*02:01g~C*03:04g~B*40:01g~DRB1*04:01~DQB1*03:02g,A*30:04~C*14:02g~B*51:22~DRB1*14:02~DQB1*03:02g,A*24:02g~C*12:03g~B*38:01~DRB1*13:01~DQB1*06:03g,A*02:01g~C*06:02g~B*57:01g~DRB1*07:01g~DQB1*03:03g
000002,000003,000004,A*11:01g~C*04:01g~B*44:03g~DRB1*15:01g~DQB1*06:02,A*11:01g~C*04:01g~B*35:01g~DRB1*01:01~DQB1*05:01g,A*03:01g~C*07:02g~B*07:02g~DRB1*15:01g~DQB1*06:02,A*02:01g~C*06:02g~B*13:02g~DRB1*07:01g~DQB1*02:01g
000003,000005,000006,A*01:01g~C*12:03g~B*18:01g~DRB1*11:04~DQB1*03:01g,A*02:01g~C*03:03g~B*15:01g~DRB1*15:01g~DQB1*06:02,A*24:02g~C*06:02g~B*13:02g~DRB1*07:01g~DQB1*02:01g,A*03:01g~C*06:02g~B*47:01g~DRB1*13:01~DQB1*06:03g
000004,000007,000008,A*02:01g~C*02:02g~B*27:05g~DRB1*08:01g~DQB1*04:02,A*68:01g~C*03:04g~B*40:01g~DRB1*04:01~DQB1*03:01g,A*03:01g~C*07:02g~B*07:02g~DRB1*13:01~DQB1*06:03g,A*02:01g~C*07:01g~B*08:01g~DRB1*04:01~DQB1*03:02g
```

For cross-validation study, this file will not be shared.



#### PEDIGREE master table:

Generated by concatenating PARENT and CHILD files together.

Contains true haplotype pairs for every individual.

`Family_ID,Individual_ID,Paternal_ID,Maternal_ID,Sex,Phenotype,Hap1,Race1,Hap2,Race2`



#### RECOMB table:

For each child, this table records where the recombination occurred for haplotypes inherited from paternal and maternal side.

If no recombination occurred, "NA" is recorded.

```
INDIVIDUAL_ID,RECOMB_PAT,RECOMB_MAT
100001,A~C,NA
100002,NA,NA
100003,C~B,NA
100004,NA,NA
```





## Pedigree Tool Objectives:

#### 1. Recover true haplotype phase for each individual

Validation procedure - for each row in `PEDIGREE` master file, check that inferred haplotypes matches the true haplotypes `Hap1` and `Hap2`



#### 2. Determine which haplotypes the children inherited from parents

Validation procedure - For each child in `PEDIGREE` master file, assign each haplotype as coming from father or mother.

Check true parental family haplotypes in `FAMILY` file - `Haplo_P1`, `Haplo_P2`, `Haplo_M1`, and `Haplo_M2`



#### 3. Determine where any recombination events occurred between parental haplotypes before they inherited in particular children.

Because recombination during meiosis would not change the HLA genotype of the children, we had to model changes in haplotypes before inheritance from the parents. 

Validation procedure - Recombination events can be determined by lack of overlap between true non-recombined parental haplotypes in `FAMILY` file (`Haplo_P1`, `Haplo_P2`, `Haplo_M1`, and `Haplo_M2`) and true potentially-recombined child haplotypes in `PEDIGREE` master file (`Hap1` and `Hap2` )

For children that cannot be assigned an inherited haplotype from a parent, determine where recombination event occurred between `Haplo_P1` and `Haplo_P2` and/or  `Haplo_M1` and `Haplo_M2`. 



#### 4. FUTURE - Assign race/ethnic labels per haplotype (when multi-race individuals are introduced)



#### 5. FUTURE - Determine false paternity (when that feature is introduced)



#### 6. FUTURE - Determine parental allele from which mutation or gene conversion occurred (when that feature is introduced)





 

## Pedigree Typing Simulation Algorithm



#### 1. Simulate Parents by Sampling from Simulated HLA Typing Files

Haplotype phase is hidden in `Phenotype` field.

Missing / ambiguous typing at some loci can be introduced in `Phenotype`.  Ambiguity string can be found in typing simulation files input to pedigree simulator.



#### 2. Simulate Children by Sampling 1 Haplotype From Each Parent

Number of children per family is a parameter that can vary. Default is 5.

Recombination can optionally occur at most once per haplotype. Two loci between which the recombination occurs is randomly chosen.

Sex of children is randomized and should not impact pedigree analysis.

Can't easily simulate SSO typing for children because would have to run Java imputation server for batch of new haplotype combinations in children.

Mutation is not yet modeled.




## How were the simulations created? Simulated HLA Typing Data:

Simulated HLA typing files are generated by sampling from published NMDP full registry haplotype distributions for US populations. DOI: 10.1016/j.humimm.2013.06.025

The first step was generating 5-locus simulated high resolution genotypes for the Caucasian population. High resolution is defined as two-field IMGT/HLA allele nomenclature, combining alleles with identical nucleotide sequence in the antigen recognition domain (ARD)

True phased haplotypes that were sampled are retained. 

File format examples:

```
ID%Pheno%Hap1%Hap2

D00000001%A*02:01g+A*30:04^C*03:04g+C*14:02g^B*40:01g+B*51:22^DRB1*04:01+DRB1*14:02^DQB1*03:02g+DQB1*03:02g;A*02:01g+A*30:04^C*03:04g+C*14:02g^B*40:01g+B*51:22^DRB1*04:01+DRB1*14:02^DQB1*03:02g+DQB1*03:02g%A*02:01g~C*03:04g~B*40:01g~DRB1*04:01~DQB1*03:02g%A*30:04~C*14:02g~B*51:22~DRB1*14:02~DQB1*03:02g
```



## Pedigree Simulation Scripts:

WARNING: Perl scripts are provided for information only, but cannot be run to ARS.pm dependency on NMDP database. 

Needs to be ported to Python and switch to py-ARD to eliminate ARS.pm dependency on NMDP database. 




## Pedigree Simulation Experiment Parameters

Parameters for parent simulation with `pedigree_sim_parents.pl`:

`exp_name` - name of experiment

`fully_typed` - will randomly mask C and DQB1 if set to 1 - Probability of missing C and DQB1 is 0.25 and 0.1.

`vary_nparents_output` - number of parents will randomly vary between 0 and 2 if set to 1, otherwise will always have 2 parents

`resolution` - ambiguous genotypes for `PARENT` file - can be high resolution (HR),  first-field DNA (DNA2), or SSO.



Parameters for child simulation with `pedigree_sim_children.pl`:

`exp_name` - name of experiment

`fully_typed` - will randomly mask C and DQB1 if set to 1 - Probability of missing C and DQB1 is 0.25 and 0.1.

`vary_family_size` - number of children in each family will randomly vary between 2 and 5 if set to 1, otherwise will have 5 children per family

`resolution` - either high resolution (HR) or can be set to first-field DNA (DNA2) - can't use simulated typing files to run SSO

`recomb_rate` - set to 0.0 to have no recombination


Script to generate simulated pedigree experiments from simulated HLA typings:

` bash run_pedigree_sim.sh` 



#### List of experiments:

`EASY` - Fully Typed, 2 parents, 5 children, No ambiguity, No recombination

`MM` - Fully Typed, 2 parents, 5 children, No ambiguity, 5% recombination

`LESSEASY` - Fully Typed, 2 parents, 5 children, Typing ambiguity, No recombination

`HARD` -  Fully Typed, No parents, 5 children, Typing ambiguity, No recombination

`HARDER` - Fully Typed, No parents, 5 children, Typing ambiguity, 5% recombination

`HARDEST` - Some Missing C/DQB1 Typing, 0-2 parents, 2-5 children, Typing ambiguity, 5% recombination



```
perl pedigree_sim_parents.pl EASY 1 0 HR
perl pedigree_sim_children.pl EASY 1 0 HR 0.0

perl pedigree_sim_parents.pl MM 1 0 HR
perl pedigree_sim_children.pl MM 1 0 HR 0.05

perl pedigree_sim_parents.pl LESSEASY 1 0 SSO
perl pedigree_sim_children.pl LESSEASY 1 0 DNA2 0.00

perl pedigree_sim_parents.pl HARD 1 0 SSO
perl pedigree_sim_children.pl HARD 1 0 DNA2 0.00

perl pedigree_sim_parents.pl HARDER 1 0 SSO
perl pedigree_sim_children.pl HARDER 1 0 DNA2 0.05

perl pedigree_sim_parents.pl HARDEST 0 1 SSO
perl pedigree_sim_children.pl HARDEST 0 1 DNA2 0.05
```

