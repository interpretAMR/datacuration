# AMRrules data curation

This repo will be a space for sharing R code to process & analyse matched AMR genotype & phenotype data, to support the ESGEM-AMR working groups to [develop AMRrules](https://github.com/interpretAMR/AMRrulesCuration) for the interpretation of AMR genotype data in terms of expected phenotypes.

## Example code

R code with examples using public NCBI AST data, and AllTheBacteria AMRfinderplus results, to explore core genes and the association of AMR genotypes with phenotypes.

* `coreGenes.Rmd` - identify AMR gene nodes that have high within-species frequency (alone or combined with its child nodes), so need core gene interpretation rules
* `Enterobacter/Enterobacter.Rmd` (output in `Enterobacter/Enterobacter.html`)
* `Pseudomonas/Pseudomonas_aeruginosa.Rmd` (output in `Pseudomonas/Pseudomonas_aeruginosa.html`)

Helper functions to parse and analysis AST/AMRfp data:
* `AllTheBacteria_functions.R` (function to extract AMRfinderplus results for a given species)
* `NARMS_functions.R` (function to parse AST data from the US National Antimicrobial Resistance Monitoring System, [NARMS](https://www.fda.gov/animal-veterinary/national-antimicrobial-resistance-monitoring-system/narms-now-integrated-data)
* `functions.R` (functions for comparing AST vs AMRfp genotype data, including assessing solo positive-predictive value per marker, and fitting and plotting logistic regression for a given drug and associated markers)

## Data files

### Antimicrobial susceptibility testing (AST) phenotypes

`/AST`

public AST data downloaded from [NCBI AST browser](https://www.ncbi.nlm.nih.gov/pathogens/ast#scientific_name:Escherichia%20coli)

* `Ecoli_AST.tsv.gz` -  (9094 unique biosamples with AST data on at least one drug (median 15 drugs, IQR 14-18 drugs)) (<2 MB)
* `AST_Enterobacter.tsv.gz` (221 unique biosamples with AST data on at least one drug (median 21 drugs, IQR 15-26 drugs) (<100 KB)
* `AST_Pseudomonas_aeruginosa.tsv.gz` (607 unique biosamples with AST data on at least one drug (median 1 drugs, IQR 6-14 drugs) (<100 KB)

### NCBI refgenes DB
* `refgenes_20241003.tsv` - Local download of NCBI [refgene DB](https://www.ncbi.nlm.nih.gov/pathogens/refgene/) of AMR determinants (<2 MB)
* `ReferenceGeneHierarchy.txt` - Local download of NCBI [node hierarchy](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneHierarchy.txt) (1 MB)

### AMRFinderPlus results - via AllTheBacteria
AMRFinderPlus (v3.12.8, DB 2024-01-31.1) output for all samples in [AllTheBacteria](https://github.com/AllTheBacteria/AllTheBacteria/tree/main/reproducibility/All-samples/AMR/AMRFinderPlus) can be downloaded from here: [https://osf.io/zgexh](https://osf.io/zgexh) (0.5GB .tsv.gz file, expands to 6.2GB)

This has >26 M lines, one per genome-gene hit, but as it's the raw AMRfinderplus output it doesn't have any column to indicate the species for each row.

Species calls were made separately for AllTheBacteria using sylph, and can be downloaded from [here](https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/metadata/species_calls.tsv.gz) (13.4 MB .tsv.gz file)

This repo contains files with the AMRfinderplus output extracted for all the ESGEM-AMR target organism subgroups:

`/AllTheBacteria`

**Organism, lines, strains**
* Achromobacter 1606 440
* Acinetobacter baumannii 244302 14460
* Aeromonas 9954 1837
* Bordetella 715 519
* Brucella 12515 2514
* Burkholderia cepacia 1061 238
* Burkholderia mallei 20340 5119  [Note, this is how GTDB labels Burkholderia pseudomallei]
* Campylobacter 314666 107848
* Chryseobacterium 107 74
* Clostridioides difficile 373507 26347
* Corynebacterium diphtheriae 2373 822
* Escherichia coli 9179705 314978 [Note, only selected columns are included due to file size restrictions on github, parse your own copy if you want the complete file]
* Edwardsiella 97 26
* Enterobacter 334367 9545
* Enterococcus 551141 38920
* Haemophilus influenzae 16867 10990
* Klebsiella pneumoniae 2033309 57070
* Klebsiella quasipneumoniae 87992 3034
* Klebsiella variicola 48994 2423
* Legionella 1986 1803
* Listeria 236353 61136
* Mycobacterium tuberculosis 402953 134246
* Mycoplasma 280 173
* Neisseria gonorrhoeae 667163 43643
* Neisseria meningitidis 198052 34561
* Pasteurella 751 122
* Proteus mirabilis 11032 1125
* Pseudomonas_aeruginosa 389574 25057
* Salmonella 6133567 534667 [Note, only selected columns are included due to file size restrictions on github, parse your own copy if you want the complete file]
* Serratia 44719 2534
* Shewanella 723 164
* Staphylococcus 3116441 119669
* Stenotrophomonas 23599 1998
* Streptococcus 328469 129068
* Treponema 9 8
* Vibrio 146485 19721
* Yersinia 65021 4370

Note that there are no AMRFP results for Ureaplasma
  
The script `ATB_ESGEM_orgs.Rmd` shows how we used R to pull out AMRfinderplus results for genomes belonging to a particular species
You can use it to extract data for other species, like this:

```
library(tidyverse)
library(dplyr)
source("AllTheBacteria_functions.R")

# set location of the full set of AMRfinderplus results (can be downloaded from here: [https://osf.io/zgexh](https://osf.io/zgexh) (0.5GB .tsv.gz file)
ATB_amrfp <- "AMRFP_results.tsv.gz" 

# name of organism to extract data for (we will search for this string in the organism names, so it can be a species or genus)
taxa <- "Enterobacter"

# extract AMRfinderplus results for this taxa of interest
amrfp <- atb_amrfp_filter_by_taxa(user_taxa=taxa, atb_armfp_results_path=ATB_amrfp)

# write out the resulting file
write_tsv(amrfp, file="Enterobacter_AFP.tsv.gz"))
```
