# AMRrules data curation

This repo will be a space for sharing R code to process & analyse matched AMR genotype & phenotype data, to support the ESGEM-AMR working groups to [develop AMRrules](https://github.com/interpretAMR/AMRrulesCuration) for the interpretation of AMR genotype data in terms of expected phenotypes.

## Example code

R code with examples using public NCBI AST data, and AllTheBacteria AMRfinderplus results, to explore core genes and the association of AMR genotypes with phenotypes.

* `Enterobacter/Enterobacter.Rmd` (output in `Enterobacter/Enterobacter.html`)
* `Pseudomonas/Pseudomonas_aeruginosa.Rmd` (output in `Pseudomonas/Pseudomonas_aeruginosa.html`)

Helper functions to parse and analysis AST/AMRfp data:
* `AllTheBacteria_functions.R` (function to extract AMRfinderplus results for a given species)
* `NARMS_functions.R` (function to parse AST data from NARMS data)
* `functions.R` (functions for comparing AST vs AMRfp genotype data, including assessing solo positive-predictive value per marker, and fitting and plotting logistic regression for a given drug and associated markers)

## Data files

### Antimicrobial susceptibility testing (AST) phenotypes

`/AST`

* `Ecoli_AST.tsv.gz` - public AST data downloaded from [NCBI AST](https://www.ncbi.nlm.nih.gov/pathogens/ast#scientific_name:Escherichia%20coli) (9094 unique biosamples with AST data on at least one drug (median 15 drugs, IQR 14-18 drugs)) (<2 MB)
* `AST_Enterobacter.tsv.gz` (221 unique biosamples with AST data on at least one drug (median 21 drugs, IQR 15-26 drugs) (<100 KB)
* `AST_Pseudomonas_aeruginosa.tsv.gz` (607 unique biosamples with AST data on at least one drug (median 1 drugs, IQR 6-14 drugs) (<100 KB)

### NCBI refgenes DB
* `refgenes_20241003.tsv` - Local download of NCBI [refgene DB](https://www.ncbi.nlm.nih.gov/pathogens/refgene/) of AMR determinants (<2 MB)
* `ReferenceGeneHierarchy.txt` - Local download of NCBI [node hierarchy](https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneHierarchy.txt) (1 MB)

### AMRFinderPlus results - via AllTheBacteria
AMRFinderPlus (v3.12.8, DB 2024-01-31.1) output for all samples in [AllTheBacteria](https://github.com/AllTheBacteria/AllTheBacteria/tree/main/reproducibility/All-samples/AMR/AMRFinderPlus) can be downloaded from here: [https://osf.io/zgexh](https://osf.io/zgexh) (0.5GB .tsv.gz file, expands to 6.2GB)

This has >26 M lines, one per genome-gene hit, but as it's the raw AMRfinderplus output it doesn't have any column to indicate the species for each row.

Species calls were made separately for AllTheBacteria using sylph, and can be downloaded from [here](https://ftp.ebi.ac.uk/pub/databases/AllTheBacteria/Releases/0.2/metadata/species_calls.tsv.gz) (13.4 MB .tsv.gz file)

This repo contains files with the AMRfinderplus output extracted for selected species:

`/AllTheBacteria`

* _Escherichia coli_: `ATB_Ecoli_AFP.tsv.gz` - 9,179,705 lines for 314,978 unique genomes (43 MB)
* _Salmonella enterica_: `ATB_Salmonella_AFP.tsv.gz` - 6,133,567 lines for 534,667 unique genomes (26 MB)
* _Pseudomonas aeruginosa_: `ATB_Pseudomonas_aeruginosa_AFP.tsv.gz` - 389,574 lines for 24,854 unique genomes (7.6 MB)
* _Enterobacter_: `ATB_Enterobacter_AFP.tsv.gz` - 332,099 lines for 9,465 unique genomes (7.4 MB) 
  
The script `ATB.Rmd` shows how we used R to pull out AMRfinderplus results for genomes belonging to a particular species, which can be used to extract data for other species, like this:

```
# get list of E. coli samples
species_calls <- read_tsv("AllTheBacteria/ATB_species_calls.tsv.gz")
ecoli <- species_calls %>% filter(Species=="Escherichia coli") %>% pull(Sample)

# read in only those lines matching these samples
f <- function(x, pos) subset(x, Name %in% ecoli)
ecoli_AFP <- read_tsv_chunked("AMRFP_results.tsv.gz", DataFrameCallback$new(f), chunk_size = 10000)

# select key columns to keep output file size small-ish
ecoli_AFP %>% select(Name, `Gene symbol`, `Hierarchy node`, Class, Subclass, `% Coverage of reference sequence`, `% Identity to reference sequence`) %>%
  write_tsv(file="ATB_Ecoli_AFP.tsv.gz")
```

Or you can use the function `atb_amrfp_filter_by_taxa` in the file `AllTheBacteria_functions.R`
