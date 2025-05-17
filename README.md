# InterOmics App: Integration Scripts
Brandao et al.
Last updated on: 17-05-2025

## Description

InterOmics is a multiOmic App to integrate, in biological terms, Omics
data. In this section, We describe a R script to integrate Genomic and
Transcriptome data sourced from Next generation sequencing.

## Install and Requirements

1.  \*\* Input File \*\* You need to have previously processed both WES
    and RNAseq data in csv format. The file must to have at least the
    following columns:

| sample | omics | genotype | Chr | Start | End | Ref | Alt | Gene.refGene | Func.refGene | ExonicFunc.refGene | AF | avsnp154 | QUAL | INFO_DP | ref.count | count |
|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
| Sample ID | exome or rnaseq | het or hom | Chromosome (e.g., chr1) | Start base position | End base position | Reference allele (e.g., A) | Alternative allele (e.g., C) | Gene symbol (e.g., MBL2) | Variant location in gene (e.g., exonic, upstream) | Variant consequence (e.g., nonsynonymous SNV) | GnomAD allele frequency | dbSNP 154 ID (e.g., rs…) | Base quality score (e.g., 50) | Depth (e.g., 142) | GTEx expression (TPM) | Sample expression (TPM) |

1.  \*\* R version 4.4.2+ (<https://www.r-project.org>) \*\*

-   the package [tidyverse v.2.0.0](https://www.tidyverse.org/).

## Usage of InterOmics:

### Read files in R:

We present an object called *multiomics_individual_data* with an ideal
input for the R commands that follows. In order to execute de commands,

<details class="code-fold">
<summary>Code</summary>

``` r
library(tidyverse)
multiomics_individual_data = data.table::fread("multiomics_individual_data.zip")
```

</details>

### Acessing ASE variants

Allele Specific Expression (ASE) variants can be retrieved from the
dataset using the command:

<details class="code-fold">
<summary>Code</summary>

``` r
ase = multiomics_individual_data %>% 
  group_by(sample, avsnp154) %>%
  filter(n_distinct(omics) == 2) %>%
  filter(any(omics == "exome" & genotype == "het") & 
         any(omics == "rnaseq" & genotype == "hom")) %>%
  ungroup() %>% 
  filter(Chr %in% c(paste0("chr", seq(1,22)), "chrX", "chrY")) %>% 
  filter(omics == "rnaseq") %>% 
  filter(!Func.refGene %in% c("intronic", "intergenic", "upstream",
                              "downstream")) %>%
  filter(ExonicFunc.refGene != "synonymous SNV") %>%
  filter(AF < 0.01 | is.na(AF) == TRUE) %>%
  distinct(avsnp154, sample, .keep_all = TRUE)
```

</details>

### Acessing RNAe variants

RNA edited variants, i.e., those that have changed during the mRNA
processing stage, can be found by using the command:

<details class="code-fold">
<summary>Code</summary>

``` r
rna_e = multiomics_individual_data %>%
  group_by(sample, avsnp154) %>%
  filter(n_distinct(omics) == 2) %>%
  filter(any(omics == "exome" & genotype == "hom") & 
         any(omics == "rnaseq" & genotype == "het")) %>%
  ungroup() %>% 
  filter(Chr %in% c(paste0("chr", seq(1,22)), "chrX", "chrY")) %>% 
  filter(omics == "rnaseq") %>% 
  filter(!Func.refGene %in% c("intronic", "intergenic", "upstream",
                              "downstream")) %>%
  filter(ExonicFunc.refGene != "synonymous SNV") %>%
  filter(AF < 0.01 | is.na(AF) == TRUE) %>%
  distinct(avsnp154, sample, .keep_all = TRUE)
```

</details>

### Acessing NMD variants

Nonsense-mediated messenger RNA (mRNA) decay (NMD) can be seen using
this command line:

<details class="code-fold">
<summary>Code</summary>

``` r
nmd = multiomics_individual_data %>% 
  filter(Chr %in% c(paste0("chr", seq(1,22)), "chrX", "chrY")) %>% 
  filter(omics == "exome") %>% 
  filter(grepl("stop|start|^frameshift", ExonicFunc.refGene) == TRUE) %>% 
  filter(AF < 0.01 | is.na(AF) == TRUE) %>%
  mutate_at(.vars = c("ref.count", "counts"),
            .funs = ~ if_else(is.na(.x), 0, .x)) %>% 
  mutate(ref.count_categ = case_when(ref.count <  0.5 ~ "below cutoff",
                                     ref.count > 0.5 & ref.count <= 10   ~ "low",
                                     ref.count > 10  & ref.count <= 1000 ~ "medium",
                                     ref.count > 1000 ~ "high"
                                     )) %>% 
  mutate(counts_categ = case_when(counts <  0.5 ~ "below cutoff",
                                  counts > 0.5 & counts <= 10   ~ "low",
                                  counts > 10  & counts <= 1000 ~ "medium",
                                  counts > 1000 ~ "high"
  )) %>%
  mutate(ref.count_categ = fct_relevel(ref.count_categ, 
                                    c("below cutoff", "low",
                                      "medium", "high"))) %>% 
  mutate(counts_categ = fct_relevel(counts_categ, 
                                    c("below cutoff", "low",
                                      "medium", "high"))) %>% 
  filter(counts_categ %in% c("below cutoff", "low")) %>% 
  filter(ref.count_categ %in% c("medium", "high")) 
```

</details>

### Acessing Gain of Function (GoF) variants

GoF mutations can be fetched using this command:

<details class="code-fold">
<summary>Code</summary>

``` r
gof = multiomics_individual_data %>% 
  filter(Chr %in% c(paste0("chr", seq(1,22)), "chrX", "chrY")) %>% 
  filter(omics == "exome") %>% 
  filter(!Func.refGene %in% c("ncRNA_intronic", "intronic", "intergenic")) %>%
  filter(grepl("synonymous SNV|^start|^stop|^frameshift|unknown", 
               ExonicFunc.refGene) == FALSE) %>%
  filter(AF < 0.01 | is.na(AF) == TRUE) %>%
  mutate_at(.vars = c("ref.count", "counts"),
            .funs = ~ if_else(is.na(.x), 0, .x)) %>% 
  mutate(ref.count_categ = case_when(ref.count <  0.5 ~ "below cutoff",
                                     ref.count > 0.5 & ref.count <= 10   ~ "low",
                                     ref.count > 10  & ref.count <= 1000 ~ "medium",
                                     ref.count > 1000 ~ "high"
  )) %>% 
  mutate(counts_categ = case_when(counts <  0.5 ~ "below cutoff",
                                  counts > 0.5 & counts <= 10   ~ "low",
                                  counts > 10  & counts <= 1000 ~ "medium",
                                  counts > 1000 ~ "high"
  )) %>%
  mutate(ref.count_categ = fct_relevel(ref.count_categ, 
                                       c("below cutoff", "low",
                                         "medium", "high"))) %>% 
  mutate(counts_categ = fct_relevel(counts_categ, 
                                    c("below cutoff", "low",
                                      "medium", "high"))) %>% 
  filter(counts_categ %in% c("medium", "high")) %>% 
  filter(ref.count_categ %in% c("below cutoff", "low"))
```

</details>

### Generate tables

To export the tables, we can iterate a function to create .tsv files
that can be easily read by and spreadsheet software.

<details class="code-fold">
<summary>Code</summary>

``` r
for(i in c("ase", "gof", "nmd", "rna_e")){
  write.table(get(i), 
            file = paste0(i, ".tsv"), 
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE, 
            quote = FALSE)
}
```

</details>
