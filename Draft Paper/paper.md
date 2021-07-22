---
output:
  html_document: default
  pdf_document: default
---



---
title: 'RNAheteroplasmy'
aas-journal: Astrophysical Journal <- The name of the AAS journal.
date: "17 June 2021"
bibliography: paper.bib
authors:
- affiliation: 1, 2, 3
  name: Gabriele Lubatti^[first author]
- affiliation: 1,2, 3
  name: Antonio Scialdone^[corresponding author]
affiliations:
- index: 1
  name: Institute of Epigenetics and Stem Cells, Helmholtz Zentrum München, Munich, Germany
- index: 2
  name: Institute of Functional Epigenetics, Helmholtz Zentrum München, Neuherberg, Germany
- index: 3
  name: Institute of Computational Biology, Helmholtz Zentrum München, Neuherberg, Germany
tags:
- R
- bioinformatics
- single cell RNA seq
- heteroplasmy
aas-doi: 10.3847/xxxxx 
---

# Summary
A given base on the DNA can be occupied by 4 different alleles (A, C, T and G). The DNA heteroplasmy measures the frequency of different DNA variants in a given base of the genome. For the mitochondrial DNA (mtDNA), mutation rates are estimated to be 10 to 100 fold higher than for nuclear DNA and therefore mtDNA variants can achieve an higher level of heteroplasmy [@Khrapko13798]. Heteroplasmy from mtDNA data helps to distinguish cells that are in different states or fate commitment as well as for lineage tracing analysis [@LUDWIG20191325], [@Miller2021]. It has been recently shown that single cell RNA seq (scRNA-seq) can be used to reliably identify mtDNA variants, although with a lower statistical power compared to more direct approaches, like mtDNA sequencing [@LUDWIG20191325]. Here we provide RNAheteroplasmy, an highly user friendly R library that allow the quantification of heteroplasmy from RNA seq data starting directly from the raw fastq files.

# Statement of need

Although some tools (based on C or python  [@Huang2021],[@Prashant2020]) have been developed recently to asses heteroplasmy from scRNA-seq data, a unified approach directly from raw fastq files all implemented in R that include both the quantification of heteroplamsy with some more down stream analysis (i.e statistical tests, plotting and clustering based on heteroplasmy values) is not yet available.
`RNAheteroplasmy` is an R package developed in order to fill this gap by providing a user friendly set of functions for the detection of allele frequencies, heteroplasmy and for additional down-stream analysis: 
The two main functions of `RNAheteroplasmy` are:

1. `get_raw_counts_allele`: an highly parallelized function that rely on Rsamtools, used to generate the raw counts matrix starting from fastq files, with cells as rows and bases with the four possible allele as columns.
2. `get_heteroplasmy`: Starting from the ouput of `get_raw_counts_allele`, it computes the matrix with heteroplasmy values (defined as 1 minus the frequency of the most common allele) and the matrix with allele frequencies values,for all the cells and bases that pass a filtering step procedure.

Among the down stream analysis implemented in the package there are:
1. Several statistical tests for identification of most different bases according to heteroplasmy, between group of cells or along a trajectory of cells (i.e cells sorted according to a diffusion pseudo time).
2. Plotting functions for the visualization of heteroplamsy and corresponding allele frequencies values among cells
3. Hierarchical clustering on cells based on a distance matrix defined from angular distance of allele frequencies, that could be relevant for lineage tracing analysis 

The main concepts implemented in the current package were proven to be extremely powerful in our recent paper [@Lima2020], where we revealed that cells eliminated by apoptosis in the mouse embryo at the beginnig of gastrulation showed an higher level of mt-DNA heteroplasmy.



# References


