# MitoHEAR
MitoHEAR (**Mito**chondrial **HE**teroplasmy **A**nalyze**R**) is an R package that allows the estimation as well as downstream statistical analysis of the mtDNA heteroplasmy calculated from single-cell datasets.

## Installation
Before installing MitoHEAR, the following packages should be installed first:
```

if (!requireNamespace("BiocManager", quietly = TRUE)){
install.packages("BiocManager")}
if (!require("Rsamtools")) {BiocManager::install("Rsamtools")}
if (!require("Biostrings")) {BiocManager::install("Biostrings")}
if (!require("GenomicRanges")) {BiocManager::install("GenomicRanges")}
if (!require("regioneR")) {BiocManager::install("regioneR")}
if (!require("IRanges")) {BiocManager::install("IRanges")}
if (!require("karyoploteR")) {BiocManager::install("karyoploteR")}
if (!require("ComplexHeatmap")) {BiocManager::install("ComplexHeatmap")}
if (!require("ggplot2")) { install.packages("ggplot2", dependencies = TRUE, repos="http://cran.r-project.org")}
if (!require("data.table")) {install.packages("data.table", dependencies = TRUE, repos="http://cran.r-project.org")}
if (!require("gam")) {install.packages("gam", dependencies = TRUE, repos="http://cran.r-project.org")}
if (!require("rdist")) {install.packages("rdist", dependencies = TRUE, repos="http://cran.r-project.org")}
if (!require("dynamicTreeCut")) {install.packages("dynamicTreeCut", dependencies = TRUE, repos="http://cran.r-project.org")}
if (!require("circlize")) {install.packages("circlize", dependencies = TRUE, repos="http://cran.r-project.org")}
if (!require("rlist")) {install.packages("rlist", dependencies = TRUE, repos="http://cran.r-project.org")}
if (!require("gridExtra")) {install.packages("gridExtra", dependencies = TRUE, repos="http://cran.r-project.org")}
if (!require("clustree")) {install.packages("clustree", dependencies = TRUE, repos="http://cran.r-project.org")}
if (!require("mcclust")) {install.packages("mcclust", dependencies = TRUE, repos="http://cran.r-project.org")}
if (!require("fmsb")) {install.packages("fmsb", dependencies = TRUE, repos="http://cran.r-project.org")}
if (!require("plotly")) {install.packages("plotly", dependencies = TRUE, repos="http://cran.r-project.org")}

```

To install MitoHEAR, please run the following:
```
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("https://github.com/ScialdoneLab/MitoHEAR/tree/master",auth_token="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")
library(MitoHEAR)
```
The package has two main functions: **get_raw_counts_allele** and **get_heteroplasmy**.

## Troubleshooting installation errors
The package was developed within R.4.0.2. A different version of R (higher or lower) may lead to installation failed. If is not possible to install the package following the command in the section **Installation**, then download this GitHub repository (main and master) at the location *your_path* and then run the following:

```
#load required packages and functions
setwd(paste0(your_path,"/MitoHEAR-master/vignettes"))
source("MitoHEAR_packages.R")
source("get_heteroplasmy.R")
source("get_raw_counts_allele.R")
source("lineage_tracing.R")
source("other_functions.R")
source("plotting_functions.R")
source("statistical_test.R")
```
With the command above, all the functions of the library MitoHEAR will be available to use.
The input data in order to reproduce the vignettes in the section **Tutorials** are available at the paths:
1. **/your_path/MitoHEAR-master/inst/extdata**

If MitoHEAR is not installed, it is possible to fully reproduce all the vignettes in the section **Tutorials** using the scripts available in the folder 
**[Tutorials_no_installation](https://github.com/ScialdoneLab/MitoHEAR/tree/master/Tutorials_no_installation):**

## get_raw_counts_allele

```
get_raw_counts_allele(bam_input,path_fasta,cell_names,cores_number=1) 
```

1. **bam_input**: character vector of sorted bam files (one for each sample) with full path.
2. **path_fasta**: fasta file of the genomic region of interest with full path.
3. **cell_names**: character vector with the names of the samples.
4. **cores_number**: Number of cores to use.

In the same location of the sorted bam file, also the corresponding index bam file (.bai) should be present

An example of input could be:
```
load(system.file("extdata", "after_qc.Rda", package = "MitoHEAR"))
cell_names=as.vector(after_qc$new_name)
cell_names[1:5]
[1] "24538_8_14" "24538_8_23" "24538_8_39" "24538_8_40" "24538_8_47"
path_to_bam="/home/ies/gabriele.lubatti/revision_heteroplasmy/Cell_Competition_data/all_unique_bam_files/"
bam_input=paste(path_to_bam,cell_names,".unique.bam",sep="")
path_fasta="/home/ies/gabriele.lubatti/revision_heteroplasmy/heteroplasmy_mt/Genome/Mus_musculus.GRCm38.dna.chromosome.MT.fa"
output_SNP_mt=get_raw_counts_allele(bam_input,path_fasta,cell_names)
```
where **after_qc** is a dataframe with number of rows equal to the number of samples and with columns related to meta data information (i.e. cluster and batch).

The output of **get_raw_counts_allele** is a list with three elements:
1. **matrix_allele_counts**: matrix with rows equal to **cell_names** and with columns equal to the bases in the **path_fasta** file with the four possible alleles. For each pair sample-base there is the information about the counts on the alleles A,C,G and T
2. **name_position_allele**: character vectors with length equal to *n_cols* of **matrix_allele_counts** with information about the name of the bases and the corresponding allele
3. **name_position**: character vectors with information about the name of the bases.

```
load(system.file("extdata", "output_SNP_mt.Rda", package = "MitoHEAR"))
matrix_allele_counts=output_SNP_mt[[1]]
## In this example we have 723 cells and 65196 columns (4 possible alleles for the 16299 bases in the mouse MT genome)
dim(matrix_allele_counts)
[1]   723 65196
head(matrix_allele_counts[1:5,1:5])
           1_A_G_MT 1_C_G_MT 1_G_G_MT 1_T_G_MT 2_A_T_MT
24538_8_14        0        0        0        0        0
24538_8_23        0        0        0        0        0
24538_8_39        0        0        0        0        0
24538_8_40        0        0        0        0        0
24538_8_47        0        0        0        0        0

name_position_allele=output_SNP_mt[[2]]
name_position_allele[1:8]
[1] "1_A_G_MT" "1_C_G_MT" "1_G_G_MT" "1_T_G_MT" "2_A_T_MT" "2_C_T_MT" "2_G_T_MT" "2_T_T_MT"

name_position=output_SNP_mt[[3]]
name_position[1:8]
[1] "1_MT" "1_MT" "1_MT" "1_MT" "2_MT" "2_MT" "2_MT" "2_MT"
```
It is also possible to run a [command line implementation](https://github.com/ScialdoneLab/MitoHEAR/blob/master/vignettes/get_raw_counts_allele_script.R) of the function **get_raw_counts_allele**: For the command line implementation is not necessary to install MitoHEAR, but the following libraries should be installed:
1. **rlist**
2. **Rsamtools**
3. **Biostrings**
4. **data.table**
5. **GenomicRanges**
6. **IRanges**


```
Rscript --vanilla get_raw_counts_allele_script.R -b bam_input -f "Mus_musculus.GRCm38.dna.chromosome.MT.fa" -c cell_names -o output_SNP_mt.Rda -s 20

Options:
	-b CHARACTER, --bam_input=CHARACTER
		character vector of sorted bam files

	-f CHARACTER, --path_fasta=CHARACTER
		fasta file of the genomic region of interest

	-c CHARACTER, --cell_names=CHARACTER
		character vector with the names of the samples that will be used in the down-stream analysis", metavar="character

        -o CHARACTER, --output=CHARACTER
		output file name [default= out.Rda]

	-s INTEGER, --cores_number=INTEGER
		Number of cores to use[default= 1]
```

## get_heteroplasmy
```
get_heteroplasmy(matrix_allele_counts,name_position_allele,name_position,number_reads,number_positions,filtering=1,my.clusters=NULL) 
``` 
starts from the output of **get_raw_counts** and performs a two step filtering procedure, the first on the cells and the second on the bases. The aim is to keep only the cells that have more than **number_reads** counts in more than **number_positions** bases and to keep only the bases that are covered by more than **number_reads** counts in all the cells (**filtering**=1)  or in at least 50% of cells in each cluster (**filtering**=2, with cluster specified by **my.clusters**).

An example of input could be:

```
load(system.file("extdata", "output_SNP_mt.Rda", package = "MitoHEAR"))
load(system.file("extdata", "after_qc.Rda", package = "MitoHEAR"))

## We compute heteroplasmy only for cells that are in the condition "Cell competition OFF" and belong to cluster 1, 3 or 4
row.names(after_qc)=after_qc$new_name
cells_fmk_epi=after_qc[(after_qc$condition=="Cell competition OFF")&(after_qc$cluster==1|after_qc$cluster==3|after_qc$cluster==4),"new_name"]
after_qc_fmk_epi=after_qc[cells_fmk_epi,]
my.clusters=after_qc_fmk_epi$cluster

matrix_allele_counts=output_SNP_mt[[1]]
name_position_allele=output_SNP_mt[[2]]
name_position=output_SNP_mt[[3]]

epiblast_ci=get_heteroplasmy(matrix_allele_counts[cells_fmk_epi,],name_position_allele,name_position,number_reads=50,number_positions=2000,filtering = 2,my.clusters)
```
The output of **get_heteroplasmy** is a list with five elements.
The most relevant elements are the matrix with heteroplasmy values (**heteroplasmy_matrix**) and the matrix with allele frequencies (**allele_matrix**), for all the cells and bases that pass the two step filtering procedures. 
The heteroplasmy is computed as **1-max(f)**, where **f** are the frequencies of the four alleles for every sample-base pair.
For more info about the ouput see **?get_geteroplasmy**.

```
heteroplasmy_matrix=epiblast_cell_competition[[3]]
## 261 cells and 5736 bases pass the two step filtering procedures
dim(heteroplasmy_matrix_ci)
[1]  261 5376
head(heteroplasmy_matrix_ci[1:5,1:5])
           136_MT    138_MT 140_MT      141_MT 142_MT
24538_8_14      0 0.0000000      0 0.000000000      0
24538_8_39      0 0.0000000      0 0.000000000      0
24538_8_40      0 0.0000000      0 0.008695652      0
24538_8_47      0 0.0952381      0 0.000000000      0
24538_8_69      0 0.0000000      0 0.000000000      0

allele_matrix=epiblast_cell_competition[[4]]
## 261 cells and 21504 allele-base (4 possible alleles for the 5736 bases).
dim(allele_matrix_ci)
[1]   261 21504
head(allele_matrix_ci[1:4,1:4])
           136_A_G_MT 136_C_G_MT 136_G_G_MT 136_T_G_MT
24538_8_14          0          0          1          0
24538_8_39          0          0          1          0
24538_8_40          0          0          1          0
24538_8_47          0          0          1          0


```
## Down-stream analysis
**MitoHEAR** offers several ways to extrapolate relevant information from heteroplasmy measurement. 
For the identification of most different bases according to heteroplasmy between two group of cells (i.e. two clusters), an unpaired two-samples Wilcoxon test is performed with the function **get_wilcox_test**.  The heteroplasmy and the corresponding allele frequencies for a specific base can be plotted with **plot_heteroplasmy** and **plot_allele_frequency**. 
If for each sample a diffusion pseudo time information is available, then it is possible to detect the bases whose heteroplasmy changes in a significant way along pseudo-time with **dpt_test** and to plot the trend with **plot_dpt**.
It is also possible to perform a cluster analysis on the samples based on distance matrix obtained from allele frequencies with **clustering_dist_ang** and to visualize an heatmap of the distance matrix with samples sorted according to the cluster result with **heatmap_plot**. This approach could be usufel for lineage tracing analysis.
For more exaustive information about the functions offered by **MitoHEAR** see **Tutorials** section below and the help page of the single functions. (**?function_name**).

## Tutorials

The following tutorials are completely reproducible within the package **MitoHEAR**:

### **[cell_competition_mt_example_notebook.Rmd](https://github.com/ScialdoneLab/MitoHEAR/blob/master/vignettes/cell_competition_mt_example_notebook.Rmd):**
This tutorial uses single cell RNA seq mouse embryo data ([Lima *et al.*, 2021](https://www.nature.com/articles/s42255-021-00422-7?proof=t))(Smart-Seq2 protocol).
The heteroplasmy is computed for the mouse mitochondrial genome.
Identification and plotting of most different bases according to heteroplasmy between clusters (with **get_wilcox_test**, **plot_heteroplasmy** and **plot_allele_frequency**) and along pseudo time (with **dpt_test** and **plot_dpt**) are shown.
The top 10 bases with highest variation in heteroplasmy belong to the genes mt-Rnr1 and mt-Rnr2 and in these positions the heteroplasmy always increases with the diffusion pseudo time. 


### **[cell_competition_ercc_example_notebook.Rmd](https://github.com/ScialdoneLab/MitoHEAR/blob/master/vignettes/cell_competition_ercc_example_notebook.Rmd):**
This tutorial uses single cell RNA seq mouse embryo data ([Lima *et al.*, 2021](https://www.nature.com/articles/s42255-021-00422-7?proof=t))(Smart-Seq2 protocol).
The heteroplasmy is computed for ERCC-Spike-In bases, as a technical control. We should not see a significant variation in heteroplasmy between clusters or along pseudo time, since the only source of variation for ERCC is technical and not due to biological reasons.
Identification and plotting of most different bases according to heteroplasmy between clusters (with **get_wilcox_test**, **plot_heteroplasmy** and **plot_allele_frequency**) and along pseudo time (with **dpt_test** and **plot_dpt**) are shown.

### **[cell_competition_bulk_data_mt__example_notebook.Rmd](https://github.com/ScialdoneLab/MitoHEAR/blob/master/vignettes/cell_competition_bulk_data_mt_example_notebook.Rmd):**
This tutorial uses bulk RNA seq data from data from two mtDNA cell lines( [Lima *et al.*, 2021 ](https://www.nature.com/articles/s42255-021-00422-7?proof=t)). 
Since the mt DNA sequence of the two cell lines BG(Loser-95%) and HB(Winner-24%) is available (https://www.ncbi.nlm.nih.gov/nuccore/KC663619.1  and https://www.ncbi.nlm.nih.gov/nuccore/KC663620.1), we identify the bases that are different from the reference genome and then check that in these positions the level of heteroplasmy given by **MitoHEAR** is very close to the one expected.
Cluster analysis among samples based on allele frequency values (done with **clustering_dist_ang**) reveals that we can perfectly distinguish between the two cell linegaes only by looking at the heteroplasmy values of the mitochondrial bases.

### **[lineage_tracing_example_notebook.Rmd](https://github.com/ScialdoneLab/MitoHEAR/blob/master/vignettes/lineage_tracing_example_notebook.Rmd):**
This tutorial uses single cell RNA seq mouse embryo data from  [Goolam *et al.*, Cell, 2016 ](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-3321/?query=antonio+scialdone)(Smart-Seq2 protocol). 
There are embryos at different stages from 2-cells to 8-cells stage. At each stage, for every cell it is known the embryo of origin.
We illustrate how unsupervised cluster analysis based on allele frequencies information (performed with **clustering_dist_ang**) could be used in order to perform a lineage tracing analysis, by grouping together cells which are from the same embryo.

### **[Ludwig_et_al_example_notebook.Rmd](https://github.com/ScialdoneLab/MitoHEAR/blob/master/vignettes/Ludwig_et_al_example_notebook.Rmd):**
This tutorial uses two single cell RNA seq human cells dataset from  [Ludwig *et al.*, Cell, 2019 ](https://doi.org/10.1016/j.cell.2019.01.022). 
We illustrate how unsupervised cluster analysis based on allele frequencies information (performed with **clustering_dist_ang**) can be used in order to aggregate cells. The result from unsupervised cluster analysis are consistent with previously available information (colonies of cells, donors).



