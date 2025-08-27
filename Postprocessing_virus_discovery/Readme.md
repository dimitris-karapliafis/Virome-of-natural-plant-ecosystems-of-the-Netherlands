# Virus Discovery Post-Processing Pipeline

This directory contains a collection of scripts and workflows designed for the post-processing and analysis of viral contigs identified from metagenomic data. The pipeline takes raw viral operational taxonomic units (vOTUs) and enriches them with taxonomic, host, and abundance information, culminating in a preliminary ecological analysis.

---

### 1. `Contamination_check/`

This module is responsible for evaluating potential contamination from the positive control to the rest of the sequencing libaries.

*   `contamination_control.Rmd`: An R Markdown script that  implements the logic for checking for the presence of Bacopa 
Chlorosis Virus RNA2 in the libraries of this study.

### 2. `Ribodepletion/`

This directory contains scripts to assess the effectiveness of the ribosomal RNA (rRNA) depletion step performed during library preparation.

*   `ribodepletion_markdown.Rmd`: An R Markdown file for analyzing the proportion of rRNA reads remaining in the trimmed data.
*   `ribodepletion_summary.tsv`: A summary table with statistics on rRNA content per sample.

### 3. `Read_coverage_normalisation/`

This module focuses on quantifying the abundance of each vOTU by mapping reads back to the contigs and normalizing the coverage values. This allows for fair comparison of vOTU abundances across different samples.

*   `normalisation_markdown.Rmd`: An R Markdown script that details the normalization procedures.
*   The `.tsv` file `master_clusters_custom_rpkm_bwa2.tsv` contains the normalized abundance matrices using RPKM (Reads Per Kilobase of transcript per Million reads) 

### 4. `CoverM/`

This directory likely contains brief explanations related to the [CoverM](https://github.com/wwood/CoverM) tool, which is used for calculating the coverage and relative abundance of genomes or contigs from metagenomic data. The results from this tool are a key input for the `Read_coverage_normalisation` step.

### 5. `easy-taxonomy/`

This directory is likely used for scripts and results related to the taxonomic classification of vOTUs, probably using the `MMseqs2 easy-taxonomy` workflow. The output from this step is essential for both the host association and ecological analysis.

### 6. `Host_association/`

After taxonomic classification, this module aims to predict the hosts of the identified vOTUs. It integrates information from multiple annotation sources and ICTV to infer virus-host relationships. See the `Host_association/Readme.md` for a detailed explanation of this step.

### 7. `Preliminary_ecological_analysis/`

An initial ecological investigation of the alpha- and beta- diversity of the natural plant ecosystems of the Netherlands

*   `ecological_analysis_markdown.Rmd`: An R Markdown script for performing analyses such as alpha/beta diversity, differential abundance, and exploring the relationships between viral communities and environmental variables.
*   `master_clusters_custom_rpkm_bwa2.tsv`: The abundance data used as input for the ecological analysis.
