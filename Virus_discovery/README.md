# Virus Discovery Pipeline

This pipeline is designed to discover and annotate RNA viral contigs from metatranscriptomic data. It takes assembled contigs as input and performs a series of steps including viral domain searching, annotation against various databases, and filtering to produce a high-confidence set of viral operational taxonomic units (vOTUs).

## Dependencies

### Software
- [SeqKit](https://bioinf.shenwei.me/seqkit/)
- [EMBOSS Transeq](https://www.ebi.ac.uk/Tools/st/emboss_transeq/)
- [HMMER](http://hmmer.org/)
- [MMseqs2](https://github.com/soedinglab/mmseqs2)
- [DIAMOND](https://github.com/bbuchfink/diamond)
- [Palmscan](https://github.com/rcedgar/palmscan)
- [GeNomad](https://portal.nersc.gov/genomad/)

### Python Packages
- pandas
- upsetplot
- matplotlib
- seaborn
- rich
- pyhmmer
- biopython

### Databases
The pipeline requires several databases for HMM searches and blast-based annotations. The paths to these databases are hardcoded in the scripts and need to be configured before running.

- RdRp-scan HMM profiles
- RVMT HMM profiles
- RdRp-scan diamond database
- Genomad database
- NCBI NR diamond database
- Plant Virus Database (for Palmscan)
- PalmDB
- Embryophyta database (nucleotide)
- MMseqs NT database

## How the Scripts Work

The pipeline is modular, with different scripts responsible for specific stages of the analysis.

### `virome_analysis.py`
This is the main script that drives the entire pipeline. It calls functions from the other modules in the correct order:
1.  **Virus Discovery**: Identifies potential viral contigs.
2.  **Annotation**: Annotates the candidate viral contigs.
3.  **Filtering**: Filters the annotated contigs to produce the final vOTU set.

### `vir_disc_module.py`
This module focuses on the initial discovery of viral contigs.
- **Length Filtering**: Removes contigs shorter than 400bp using `seqkit`.
- **Translation**: Translates nucleotide contigs into all six reading frames using `transeq`.
- **HMM Search**: Searches for RdRp protein domains (RdRp, RVMT) in the translated sequences using `pyhmmer`.
- **Parsing**: The `hmmsearch_parser` class processes the HMM search output, calculating coverage and normalized bitscores.
- **Extraction**: Contigs with significant HMM hits are extracted for further analysis.
- **Master File Creation**: Combines all high-confidence viral contigs into a single "master_contigs.fasta" file.

### `annotation_module.py`
This module annotates the viral contigs found in the discovery phase.
- **Clustering**: Groups similar contigs using `mmseqs easy-cluster`.
- **Homology Search**: Performs `diamond blastx` searches against multiple protein databases:
    - Plant Virus Database Palmprints
    - PalmDB
    - NCBI NR
- **Palmscan**: Runs `palmscan` to identify RdRp domains.
- **GeNomad**: Runs `genomad` to identify viruses and plasmids.
- **MMseqs Search**: Searches against nucleotide databases (Embryophyta, MMseqs NT) to identify potential host contamination.
- **Aggregation**: Combines all annotation results into a single comprehensive table.

### `filtering_module.py`
This module applies a series of stringent filters to the annotated contigs to remove false positives.
- **Control Removal**: Excludes any control sequences from the analysis.
- **EVE Filtering**: Identifies and removes potential Endogenous Viral Elements (EVEs) by checking for homology to plant genomes and cellular sequences in the NR.
- **Evidence-based Filtering**: Filters contigs based on multiple lines of evidence:
    - Presence of a Palmscan hit.
    - High HMM profile coverage.
    - Significant hits to the PalmDB.
- **vOTU generation**: Produces the final set of vOTUs in FASTA format and generates summary tables.
- **Visualization**: Can generate Upset plots to visualize the overlap between different annotation methods.

### `logger.py`
A utility script that provides a logging class for writing progress and error messages to both the console and a log file (`vir_disc_pipeline.log`).

## Usage

The pipeline is run from the command line using `virome_analysis.py`.

```bash
python virome_analysis.py -i /path/to/input_directory
```

The input directory should contain subdirectories for each sample, and each sample subdirectory must contain a `contigs.fasta` file.

## Outputs

The main results are organized in a `results/` directory created in the parent directory of your input. Inside, an `overview_results/` directory will contain the most important outputs:

-   `overview_results/viral_discovery/master_contigs.fasta`: FASTA file of all putative viral contigs identified by HMM searches.
-   `overview_results/viral_discovery/summary_filtered_1e-05.tsv`: Table of HMM search results.
-   `overview_results/annotate/master_contigs.annotated.tsv`: A comprehensive annotation table for all putative viral contigs.
-   `overview_results/filtering/vOTUs.fasta`: The final, filtered set of high-confidence viral contigs (vOTUs).
-   `overview_results/filtering/vOTUs.total_contigs.fasta`: All contigs belonging to the vOTU clusters.
-   `overview_results/filtering/eve_candidates.tsv`: Contigs that were filtered out as potential EVEs.
-   `overview_results/filtering/master_contigsno_eves.palmscan.prof_cov.palmdb.tsv`: Contigs that passed all the filtering steps.
