# Differential Expression Analysis and Gene Set Enrichment Analysis Pipeline

A Differential Expression and Gene set Enrichment analysis workflow for paired-end and single-end sequence data, built using Python, R, Snakemake, and Conda, with a Click-based CLI.

This workflow is designed such that it's reproducible and scalable to multiple samples, conditions, and configurable to different references, gene sets, and organisms.

## Overview

This project takes single-ended and paired-end RNA-seq reads from a sample of your choice through QC, genome alignment, BAM processing, and gene level quantification, DESeq2 analysis, and pathway enrichment (GSEA). Applied here to compare human bronchial epithelial gene expression in cystic fibrosis (CF) vs. non-CF donors.

### Key results

#### Differential expression (DESeq2, CF vs. control):
<img width="4200" height="2400" alt="volcano" src="https://github.com/user-attachments/assets/a2a73083-3b91-4850-a2eb-905a0ad002d6" />

#### Pathway enrichment (GSEA, MSigDB Hallmark gene sets):
<img width="3000" height="2400" alt="barplot" src="https://github.com/user-attachments/assets/12cfc300-2a07-4aa4-8c7d-4c96f62948b2" />

The top enriched-in-CF pathways show a coherent inflammatory/immune signature; TNFa signaling via NFκB (Bodas, M., & Vij, N., _Discovery medicine_, 2010), inflammatory response, interferon-gamma response, complement and allograft rejection along with downregulated oxidative phosphorylation are all consistent with CF airway epithelium's known chronic inflammatory phenotype. 
A few tissue-mismatched hits (e.g. spermatogenesis) also appear in the top-20; rather than being genuine biological findings, they are flagged. Gene sets unrelated to airway epithelium showing up in bulk RNA-seq GSEA is a known artifact of the hallmark collection, and should not be treated as real evidence.

## Workflow
Snakemake-based
Input : Fastq (from CLI fetch/fetchall/local) -> samples.tsv
-> ReadQC and Trimming (Fastp)
-> Alignment (HISAT2 and SAMtools)
-> Quantification (featureCounts)
-> Differnential Expression (DESeq2) -> Volcano plot
-> Enrichment (clusterProfiler, msigdbr) -> Bar plot

The pipeline automatically branches when dealing with single-end vs. paired-end reads. Therefore mixed-layout datasets can be run without manual intervention.

Also, the reference genome, strandedness (for featureCounts), species, gene collection, graph titles, and GSEA ranking metric can all be directly configured through the provided config.yaml in config/.

## Setup
```bash
# Python/alignment environment
conda env create -f envs/RNASeqPipelineProject.yml
conda activate RNASeqPipelineProject

# R/DESeq2 environment (built separately, activated per-rule by Snakemake)
conda env create -f envs/deseq2.yml
```

## Usage
### Add samples:
```bash
# Single accession
python add_sample.py fetch SRR31795696 --name cf_donor1

# Multiple accessions at once
python add_sample.py fetchall SRR31795672 SRR31795673 --name cf_donor1 --name cf_donor2

# Already-downloaded local FASTQs
python add_sample.py local data/raw/

# Assign conditions (positional: first N rows = control, next groups = named treatments)
python add_sample.py condition --control 6 --treatment 6 --names cf
```
### The pipeline run with the following command:
```bash
snakemake --cores 4 --use-conda
```
The HISAT2 index is not provided nor is the annotation file.

## Commands

Samples are registered in `config/samples.tsv` using `add_sample.py`, a Click-based CLI with 4 subcommands.

### fetch

Downloads a single accession via `prefetch` + `fasterq-dump` and registers it as a sample.

```bash
python add_sample.py fetch SRR13970441
python add_sample.py fetch SRR13970441 --name ecoli_ctrl_1
python add_sample.py fetch SRR13970441 --force
```

- `--name` — sample name to use in `samples.tsv` (defaults to the accession)
- `--force` — overwrite an existing sample with the same name

### fetchall

Downloads two or more accessions in one call. Failures (a failed download, missing paired output, or a duplicate sample name) are skipped rather than stopping the batch, and a summary is printed at the end.

```bash
python add_sample.py fetchall SRR000001 SRR000002 SRR000003
python add_sample.py fetchall SRR000001 SRR000002 --name ctrl_1 --name treat_1
```

- `--name` — one name per accession, given in the same order (defaults to each accession)
- `--force` — overwrite existing samples with matching names

For a single accession, use `fetch` instead.

### local

Registers FASTQ files that are already downloaded, either as individual files or a directory. R1/R2 mates are paired automatically by filename (`_1`/`_2` or `_R1`/`_R2`).

```bash
python add_sample.py local data/raw/mysample_R1.fastq.gz data/raw/mysample_R2.fastq.gz
python add_sample.py local data/raw/
```

- `--force` — overwrite existing samples with matching names

### condition

Adds a separate column to the samples.tsv file for control(s) and treatment(s)/condition(s), this column allows the DESeq2 rule to automatically match and perform analysis.

```bash
python add_sample.py condition --control 3 --treatment 3 --treatment 3 --names penicillin --names ampicillin
```
The above command will add control to the condition column for the first 3 rows, and penicillin to rows 4-6, and ampicillin to rows 7-9.

```bash
python add_sample.py condition --control 2 --treatment 2 --names treatment1
```
The above command will add control to the condition column for the first 2 rows, and will add treatment1 to rows 3-4.

## Tools

Snakemake | Workflow management
Fastp | Read Quality Control and Trimming
HISAT2 | Read alignment
SAMtools | BAM sorting and indexing
featureCounts | Gene-level quantification
MultiQC | QC report aggregation
DESeq2 | Differential Expression Analysis
Ncbi-tools-cli | datasets, testing and validation

## Data
Organism:
Tested with Homo sapiens and E.coli.

A note on dataset selection: an earlier pull (PRJNA1478899) initially looked usable but turned out, on closer inspection of BioSample accessions, to be a single CF patient and a single non-CF patient each sequenced twice, technical replicates, not biological ones. This led to impossibly high -log10p values (>400) and led to me switching out the dataset for current one which produced much more reasonable numbers. The other dataset also have several genes with reversed log2 fold change values (upregulated -> downregulated, vice versa) which drew critical attention and led to ultimately pulling the plug on that dataset. 

Reference genome: Homo sapiens GRCh38 v116, aligned with HISAT2's prebuilt grch38_tran index. One issue may be that HISAT2's index represents an earlier version of the reference (v84), but no apparent issues arose during test runs.

Library strandedness: The featureCounts was done using -s 2, for TruSeq Stranded mRNA since it's reverse stranded. Though this is configurable via the config file.

Example run done with GSE285099. (Szachowicz, P. J. et al, _iScience_, 2025)
6 CF donors, vs 6 non-CF donors, human bronchial epithelium cells, paired-end reads.
Exact runs:

Control:
SRR31795696
SRR31795697
SRR31795698
SRR31795699
SRR31795700
SRR31795701

Cystic Fibrosis:
SRR31795672
SRR31795673
SRR31795674
SRR31795675
SRR31795676
SRR31795677

Reads:
Single-end FASTQ,
Paired-end FASTQ

## Environment

This pipeline uses Conda environments defined in:

`envs/RNASeqPipelineProject.yml`
`envs/deseq2.yml`

Snakemake directly manages the Conda environments, so no need to activate the environment as long as you have Snakemake installed.

## Output

The pipeline generated:

- Fastp quality-control reports
- MultiQC summary report
- Sorted and indexed BAM
- SAMtools alignment QC
- featureCounts gene-level count table
- Count Matrix for input into DESeq2
- DESeq2 Analysis
- Volcano plot of DESeq2 Analysis
- Enrichment Analysis
- Bar plot of Enrichment Analysis

## Limitations
- Hallmark GSEA output should be assessed critically, not taken at face value. A couple of top-20 pathways have no plausible connection to bronchial epithelium (e.g. spermatogenesis). Gene sets not credibly linked with associated samples is a known pattern and limitation of the Hallmark collection with RNA-seq.

- Example uses N = 6, which is a modest sample size, and its results should be taken as exploratory rather as clinically definitive.

- Gene ID mapping (Ensembl -> gene symbol) was done with multivals = "first", with rare one to many mappings done to the first hit. Certain gene symbols may be arbitrarily chosen among certain valid alternatives.


## Project Structure

```text
RNASeqPipelineProject/
|- Snakefile
|- add_sample.py
|- README.md
|- config/
|	- config.yaml
|	- samples.tsv
|- data/
|	- raw/
|- envs/
|	- RNASeqPipelineProject.yml
|	- deseq2.yml
|- reference/
|	- hisat2_index/
|- results/
|	- bam/
|	- counts/
|	- fastp/
| - qc/
| - matrix/
| - trimmed/
| - matrix
|	- multiqc/
| - deseq/
| - volcano/
| - gsea/ 
|-scripts/
|	- DE.R
|	- make_matrix.py
| - volcano.R
| - gsea.R
```
## Future Work

In the future, I would like implement measures to improve accuracy and switch to (or at least provide the option for different alignment strategies).
### Plans
- Transcript-level Pseudo-Alignment using Salmon
- UMI Dedup
- Automatic Strandedness Inference
- Post-processing using Picard

## Citations

Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome biology, 15(12), 550. https://doi.org/10.1186/s13059-014-0550-8

Chen, J. W., Shrestha, L., Green, G., Leier, A., & Marquez-Lago, T. T. (2023). The hitchhikers' guide to RNA sequencing and functional analysis. Briefings in bioinformatics, 24(1), bbac529. https://doi.org/10.1093/bib/bbac529

Szachowicz, P. J., Wohlford-Lenane, C., Donelson, C. J., Ghimire, S., Thurman, A., Xue, B., Boly, T. J., Verma, A., MašinoviĆ, L., Bermick, J. R., Rehman, T., Perlman, S., Meyerholz, D. K., Pezzulo, A. A., Zhang, Y., Smith, R. J. H., & McCray, P. B., Jr (2025). Complement is primarily activated in the lung in a mouse model of severe COVID-19. iScience, 28(3), 111930. https://doi.org/10.1016/j.isci.2025.111930

Bodas, M., & Vij, N. (2010). The NF-kappaB signaling in cystic fibrosis lung disease: pathophysiology and therapeutic potential. Discovery medicine, 9(47), 346–356.

Strubberg, A. M., Liu, J., Walker, N. M., Stefanski, C. D., MacLeod, R. J., Magness, S. T., & Clarke, L. L. (2017). Cftr Modulates Wnt/β-Catenin Signaling and Stem Cell Proliferation in Murine Intestine. Cellular and molecular gastroenterology and hepatology, 5(3), 253–271. https://doi.org/10.1016/j.jcmgh.2017.11.013

