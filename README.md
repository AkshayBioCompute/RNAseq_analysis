Here's the updated `README.md` file that incorporates the newly added folder structure and explains the contents and usage of each section:

---

# RNA-Seq Pipeline Workflow

This repository provides a comprehensive RNA-seq pipeline workflow designed for high-throughput RNA sequencing data analysis. It supports multiple workflow orchestration systems, including **Nextflow**, **Snakemake**, and **WDL**. The pipeline handles all steps from quality control and read trimming to alignment and feature counting, enabling efficient and reproducible RNA-seq data analysis.

## Prerequisites

Before running the pipeline, ensure that the following tools are installed:

- **Software**:
  - [HISAT2](https://daehwankimlab.github.io/hisat2/) – for read alignment.
  - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) – for quality control of sequencing data.
  - [MultiQC](https://multiqc.info/) – to aggregate quality control results.
  - [Fastp](https://github.com/OpenGene/fastp) – for read trimming.
  - [SAMtools](http://www.htslib.org/) – for manipulating SAM/BAM files.
  - [featureCounts](http://subread.sourceforge.net/) – for counting reads per gene.

- **Python packages** (listed in `requirements.txt`):
  - `pandas==1.5.0`
  - `numpy==1.23.3`
  - `matplotlib==3.5.3`
  - `biopython==1.81`

To install the required Python dependencies, run:

```bash
pip install -r requirements.txt
```

**Perl** is also required for running ANNOVAR, but is not included in this pipeline.

## Folder Structure

This repository is organized into the following main directories:

```bash
├── pipeline.wdl         # WDL workflow for RNA-seq analysis
├── nexflow.nf           # Nextflow script for RNA-seq pipeline
├── snakefile            # Snakemake script for RNA-seq pipeline
├── InputWDL.json        # Input parameters for the WDL workflow
├── requirements.txt     # Python package dependencies
└── src/                 # Contains all shell scripts for each process in the pipeline
    ├── unzip_reference.sh
    ├── build_hisat2_index.sh
    ├── generate_md5.sh
    ├── run_fastqc.sh
    ├── run_multiqc.sh
    ├── trim_reads.sh
    ├── align_reads.sh
    ├── convert_sam_to_bam.sh
    └── count_features.sh
```

### Data Directories

```bash
/Data_reproducibility_pattern_variation
│   ├── graph
│   │   ├── PCA_plot_high_res.png
│   │   ├── correlation_heatmap_high_res.png
│   │   ├── scatterplot_matrix_high_res.png
│   │   └── scree_plot_high_res.png
│
/Differential_expression
│   ├── DEG
│   │   ├── All_Time_Specific_DEGs.csv
│   │   ├── All_Tissue_Specific_DEGs.csv
│   │   ├── Interaction_Heart_ZT12_vs_ZT0.csv
│   │   ├── Interaction_Liver_ZT12_vs_ZT0.csv
│   │   ├── Interaction_Liver_vs_Heart_ZT0.csv
│   │   └── Interaction_Liver_vs_Heart_ZT12.csv
│   ├── Graph
│   │   ├── Heatmaps
│   │   └── Volcano
│   └── TopGenes
│       ├── Top_Interaction_Heart_ZT12_vs_ZT0.csv
│       ├── Top_Interaction_Liver_ZT12_vs_ZT0.csv
│       ├── Top_Interaction_Liver_vs_Heart_ZT0.csv
│       ├── Top_Interaction_Liver_vs_Heart_ZT12.csv
│       ├── Top_Time_Specific_DEGs.csv
│       └── Top_Tissue_Specific_DEGs.csv
│
/Functional_Enrichment_Analysis
│   ├── GO
│   │   └── (Enrichment analysis results)
│   ├── Graph
│   │   └── (Graphical outputs like enrichment plots)
│   └── KEGG
│       └── (KEGG pathway analysis results)
│
/RNAseq
│   ├── Data
│   │   └── (Files such as multiqc reports, alignment statistics)
│   ├── multiqc
│   │   └── (MultiQC aggregated results)
│   ├── alignment
│   │   └── (Alignment files, BAM, SAM)
│   └── result
│       └── (Final RNA-seq analysis results)
│
/script
│   ├── Data_reproducibility_pattern_variation
│   │   └── (R scripts for reproducibility analysis)
│   ├── Differential_Expression_Analysis
│   │   └── (R scripts for differential expression analysis)
│   ├── Functional_Enrichment_analysis
│   │   └── (R scripts for enrichment analysis)
│   ├── RNAseq_pipeline_Python
│   │   └── (Python scripts for RNA-seq pipeline)
│   └── RNAseq_pipeline_Shell
│       └── (Shell scripts for RNA-seq pipeline)
│
/src
│   ├── RNAseq
│   │   └── (Main RNA-seq analysis scripts)
│
/InputWDL.json
│   └── (Input configuration file for WDL workflow)
│
/README.md
│   └── (Documentation and instructions for using the RNA-seq pipeline)
│
/RNAseq_WDL.wdl
│   └── (WDL workflow file for RNA-seq pipeline)
│
/RNAseq_nexflow.nf
│   └── (Nextflow script for RNA-seq pipeline)
│
/requirements.txt
│   └── (Python dependencies required for analysis)
│
/snakefile
│   └── (Snakemake file for RNA-seq pipeline)
```

## Workflow Overview

The RNA-seq pipeline includes the following major steps:

1. **Unzip Reference Genome**:
   - Unzips the reference genome using `gunzip`.

2. **Build HISAT2 Index**:
   - Builds a HISAT2 index using the reference genome.

3. **Generate MD5 Checksums**:
   - Generates MD5 checksums for input FASTQ files to ensure data integrity.

4. **Run FastQC**:
   - Quality control for FASTQ files using FastQC.

5. **Run MultiQC**:
   - Aggregates all FastQC reports into a single summary using MultiQC.

6. **Trim Reads with Fastp**:
   - Trims paired-end reads using Fastp to improve read quality.

7. **Align Reads with HISAT2**:
   - Aligns the trimmed reads to the reference genome using HISAT2.

8. **Convert SAM to Sorted BAM**:
   - Converts SAM files to sorted BAM files and generates alignment statistics.

9. **Count Features with featureCounts**:
   - Counts the number of reads mapped to each gene using featureCounts.

## Usage

### 1. Install Python dependencies

Install the required Python dependencies:

```bash
pip install -r requirements.txt
```

### 2. Running the Workflow

#### For Nextflow:

```bash
nextflow run nexflow.nf
```

#### For Snakemake:

```bash
snakemake --snakefile snakefile
```

#### For WDL:

Ensure Cromwell is installed and use the following command:

```bash
java -jar cromwell.jar run pipeline.wdl -i InputWDL.json
```

### 3. Shell Scripts

The `src/` folder contains shell scripts used in the pipeline. These handle specific tasks such as aligning reads, trimming them, and generating statistics.

### 4. Outputs

The pipeline generates various outputs:

- **MD5 Checksums**: Ensures data integrity.
- **FastQC Reports**: Quality control metrics for the FASTQ files.
- **MultiQC Report**: Aggregated results from FastQC.
- **Trimmed Reads**: FASTQ files with trimmed reads.
- **BAM Files**: Sorted BAM files with aligned reads.
- **Gene Counts**: A file containing counts of reads mapped to genes.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---

This `README.md` file should guide users through the pipeline and its structure. Let me know if you need further adjustments!
