# **Gobi2025miRNA**
## *2-weeks project for a group of Bioinformatics students from the GoBi course (Feb 2025) on the topic of miRNA*

### **From chimeric reads to miRNA/isomiR to target interactions**

This pipeline is based on the following two GitHub repositories:
- [SCRAP](https://github.com/Meffert-Lab/SCRAP/blob/main/README.md)
- [miRNA Project GoBi 2024](https://github.com/giulic3/mirna-project-gobi-2024)

## **Content**
```
├── adapters/                    # Adapter sequences for preprocessing
├── annotations/                 # miRBase annotation files
├── bin/                         # Source code and scripts
├── README.md                    # Project documentation
└── SCRAP_environment.yml        # Conda environment file for dependency setup
```
## bin/ Directory Contents
The `bin/` directory contains various scripts used in the analysis pipeline. Below is a brief description of each script:

- **Peak_Calling_Histo.sh** - Shell script for generating histograms related to peak calling.
- **Peak_Annotation.sh** - Annotates detected peaks with reference genome information.
- **Peak_Calling.sh** - Main script for performing peak calling from sequencing data.
- **Reference_Installation.sh** - Installs and configures reference genome files for analysis.
- **SCRAP.sh** - Runs the SCRAP pipeline for miRNA/isomiR detection.
- **SCRAP_original.sh** - Original, unmodified version of the SCRAP pipeline.

### **bin/downstream_analysis/ Directory Contents**
This subdirectory contains scripts for downstream data analysis after peak calling.

- **2a_miRNA_distribution.py** - Python script for analyzing miRNA distribution.
- **2b_Peak_Calling_Histo.sh** - Shell script for generating histograms related to miRNA peak calling.
- **2c_Find_Shared_3UTR_Targets.sh** - Identifies miRNA target genes with shared 3' UTR regions.
- **2c_find_shared_targets.sh** - Detects shared miRNA target genes.
- **2c_plot_histogram.py** - Python script for visualizing distributions of shared miRNA targets.

## **Setup**
```bash
conda env create -f SCRAP_environment.yml
conda activate scrap_env  # Adjust the name if different
```

## **Usage**
```bash
# Clone the repository
git clone https://github.com/your-repo/Gobi2025miRNA.git
cd Gobi2025miRNA

# Activate the environment
conda activate scrap_env

# Run the pipeline (modify command as needed)
bash run_pipeline.sh
```