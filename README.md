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
├── cortex_samples/              # Mouse cortex sample data
├── keratinocyte_samples/        # Mouse keratinocyte sample data
├── README.md                    # Project documentation
└── SCRAP_environment.yml        # Conda environment file for dependency setup
```

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