# **CHROMACS_V1**
> **Graphical user interface for chromatin accessibility studies from ATAC-seq datasets**



## 📖 About

CHROMACS (Chromatin Accessibility Studies) provides a streamlined, GUI-based workflow for analyzing ATAC-seq data, including:

- Quality Control (`FastQC, MultiQC`) ==> of both raw reads and trimmed reads;

- Trimming (`Trim Galore`) ==> optional;

- Genome Alignment and Coverage Analysis(`Bowtie2, Samtools, DeepTools`) ==> flexible with several model organisms and beyond, (uses Ensembl toplevel reference);

- Peak Calling (`Genrich or MACS3`) ==> employes two different peak calling tools, users can select any of the two, with several options for parameter setup;

- Peak Annotation (`ChIPseeker`) ==> builds own txDB by fetching GTF files from ENSEMBL, uses default orgdb from bioconductor;

- Differential Peak Analysis (`Diffbind`) ==> employs diffbind to visualize differential studies and also generates essential up- and downgraded regions; includes annotation of those regions as well


Designed to assist both beginners and experienced users in analyzing chromatin accessibility with minimal command-line usage



## 📥 Installation

#### Step 1: Clone the repository (Recommended) :
```bash
git clone https://github.com/MPrince2703/CHROMACS_V1.git
cd CHROMACS_V1
```

#### Step 1 (alternate): Download as ZIP
If you download the zip file from the repository, ensure to avoid the above step

1. Go to https://github.com/MPrince2703/CHROMACS_V1
2. Click the green Code button → Download ZIP
3. Extract the ZIP file.
The folder will be named CHROMACS_V1-main. Configure to the CHROMACS_V1-main
```bash
cd CHROMACS_V1-main
```

#### Step 2: Create and activate the environment :
```bash
conda env create -f environment.yml
conda activate chromacs
```

#### Step 3: Install the package :
```bash
pip install .
```

#### Step 4: Launch the application :
```bash
chromacs
```


## Acknowledgements

- FastQC, MultiQC, Trim Galore, Bowtie2, MACS3, Genrich, ChIPseeker, DiffBind; along with their dependencies
- Bioconda, Conda-Forge community
