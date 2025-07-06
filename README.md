# **ChromAcS**
> **Graphical user interface for Chromatin Accessibility Studies from ATAC-seq datasets**



## 📖 About

ChromAcS (Chromatin Accessibility Studies) provides a streamlined, GUI-based workflow for analyzing ATAC-seq data, including:

- Quality Control (`FastQC, MultiQC`) ==> of both raw reads and trimmed reads;

- Trimming (`Trim Galore`) ==> optional (users can opt for working with raw data or their own trimmed data saved in "raw data directory");

- Genome Alignment and Coverage Analysis(`Bowtie2, Samtools, DeepTools`) ==> flexible with several model organisms and beyond, (uses Ensembl toplevel reference);

- Peak Calling (`Genrich or MACS3`) ==> employes two different peak calling tools, users can select any of the two, with several options for parameter setup;

- Peak Annotation (`ChIPseeker`) ==> builds own txDB by fetching GTF files from ENSEMBL, uses default orgdb from bioconductor;

- Differential Peak Analysis (`Diffbind or NOISeq`) ==> employs Diffbind to visualize differential studies and also generates essential up- and downgraded regions; includes annotation of those regions as well; in absence of biological replicates, users can opt for NOISeq which allows the same differential studies (however, if you have more than one biological replicates, it is recommended to use Diffbind)


Designed to assist both beginners and experienced users in analyzing chromatin accessibility with minimal command-line usage



## 📥 Installation

#### Step 1: Clone the repository (Recommended) :
```bash
git clone https://github.com/epigen-bioinfolab/CHROMACS.git
cd CHROMACS
```

#### Step 1 (alternate): Download as ZIP
If you download the zip file from the repository, ensure to avoid the above step

1. Go to https://github.com/epigen-bioinfolab/CHROMACS
2. Click the green Code button → Download ZIP
3. Extract the ZIP file.
The folder will be named CHROMACS-main. Configure to the CHROMACS-main
```bash
cd CHROMACS-main
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

- FastQC, MultiQC, Trim Galore, Bowtie2, MACS3, Genrich, ChIPseeker, DiffBind, NOISeq; along with their dependencies
- Bioconda, Conda-Forge community


## 📝 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

