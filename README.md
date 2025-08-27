
<div style="display: flex; align-items: center;">
  <img src="chromacs/assets/ChromAcS.png" alt="ChromAcS Logo" width="120" style="margin-left: 40px;" />
  <div>
  
# **ChromAcS**

> **An Automated, Flexible GUI for End-to-End Reproducible ATAC-seq Analysis Across Multiple Species**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.9%2B-blue)
![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20MacOS-lightgrey)
![Conda](https://img.shields.io/badge/Conda-Required-yellowgreen)
[![GitHub release](https://img.shields.io/github/v/release/epigen-bioinfolab/CHROMACS)](https://github.com/epigen-bioinfolab/CHROMACS/releases)

  </div>
</div>



## 📖 About

ChromAcS (***Chrom***atin ***Ac***cessibility Analysis ***S***uite) provides a streamlined, GUI-based workflow for analyzing ATAC-seq data. Shown below are the major analysis steps which have been incorporated with several related dependencies in a hierarchical and modular manner designed for the best user experience and easy navigation:

- Quality Control (`FastQC, MultiQC`) ==> of both raw reads and trimmed reads;

- Trimming (`Trim Galore`) ==> optional (users can opt for working with raw data or their own trimmed data saved in "raw data directory");

- Genome Alignment and Coverage Analysis(`bowtie2, SAMtools, deepTools`) ==> flexible with several model organisms and beyond, (uses Ensembl toplevel reference);

- Peak Calling (`Genrich or MACS3`) ==> employes two different peak calling tools, users can select any of the two, with several options for parameter setup;

- Peak Annotation (`ChIPseeker`) ==> builds own txDB by fetching GTF files from ENSEMBL, uses default orgdb from bioconductor;

- Differential Peak Analysis (`DiffBind or NOISeq`) ==> employs Diffbind to visualize differential studies and also generates essential up- and downgraded regions; includes annotation of those regions as well; in absence of biological replicates, users can opt for NOISeq which allows the same differential studies (however, if you have more than one biological replicates, it is recommended to use Diffbind)

- Motif Enrichment Analysis (`MEME Suite (FIMO)`) ==> performs transcription factor motif enrichment on differential peak sets using a user-provided .meme file (e.g., from JASPAR [https://jaspar.elixir.no/]). For significant motif hits, it further annotates predicted motif-to-gene associations, linking enriched motifs to nearby genes.

- Motif Footprinting (`TOBIAS`) ==> uses TOBIAS to process both single- or multiple replicate samples allowing footprinting across two different conditions. It involves signal correction (ATACorrect), footprint score calculation, differential footprint detection (BINDetect), along with visualization of aggregate plots of selected motifs with condition-specific support. 

- Peak Overlap and Expression Data Overlap  (`BEDTools`) ==> provided as a separate add-on (chromacs-addon), this utility overlaps user-supplied BED files with ChromAcS peak results to identify shared genomic regions. Also integrates expression data (e.g., RNA-seq) by matching identifiers, enabling the user to relate chromatin accessibility peaks to gene expression profiles.


Designed to assist both beginners and experienced users in analyzing chromatin accessibility with minimal command-line usage

## 📖 Documentation

[Download ChromAcS User Manual (PDF)](https://github.com/epigen-bioinfolab/CHROMACS/blob/main/User_Manual.pdf)


## ⚠️ Pre-requisites

Before installing ChromAcS, ensure you have the following:

- **For Windows users:** [Windows Subsystem for Linux (Ubuntu)](https://learn.microsoft.com/en-us/windows/wsl/install)
  
- **For Linux users:** ChromAcS Tested on [Ubuntu](https://ubuntu.com/download) and [Fedora](https://getfedora.org/); Other Linux/Mac users may need to install dependent tools from the source. For Mac users, an automated installer will be released soon
  
- **Git (optional):** [Install Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
  
- **Conda:** [Install Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

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
The folder will be named CHROMACS-main. (Using Linux terminal command) Configure to the CHROMACS-main
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

#### Step 4a (Add-on): Launch the ChromAcS-AddOn: Additional Analysis Toolkit :
```bash
chromacs-addon
```

## 📬 Contact

For any queries and support, please reach us at **epigen.bioinfolab@gmail.com**  
Visit our lab page: [www.epigen-bioinfolab.com](https://www.epigen-bioinfolab.com/)

## 🧩 Acknowledgements

- FastQC, MultiQC, Trim Galore, bowtie2, MACS3, Genrich, ChIPseeker, DiffBind, NOISeq, deepTools, MEME Suite, TOBIAS; along with their dependencies
- Bioconda, Conda-Forge community, Python Software Foundation and the Python community.

## 📝 License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
