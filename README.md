## Analysis of Orthologous Collections (AOC) on REM2, REM2: a member of the RGK (Rem, Rem2, Rad, and Gem/Kir) family of small GTPases.

This repository contains all code, data, and results for the evolutionary analysis of mammalian REM2.

### Clone repository
`git clone https://github.com/aglucaci/AOC-REM2.git`
`cd AOC-REM2`

### Create environment (conda)
`conda env create -f environment.yml`

### Configuration (HPC-environment)
Set your cluster configurations in the `cluster.json` file.

### Data
The `data` folder contains data from NCBI. 

```
── data/mammalian_REM2
│   ├── Rem2_orthologs.csv
│   ├── Rem2_refseq_protein.fasta
│   ├── Rem2_refseq_transcript.fasta

```

### Run the analysis

When ready to use the Snakemake pipeline type, from the working directory:

`bash run_HPC_AOC.sh` 

### Results

The results of our analyses are populated in the `results/mammalian_REM2` folder.

### Scripts

This folder `scripts` holds a number of useful scripts used for our analyses.

### Notebooks

`executiveSummary.ipynb` this notebook creates a high-level overview of the selection analysis results.

`SummarizeBUSTED_ModelTest.ipynb` summaries the results from the BUSTED-family of analyses and outputs the table: `mammalian_REM2_BUSTED_ModelTest.csv` in the `tables` folder

`SummarizeFEL-FWER.ipynb` summaries the results from the FEL analysis, and outputs the following table: `mammalian_REM2_FEL_Results_FWER_adjusted_mapped.csv`, also performs multiple test correction (FWER) and maps alignment sites the Human reference sequence.

`SummarizeMEME-FWER.ipynb` summaries the results from the MEME analysis, and outputs the following table: `mammalian_REM2_MEME_Results_FWER_adjusted_mapped.csv`, also performs multiple test correction (FWER) and maps alignment sites the Human reference sequence.

`SummarizeBGM.ipynb` summarizes results from the BGM analysis. It creates visualizations and outputs a table to `tables`, `mammalian_REM2_BGM_Results_Mapped.csv`

`AlignmentMapper.ipynb` creates a mapping between the Human sites in the alignment (ungapped) and the gapped alignment. This is output to the `results folder as `mammalian_REM2_codons.SA.FilterOutliers.fasta_AlignmentMap.csv`

`AlignmentProfiler.ipynb` creates a number of useful figures and calculates descriptive statistics for a multiple sequence alignment.

`ViewNewick.ipynb` A useful notebook for examining Newick trees

`View_TN93_GeneticDistances.ipynb` A useful notebook for examining TN93 genetic distances.

`Model_Structure.ipynb` is used to map sites from MEME and FEL analyses to the Human references and creates an output table used for visualization in `tables`, `mammalian_REM2_StructureView.csv`

