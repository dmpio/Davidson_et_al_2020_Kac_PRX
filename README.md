# Davidson_et_al_2019_Kac_PRX

Proteomic Data Analysis for the following publication:

> Davidson, M.T., Grimsrud, P.A., Lai, L., Draper, J.A., Fisher-Wellman, K.H., Narowski, T.M., Abraham, D.M., Koves, T.R., Kelly, D.P., and Muoio, D.M. (2019). Disruption of acetyl group balance in cardiomyocytes augments the mitochondrial acetylproteome without affecting respiratory function or heart susceptibility to pressure overload.  *In review* [link placeholder]

## Raw Data
The raw data for the proteomics experiments was submitted to the Proteome Xchange Consortium via the PRIDE partner Repository with the dataset identifier [PXD013935](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD013935)

## Directory Structure
The directory contains the analyses for 3 individual proteomic experiments. The structure of these folders is summarized below.

The 4th folder (`Comparisons between Experiments`) uses data processed from experiments 2 and 3 for additional analysis. 

## Structure for Experiments Folders

### `Exp#_..._Analysis.ipynb`
A jupyter notebook containing the python code used to analyze the data produced from Proteome Discoverer 2.2 for the specific experiment

### `pd22_exports`
Contains the tab-delimited data exported from Proteome Discoverer 2.2. One file contains the peptide information, the other contains protein information.

### `limma`
Contains the R code used to analyze the data with the *limma* library, along with the input for and the output from this analysis. The *limma* analysis is an intermediate step in processes detailed/performed in the `Exp#_..._Analysis.ipynb` file.

### `processed_files`
Contains the analysis results from `Exp#_..._Analysis.ipynb` in both Excel files and Comma Separated Value format. See below for more information

## Information about files within `processed_files` 

There are 2 Excel files within the directory:
    
1. The File ending with `plex` contains the processed data for the acetylpeptides corrected for the "relative occupancy" (abb. "RO"), acetylpeptides only corrected for loading content, and the proteins

2. The file ending with `Top25` contains processed data derived from the data in file #1 that shows the top 25 most differentially abundant acetylpeptides (corrected for the relative occupancy) along with basic information about the acetylated residues involved and the protein information
    
All of the sheets within the Excel files are also stored in ".csv" formatted files.

___

Within the first Excel file (and corresponding .csv files) for each experiment, a basic nomenclature is utilized. 

   There are 2 sheets for the acetylpeptide data. One contains the "relative occupancy" data which is subsequently analyzed within the *limma* software. The other doesn't contain the "relative occupancy data" and only analyzes the load-normalized data. In the paper, for consistency all numbers and graphs utilize the "relative occupancy" data. 

Nomenclature in the Acetylpeptide sheets: 
    
- columns that end in `_Acetyl` are the "raw" quantification values exported from Proteome Discoverer 2.2
- columns that end in `_norm` are the quantification values corrected for differences in loading between the individual channels
- columns that end in `_norm_log2` are simply the log2 transformation of the load-normalized columns
- columns that end in `_norm_log2_ro` have used the log2 transformation data of the **protein** (not peptide) to correct for any protein expression changes that may account for differences in the acetyl-peptide abundance. The `ro` stands for "relative occupancy". The relative occupancy of the post-translational modification is calculated by subtracting changes in the protein abundance from the acetyl-peptide abundance. For example, if a PTM site is found to be increasing by 2 fold, but the protein abundance is also increasing by 2 fold, then the relative occupancy would be 0.
- The results from the *limma* analysis start after the last `_ro` column until the end of each table

The Proteins Sheet follows a similar nomenclature. The difference is that the columns ending in `_Protein` are the "raw" quantification values exported from Proteome Discoverer 2.2


In the *limma* analysis, columns ending in `_significant` indicate the statistical significance of the specified comparison after multiple hypothesis correction with the Benjamini-Hochberg method and an FDR of <0.05.
    
- 0 indicates not significantly different
- 1 indicates significantly different, increased in the listed comparison
- -1 (negative one) indicates significantly different, decreased in the listed comparison
    
    
## Abbreviations:

- `DKO` Dual Knock-out. (Crat and Sirt3 ablation via cre-recombinase induced under the muscle creatine kinase promoter)
- `DFC`  Dual Flox control. (Floxed controls for the DKO animals)
- `S3KO`  Sirt3 Knock-out. (Sirt3 ablation via cre-recombinase induced under the muscle creatine kinase promoter)
- `S3FL` Sirt3 Flox control. (Floxed controls for the S3KO animals)
- `TACMI` Trans-aortic constriction + small apical myocardial infarction. (Surgical procedure to induce heart failure)
- `TAC` Trans-aortic constriction (Surgical procedure to induce heart failure)
- `FC` Fold Change
- `Log2FC` Fold Change in the Log2 space
- `Kac` Peptide containing an acetylated lysine residue (K=abbreviation for lysine)
- `RO` relative occupancy
