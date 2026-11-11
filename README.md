# Organic amendments drive bacterial and fungal soil microbiomes over time in a subtropical banana agroecosystem

[![Manuscript status: under revision](https://img.shields.io/badge/manuscript-under%20revision-orange)](#citation)
[![Raw reads: NCBI SRA PRJNA1368706](https://img.shields.io/badge/SRA-PRJNA1368706-blue)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1368706)
[![Made with R](https://img.shields.io/badge/R-4.5-276DC3?logo=r)](https://www.r-project.org/)

This repository contains the data and R code that produce all statistical analyses, figures and supplementary tables of:

> Correa-Delgado R, Romero-Lozano I, Reinhardt S, Jaizme-Vega MC, Laich F, Bodenhausen N.
> **Organic amendments drive bacterial and fungal soil microbiomes over time in a subtropical banana agroecosystem.**
> *Environmental Microbiome* (under revision).

The scripts start from the ZOTU count tables and sample metadata and end with the figure panels in `04_Figures/` and the Word tables in `04_Tables/`.

---

## Contents

- [Study overview](#study-overview)
- [Repository structure](#repository-structure)
- [Data](#data)
- [Requirements](#requirements)
- [How to reproduce the analyses](#how-to-reproduce-the-analyses)
- [Figures and tables in the manuscript](#figures-and-tables-in-the-manuscript)
- [Upstream bioinformatics (not in this repository)](#upstream-bioinformatics-not-in-this-repository)
- [Citation](#citation)
- [Contributors and contact](#contributors-and-contact)
- [License](#license)

---

## Study overview

Banana plants (*Musa acuminata* AAA cv. 'Gruesa Palmera') were grown in 170-L pots of Andisol soil at the Canary Islands Institute of Agricultural Research (ICIA, Tenerife, Spain) from August 2021 to June 2023. Five treatments were compared:

| Code in data | Treatment |
|---|---|
| `control` | Unamended control |
| `chicken` | Chicken manure |
| `cow` | Cow manure |
| `compost` | Compost |
| `pellet` | Pelletised compost |

Soil was sampled at six time points (T0 = 0, T1 = 2.4, T2 = 4.8, T3 = 9.3, T4 = 13.7 and T5 = 20.8 months). At each time point there were three composite samples per treatment (replicates A, B and C), for **90 samples** in total. Bacterial (16S rRNA gene, 27F/1492R) and fungal (ITS region, ITS1F/ITS4) communities were sequenced with PacBio long-read sequencing.

The analyses in this repository cover:

1. **Sequencing quality and community composition**: sequencing depth, rarefaction curves and phylum-level relative abundance.
2. **Alpha diversity**: richness, Shannon and evenness, modelled with linear mixed-effects models (`lme4`/`lmerTest`) and estimated marginal means (`emmeans`).
3. **Beta diversity**: CLR transformation, PCA, PERMANOVA (`vegan::adonis2`) with permutations restricted within pots, and redundancy analysis (RDA) against soil chemistry (`ordiR2step`) and plant variables (`envfit`).
4. **Differential abundance**: DESeq2 (with `apeglm` shrinkage) and ANCOM-BC2, each amendment against the control. Only OTUs detected by both methods are reported.
5. **Culture-dependent vs. sequencing data**: correlations between CFU counts and relative abundances from sequencing.

The soil chemistry, plant and CFU data come from the same experiment, which is described in [Correa-Delgado et al. (2026), *J Soil Sci Plant Nutr*](https://doi.org/10.1007/s42729-025-02979-8).

---

## Repository structure

```
.
├── 01_Raw_data/              Input data (see "Data" below)
│   └── dada2_trial/          DADA2 taxonomy assignment of ITS ZOTUs against UNITE
├── 02_Scripts/               R Markdown / R scripts, numbered in run order
├── 03_Analysis_files/        Intermediate results: saved model outputs (DESeq2, ANCOM-BC2),
│   │                         taxonomy tables, ggplot objects (.rds/.RDa), merged tables (.csv)
│   └── Supplementary_plot/   ggplot objects for the supplementary figures
├── 04_Figures/               Figure panels (.png) written by the scripts
├── 04_Tables/                Supplementary tables written by the scripts (.docx)
├── 05_word_documents/        Manuscript and additional files (revised version)
└── Banana_clean.Rproj        RStudio project file (sets the project root)
```

Most scripts come in pairs, one for bacteria (`*_16S`) and one for fungi (`*_ITS`):

| Script | Content |
|---|---|
| `01_Introduction_{16S,ITS}.Rmd` | Data import, removal of singletons, sequencing depth, rarefaction, phylum-level composition |
| `02_Alpha_{16S,ITS}.Rmd` | Richness, Shannon and evenness; linear mixed-effects models |
| `03_Beta_{16S,ITS}.Rmd` | CLR, PCA, PERMANOVA, RDA with soil chemistry (all time points) and plant variables (T5) |
| `04_DAA_16S.Rmd`, `04_DAA_ITSRmd.Rmd` | DESeq2 and ANCOM-BC2 models, volcano plots |
| `04b_DAA_{16S,ITS}_merged_dotplot.R` | Merges the DESeq2 and ANCOM-BC2 results and draws a family-level dot plot of OTUs selected by both methods |
| `04c_DAA_results_tables.R` | Exports the differential-abundance tables to Word |
| `05_Correlations_{16S,ITS}.Rmd` | CFU counts vs. relative abundance from sequencing; *Trichoderma* |

---

## Data

### Raw sequencing reads

Raw PacBio reads are deposited in the NCBI Sequence Read Archive under BioProject **[PRJNA1368706](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1368706)**. They are not needed to run the scripts in this repository.

### Files in `01_Raw_data/`

| File | Description |
|---|---|
| `p1031_run231218_Revio_16S_ZOTU_c97_Count_Sintax.txt` | 16S ZOTU (97 %) count table with SINTAX taxonomy, QIIME format |
| `p1031_run231218_Revio_16S_ZOTU_c97.fa` | 16S representative sequences |
| `p1031_run231023_ITS_ccs_ZOTU_c97_Count_Sintax.txt` | ITS ZOTU (97 %) count table with SINTAX taxonomy, QIIME format |
| `p1031_run231023_ITS_ccs_ZOTU_c97.fa` | ITS representative sequences |
| `metadata_2025.csv` | Sample metadata: `sample_name`, `time` (T0–T5), `treatment`, `rep` (A–C), `sample_ID` |
| `metadata_2025_soil.xlsx` | Soil chemistry per sample (pH, EC, TOM, OOM, N, P, Ca, Mg, Na, K) |
| `plant_data.xlsx` | Plant biometric and nutrient data at final harvest (T5) |
| `cfu_bacteria.xlsx`, `cfu_actinobacteria.xlsx`, `cfu_fungi.xlsx`, `cfu_trichoderma.xlsx` | Culture-dependent counts (CFU g⁻¹ soil) |
| `dada2_trial/` | UNITE reference (`sh_general_release_dynamic_19.02.2025.fasta`), the DADA2 assignment script and its output |

The DADA2 taxonomy tables used for differential abundance (`03_Analysis_files/taxonomy_{16S,ITS}_dada2.RDa`) are precomputed and included in the repository.

---

## Requirements

The analyses were run in **R 4.5.3** with **Bioconductor 3.22**, in RStudio. Main packages:

| Purpose | Packages |
|---|---|
| Microbiome data handling | `phyloseq`, `microbiome` |
| Community ecology | `vegan` |
| Mixed models and post-hoc tests | `lme4`, `lmerTest`, `emmeans`, `car`, `multcomp`, `rstatix` |
| Differential abundance | `DESeq2` (1.50.2), `apeglm`, `ANCOMBC` (2.x) |
| Data wrangling and plotting | `tidyverse`, `ggpubr`, `ggrepel`, `gridExtra`, `scales`, `readxl` |
| Word tables | `flextable`, `officer` |
| Taxonomy assignment (optional) | `dada2`, `Biostrings` |

To install them:

```r
install.packages(c("tidyverse", "vegan", "lme4", "lmerTest", "emmeans", "car",
                   "multcomp", "rstatix", "ggpubr", "ggrepel", "gridExtra",
                   "scales", "readxl", "flextable", "officer", "BiocManager"))

BiocManager::install(c("phyloseq", "microbiome", "DESeq2", "apeglm",
                       "ANCOMBC", "dada2", "Biostrings"))
```

The ANCOM-BC2 models in `04_DAA_*` are the slowest step. Their results are saved in `03_Analysis_files/`, so the `04b`/`04c` scripts can be re-run without re-fitting them. Random seeds are set in the scripts that use permutations or resampling (`set.seed(123)` in `03_Beta_*`, `set.seed(5070)` in `04_DAA_*`).

---

## How to reproduce the analyses

1. Clone the repository and open `Banana_clean.Rproj` in RStudio:

   ```bash
   git clone https://github.com/FiBL/banana.git
   ```

2. File paths in the scripts are relative to the **project root** (the folder that contains `Banana_clean.Rproj`). In RStudio, set *Knit → Knit Directory → Project Directory* before knitting, or render from the R console:

   ```r
   rmarkdown::render("02_Scripts/02_Alpha_16S.Rmd", knit_root_dir = getwd())
   ```

3. Run the scripts in numerical order. The 16S and ITS scripts are independent of each other. Within each marker, the order matters only for the differential-abundance scripts:

   ```
   01_Introduction  →  02_Alpha  →  03_Beta  →  05_Correlations      (independent of each other)

   04_DAA_16S.Rmd / 04_DAA_ITSRmd.Rmd       saves DESeq2 + ANCOM-BC2 results in 03_Analysis_files/
            ↓
   04b_DAA_{16S,ITS}_merged_dotplot.R       merged table + Figure 3 dot plots
            ↓
   04c_DAA_results_tables.R                 Tables S7 and S8 (Word)
   ```

   The `.R` scripts can be run with `source()` from the project root or from `02_Scripts/`.

Each script writes its figures to `04_Figures/` (PNG) and the corresponding ggplot objects to `03_Analysis_files/` (`.rds`/`.RDa`), so panels can be restyled without re-running the models. The multi-panel figures in the manuscript were assembled from these panels.

---

## Figures and tables in the manuscript

### Main figures

| Figure | Content | Script | Output in `04_Figures/` |
|---|---|---|---|
| Fig. 1 | OTU richness over time (A bacteria, B fungi) | `02_Alpha_16S.Rmd`, `02_Alpha_ITS.Rmd` | `Figure_1_16S.png`, `Figure_1_ITS.png` |
| Fig. 2 | RDA constrained by soil chemistry | `03_Beta_16S.Rmd`, `03_Beta_ITS.Rmd` | `Figure_2_16S_leyenda.png`, `Figure_2_ITS.png` |
| Fig. 3 | Differential abundance by family, DESeq2 and ANCOM-BC2 | `04b_DAA_16S_merged_dotplot.R`, `04b_DAA_ITS_merged_dotplot.R` | `DAA_16S_family_dotplot.png`, `DAA_ITS_family_dotplot.png` |
| Fig. 4 | CFU vs. relative abundance (A aerobic bacteria, B Ascomycota) | `05_Correlations_16S.Rmd`, `05_Correlations_ITS.Rmd` | `Figure_4_Bacteria.png`, `Figure_4B_Fungi.png` |

### Additional file 2: supplementary figures and tables

| Item | Content | Script | Output |
|---|---|---|---|
| Fig. S1 | Sequencing depth per sample | `01_Introduction_*` | `SupplementaryFigure1_{16S,ITS}.png` |
| Fig. S2 | Rarefaction curves | `01_Introduction_*` | `SupplementaryFigure2_{16S,ITS}.png` |
| Fig. S3 | Phylum-level relative abundance | `01_Introduction_*` | `SupplementaryFigure3_{16S,ITS}_leyenda.png` |
| Fig. S4 | Shannon diversity and evenness | `02_Alpha_*` | `Supplemetary_FigureS4*.png` |
| Fig. S5 | PCA of CLR-transformed data and PERMANOVA | `03_Beta_*` | `Supplemetary_FigureS5_{16S,ITS}.png` |
| Fig. S6 | RDA constrained by plant variables (T5) | `03_Beta_*` | `Supplemetary_FigureS6_{16S,ITS}.png` |
| Fig. S7 | Volcano plots, DESeq2 and ANCOM-BC2 | `04_DAA_16S.Rmd`, `04_DAA_ITSRmd.Rmd` | `FigureX_16S_revisado_2_deseq.png`, `FigureXX_16S_revisado_2_ancom.png`, `FigureX_ITS_deseq_revisado_2.png`, `FigureXXX_ITS_ancom_revisado_2.png` |
| Fig. S8 | Actinobacteria and *Trichoderma*: CFU vs. relative abundance | `05_Correlations_16S.Rmd`, `05_Correlations_ITS.Rmd` | `Supplemetary_FigureS7_Actinobacteria.png`, `Figure_XB_trichoderma.png`, `Figure_X_trichoderma.png` |
| Table S1 | Chemistry of the organic amendments | — (from the companion study) | — |
| Table S2 | Mixed-model results for alpha diversity | `02_Alpha_*` | printed in the knitted report |
| Tables S3–S4 | Forward selection of soil variables for the RDA | `03_Beta_*` | printed in the knitted report |
| Tables S5–S6 | `envfit` of plant variables | `03_Beta_*` | printed in the knitted report |
| Tables S7–S8 | Differentially abundant OTUs (bacteria, fungi) | `04c_DAA_results_tables.R` | `04_Tables/DAA_16S_table.docx`, `04_Tables/DAA_ITS_table.docx` |

Some output file names are older than the final figure numbering of the manuscript (for example `Supplemetary_FigureS7_Actinobacteria.png` is now part of Fig. S8). The table above gives the current mapping.

---

## Upstream bioinformatics (not in this repository)

Reads were processed before the analyses in this repository, on the ETH Zurich Euler cluster:

1. Demultiplexing with [lima](https://lima.how) v3.1.1
2. Primer trimming and chimera removal with USEARCH v11.0.667
3. Denoising into zero-radius OTUs (ZOTUs) with UNOISE3, then clustering at 97 % identity with UPARSE
4. Taxonomy with SINTAX against SILVA SSU v138 (16S) and with DADA2 `assignTaxonomy()` against UNITE (ITS; see `01_Raw_data/dada2_trial/`)

Sequencing was done at the Next Generation Sequencing Platform of the University of Bern: Revio for 16S, Sequel IIe for ITS.

---

## Citation

If you use this code or data, please cite the article. The full reference and DOI will be added here once the article is published.

> Correa-Delgado R, Romero-Lozano I, Reinhardt S, Jaizme-Vega MC, Laich F, Bodenhausen N. Organic amendments drive bacterial and fungal soil microbiomes over time in a subtropical banana agroecosystem. *Environmental Microbiome* (under revision).

The experimental design, soil chemistry, plant and CFU data are described in:

> Correa-Delgado R, Bodenhausen N, Romero I, Brito-López P, Jaizme Vega MC, Laich F. Impact of different organic amendments on soil health and banana plant performance. *Journal of Soil Science and Plant Nutrition* (2026). https://doi.org/10.1007/s42729-025-02979-8

---

## Contributors and contact

Code by Ignacio Romero-Lozano (FiBL) and Raquel Correa-Delgado (ICIA / EEZ-CSIC), with Natacha Bodenhausen (FiBL).

Corresponding authors: Raquel Correa-Delgado (raquel.correa@eez.csic.es) and Natacha Bodenhausen (natacha.bodenhausen@fibl.org).

For questions about the code, please [open an issue](https://github.com/FiBL/banana/issues).

---

## License

No license has been chosen yet. Until one is added, the default copyright applies and reuse requires permission from the authors.
