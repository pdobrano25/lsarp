***Personalized resistant starch enriches butyrate-producing bacteria and upregulates mucosal mitochondrial proteins in therapeutically responsive children with new-onset Crohn's disease: A prospective, randomized, placebo-controlled trial.***

Peter Dobranowski
Version 2026_08_28

---

Mapping, pH, and processed meta-omics data are available in **[lsarp_supp_data](https://github.com/pdobrano25/lsarp/blob/main/lsarp_data)**

16S data are available under NCBI SRA ?.

Metaproteomics data are available under PRIDE ?.

Code and analyses for each section of the manuscript:

---

*Data processing*: **[2026_01_08_lsarp_dna_processing](https://github.com/pdobrano25/lsarp/blob/main/2026_01_08_lsarp_dna_processing.R)**


Prepares FFQ data (with minimal analyses), 16S data variables (e.g. butyrogens, predicted microbial load, functional redundancy), metagenomic data, metaproteomic data (e.g. proteins collapsed to functions), and metabolomic data. Output serves as input for most analysis scripts.

---

*Stool analyses*: **[2026_05_16_lsarp_analysis](https://github.com/pdobrano25/lsarp/blob/main/2026_05_16_lsarp_analysis_DI.R)**

Conducts major stool multi-omic analyses, including group-level, response group-level, and machine learning for CD participants.

---

*Mucosal analyses*: **[2026_01_08_lsarp_cd_proteomics_analysis](https://github.com/pdobrano25/lsarp/blob/main/2026_01_08_lsarp_cd_proteomics_analysis.R)**

Conducts host biopsy proteomics and mucosal bacteriome analyses for CD participants.

---

*RapidAIM analyses*: **[2026_01_08_lsarp_rapidaim_analysis](https://github.com/pdobrano25/lsarp/blob/main/2026_01_08_lsarp_rapidaim_test.R)**

Conducts analyses including RS selections and machine learning on RapidAIM data.

---


*UC analyses*: **[2026_05_18_lsarp_uc_analysis](https://github.com/pdobrano25/lsarp/blob/main/2026_05_18_lsarp_uc_analysis_DI.R)**

Conducts data processing and 16S analyses for UC participants.

---
