# INVA

All raw data and metadata are stored in the "raw data" folder. Below is a breakdown of the data and code files provided, along with instructions on how they can be used to reproduce the analyses and figures from the manuscript.

## Data Overview

1. **FT-ICR-MS Data**: The FT-ICR-MS data and associated metadata can be directly processed using the [MetaboAnalyst pipeline](https://metabodirect.readthedocs.io/en/latest/).
2. **LCMS, NMR, and MALDI Data**: These datasets require pre-processing before analysis.
   - **Chemical Properties**: Stored in files named `NMR_chem` or `compound_link`.
   - **Feature Intensities**: Stored in files named as `[environment]_[mode]`, such as `Rhizo_LCMS_hneg.csv`. Here, `[environment]` indicates the sample type (e.g., "Rhizo"), and `[mode]` specifies the method and/or ionization mode.
   - **Linking Chemical Properties and Intensities**: Chemical properties and intensity data are linked by chemical names or unique compound IDs, such as `[FT_m/z ratio_retention time_mode]`.

## Code and Analysis Files

1.  **0_optional_LCMS_cleaning.r**: This script was used to remove gap-filled peaks in LC MS/MS data to create supplementary data - table 6
   
3.  **1_peak_overlapping.r**: This script was used to create Figures 1a and 1e.
   - The source file `peak_comp_root/rhizo_neutral.csv` is a dataframe containing the neutral mass of all features detected by each method, with methods as column names.
   - The neutral mass of LCMS features was inferred by adjusting the actual weight of the adduct by ±1.00728 for negative/positive modes, respectively.

3. **MetaboDirect Results**:
   - Results from MetaboDirect were used to generate several figures:
     - Figures 1b and 1f (from MetaboDirect output: `5_statistics, 1.4`).
     - Figure 2f (from MetaboDirect output: `3_exploratory, comparison-plant`).
     - Figures 4a and 4e (from MetaboDirect output: `4_chemodiversity, 2.2`).
     - Figure 2a and 2e were modified from `3_exploratory, 6_Composition_by_element`.
     - Figure 3a was modified from `3_exploratory, 6_Composition_by_class`.
     - Figure 4b and 4f were modified from MetaboDirect output: `3_exploratory, comparison-plant, 1_vk_plant_all_data`.
   - Remaining information is compiled in Supplementary Note 1.

4. **Ordination and Differential Analysis Scripts**:
   - `2_NMR_rhizo.r`, `2_NMR_root.r`, `3_LCMS_rhizo.r`, `3_LCMS_root.r`: These scripts follow the same structure and perform ordination analysis, plotting NMDS plots, conducting differential analysis, creating volcano plots, and generating differential tables. Most results appear in Supplementary Figure 1, Supplementary Notes 2 and 3. R² values from NMDS analyses were used to create Figures 4c and 4g.

5. **Differential Analysis Merging**:
   - `4_differential_analysis.R`: This script merges differential analysis results from LCMS and NMR for further analysis and includes figures 3b, 3f, and 3h.

6. **Boxplot Generation for LC-MS/MS**:
   - `5_differential_plotting.R`: This script generates all boxplots from LC-MS/MS data.

7. **Pathway Plotting**:
   - `6_pathway_plotting.R`: Used to plot Figure 2b by displaying the average normalized intensity of multiple modes.

8. **MALDI-MSI Analysis**:
   - `7_MALDI.R`: Conducts differential analysis for MALDI-MSI data and generates differential tables. Only boxplots (Figures 2d and 3c) are included in the main text.

9. **Enrichment Analysis**:
   - `8_ChemRichFunction.r`: This script performs enrichment analysis for Figures 2g and 4d.

