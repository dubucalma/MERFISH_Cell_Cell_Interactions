## Cell-cell Interactions Analyses

The aim of this repository is to quantify spatial colocalization between annotated single cells as a proxy for potential functional interactions within local cortical microcircuits. By testing observed spatial patterns against randomized null models, we identify pairs of cell types that are significantly enriched or depleted in physical proximity.


*Preliminary:*<br>
To generate the cell types' distribution for the Figure 2D:
`Layers_content_analysis.R`

___

### 1. Processing the interactions

#### 1a. 1000 permutations of the cells positions / Randomization of the positions

Original script to compute permutations on the whole slice:
`Compute_frequency_randomization.R`

Derived script for permutations per layer (Upper and Deeper):
`Compute_frequency_randomization_layers.R`

Derived script for permutations per bin (electrophy-bin):
`Compute_frequency_randomization_electrophy_bins.R`


**Input requirements**  
Provide a metadata file describing the cells of each sample.  

The data frame should include the following columns/info:  
- `cell_id` : unique identifier of the cell (optional if rownames are unique)  
- `predicted.id` : predicted cell type label  
- `x_final` : x-coordinate  
- `y_final` : y-coordinate  
- `layer` : layer annotation (e.g., `"upper"`, `"deeper"`) for permutation tests  

The previous scripts compute the following actions: 

**Step 1: Randomization**<br>
Randomize 1000 times the cell's positions, within a radius of 100 micrometers. 
Each cell’s original (x,y) coordinates are jittered **within a circle of radius 100 µm**.
Purpose: preserve local density while breaking real spatial structure.

**Step 2: Interaction frequency computation**<br>
Compute the cell types interactions frequencies for each permutation.
Two cells are considered neighbors if they are **≤15 µm** apart.

**Output:**<br>
A serialized R object `(.rds)` containing a list of 1,000 cell–cell contact frequency tables (one per permutation).
Each table summarizes the number of neighbor relationships between all pairs of cell types.


In the script `Compute_frequency_randomization_layers.R` the randomization is done separately for each bin ("upper", "deeper"):
instead of one global set of permutations, `Compute_frequency_randomization()` is called once per bin, so `Randomization_per_bin` stores a list of results
split by layer (upper vs deeper).
<br>
<br>

#### 1b. Calculate the enrichement of specific cell-types interactions:

Scripts to do this analysis: 
Original CCI file (adapted a bit for bins): 
`CC_interactions_clean.R` ON PROCESS OF CLEANING

Derived CCI file for layers:
`CCI_clean_layers.R` (Check input data for randomization)

**Input requirements** <br>
This step uses as input:

- The permutation results from Part 1 (1000 permutations): `1000_Permuted_mean_frequencies.rds`
- The observed contact frequency tables (from `Observed_frequency.R`)

**Processing**<br>
Aggregate permutation results to derive the null distribution for each cell-type pair
- `Calculate_mean_matrice()` → mean expected frequencies
- `Calculate_sd_matrice() → standard deviation (null variance)

Compare observed frequencies against the null distribution
- Compute Fold Change, Z-score, and P-values
- Adjust for multiple testing (FDR)

**Output**
- Main output: a full results table with Observed, Expected mean, Fold Change, Z-score, P-value, and FDR for each cell-type pair.

- Optional simplified output: a compact file containing only the key enrichment metrics (e.g. Fold Change, P-value) for easier downstream visualization or plotting.

___
### 2. Downstream analyses

#### 2a. Correlation of CCI profiles between samples

File to generate the Pearson Correlation Heatmap:
`CCI_upper_deeper_correlations_clean.ipynb`

This analysis computes the **similarity of cell-cell interactions (CCI) fold change profiles**, across multiple brain samples, for specific subclasses of cell types.

**Input requirement** <br>
The fold change table obtained from the previous step (startfied by layer, e.g. upper vs deeper).

For each sample, there are two input files:
- fold_change_upper.tsv --> interactions in upper layers (e.g. L2/3)
- fold_change_deeper.tsv --> interactions in deeper layers (e.g L4-6)

**Filter by subclasses of interest**<br>
The data are filtered with `process_data_subclasses`, keeping only the interaction pairs that involve specific subclasses (e.g. Vip vs inhibitory).

**Processing** <br>
Compute **Pearson Correlations** between fold change vectors for all sample pairs.

**Output**<br>
A correlation matrix representing the similarity of CCI fold change profiles across samples, as well as a heatmap summarizing these results.
<br>
<br>


#### 2b. Generation of cluster plots with the fold change 

File to generate clusterplots and heatmaps with the fold change per pair (MGE versus CGE) in aSTG and pSTG slices:
`CCI_upper_deeper_barplots_clean.ipynb`

This notebook `CCI_upper_deeper_brplots_clean.ipynb` generates barplots, and heatmaps summarizing cell–cell interaction (CCI) fold-changes across samples. You can run it for CGE↔CGE or MGE↔MGE pairs and for upper (e.g., L2/3) or deeper (e.g., L4–L6) layers.

**Input requirements**<br>
The fold change table obtained from the previous step (stratified by layer, e.g. upper vs deeper).

For each sample, there are two input files:
- fold_change_upper.tsv → interactions in upper layers (e.g. L2/3)
- fold_change_deeper.tsv → interactions in deeper layers (e.g. L4–L6)

**Output**<br>
Barplots of standardized fold-change (Fold_change_z) per pair across samples for aSTG and pSTG.

The notebook also includes an optional second visualization step that produces heatmaps of the Fold change values per pair (raw or standardized), for aSTG and pSTG. Only positive fold changes are annotated inside the heatmap cells, and missing values across samples are displayed as light grey.
