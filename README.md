````md
# Mfuzz Time-course RNA-seq Pipeline

A simple, configurable R pipeline for **time-course bulk RNA-seq clustering** using **Mfuzz**.

This repository is designed for users who want to:

- cluster genes by temporal expression pattern
- compare multiple color palettes for Mfuzz plots
- export both **PDF** and **cluster-wise PNG** figures
- extract **alpha-core genes**
- reuse the same pipeline for other datasets by editing only `config.yaml`

This project is especially useful for normalized bulk RNA-seq expression matrices such as:

- `Ctrl / 1hr / 4hr`
- `0h / 2h / 8h / 24h`
- any time-course design with replicate columns grouped by time point

---

# Table of Contents

- [1. What this pipeline does](#1-what-this-pipeline-does)
- [2. Repository structure](#2-repository-structure)
- [3. Requirements](#3-requirements)
- [4. Input file format](#4-input-file-format)
- [5. Quick start for beginners](#5-quick-start-for-beginners)
- [6. Running the pipeline](#6-running-the-pipeline)
- [7. How `config.yaml` works](#7-how-configyaml-works)
- [8. Important configuration sections](#8-important-configuration-sections)
- [9. Output files](#9-output-files)
- [10. Recommended starting settings](#10-recommended-starting-settings)
- [11. Troubleshooting](#11-troubleshooting)
- [12. Notes and interpretation](#12-notes-and-interpretation)
- [13. Reusing this pipeline for other datasets](#13-reusing-this-pipeline-for-other-datasets)
- [14. Future improvements](#14-future-improvements)
- [15. License](#15-license)

---

# 1. What this pipeline does

This pipeline takes a normalized RNA-seq expression matrix and performs the following steps:

1. Reads an input CSV file
2. Groups replicate columns by time point
3. Aggregates replicates using **mean** or **median**
4. Optionally applies `log2(x + pseudo_count)`
5. Handles missing values using:
   - `filter.NA`
   - `fill.NA`
6. Removes low-variance genes using `filter.std`
7. Standardizes each gene profile using `standardise`
8. Estimates the Mfuzz fuzzifier `m` using `mestimate()` or uses a user-defined value
9. Runs `mfuzz()` clustering
10. Exports:
    - processed matrices
    - cluster assignments
    - membership matrix
    - alpha-core genes
    - full PDF plots
    - per-cluster PNG plots
    - run summary files

---

# 2. Repository structure

A minimal repository looks like this:

```text
mfuzz-timecourse-pipeline/
├─ README.md
├─ config.yaml
├─ mfuzz_pipeline_csv.R
└─ input.csv
````

After running the pipeline, an output directory will be created.

Example:

```text
mfuzz-timecourse-pipeline/
├─ README.md
├─ config.yaml
├─ mfuzz_pipeline_csv.R
├─ input.csv
└─ mfuzz_out/
   ├─ 01_group_matrix_used_for_mfuzz.csv
   ├─ 02_after_na_handling.csv
   ├─ 03_filter_std_plot.pdf                  # if enabled
   ├─ 04_after_sd_filter.csv
   ├─ 05_zscore_or_standardised_matrix.csv
   ├─ 06_mestimate.csv
   ├─ 06A_cselection_plot.pdf                 # if enabled
   ├─ 06A_cselection_matrix.csv               # if enabled
   ├─ 06B_Dmin_plot.pdf                       # if enabled
   ├─ 06B_Dmin_values.csv                     # if enabled
   ├─ 07_gene_cluster_assignment.csv
   ├─ 08_alpha_core_cluster_1.csv
   ├─ 08_alpha_core_cluster_2.csv
   ├─ ...
   ├─ 08_alpha_core_cluster_total.csv
   ├─ 09_mfuzz_plot_all_genes__viridis.pdf
   ├─ 09_mfuzz_plot_all_genes__viridis__cluster_1.png
   ├─ 09_mfuzz_plot_all_genes__viridis__cluster_2.png
   ├─ ...
   ├─ 10_mfuzz_plot_core_genes__viridis.pdf
   ├─ 10_mfuzz_plot_core_genes__viridis__cluster_1.png
   ├─ 10_mfuzz_plot_core_genes__viridis__cluster_2.png
   ├─ ...
   ├─ 11_membership_matrix.csv
   ├─ 12_cluster_centers.csv
   ├─ 13_cluster_size_summary.csv
   ├─ 14_run_summary.csv
   └─ 15_sessionInfo.txt
```

---

# 3. Requirements

You need:

* **R**
* **RStudio** (recommended, but not required)

The script will try to install the required packages automatically.

## CRAN packages

* `yaml`
* `viridisLite`

## Bioconductor packages

* `Biobase`
* `Mfuzz`

---

# 4. Input file format

The input file must be a **CSV** file.

## Required conditions

* The first row must be a header
* A `gene_id` column is required
* A `symbol` column is recommended
* Expression values should be in the remaining columns
* Missing values can be written as `NA`

## Example input

```csv
gene_id,symbol,HPctrl1,HPctrl2,HPctrl3,HPctrl4,HP1hr1,HP1hr2,HP1hr3,HP1hr4,HP4hr1,HP4hr2,HP4hr3,HP4hr4
ENSMUSG00002076977,ENSMUSG00002076977,NA,NA,1.17116204,2.941738695,0.981122523,2.307445942,NA,0.837087612,NA,4.303119708,NA,1.672421109
ENSMUSG00002076974,ENSMUSG00002076974,2.779497652,NA,NA,1.96115913,NA,NA,1.125316253,0.837087612,0.925231191,5.378899635,NA,2.508631664
ENSMUSG00002076968,ENSMUSG00002076968,2.779497652,5.586472927,3.513486119,3.92231826,0.981122523,2.307445942,2.250632506,0.837087612,8.327080721,9.682019343,5.356473929,3.344842218
```

---

# 5. Quick start for beginners

If this is your first time using R, follow these exact steps.

## Step 1. Create a folder

Create a folder anywhere on your computer, for example:

```text
Desktop/RNAseq_Mfuzz
```

## Step 2. Put these files inside the folder

* `config.yaml`
* `mfuzz_pipeline_csv.R`
* `input.csv`

Your folder should look like this:

```text
RNAseq_Mfuzz/
├─ config.yaml
├─ mfuzz_pipeline_csv.R
└─ input.csv
```

## Step 3. Open RStudio

Launch RStudio.

## Step 4. Open the R script

Open `mfuzz_pipeline_csv.R` in RStudio.

## Step 5. Set the working directory

In the RStudio top menu, click:

**Session → Set Working Directory → To Source File Location**

This makes sure R looks for `config.yaml` and `input.csv` in the same folder as the script.

## Step 6. Run the script

Click the **Source** button in the top-right corner of the script editor.

The script is written to automatically run using `config.yaml`.

## Step 7. Check the results

After the run finishes, a results folder will be created.

By default, it will be:

```text
mfuzz_out/
```

inside your project folder.

---

# 6. Running the pipeline

## Option A. Run from RStudio

Open `mfuzz_pipeline_csv.R`, set the working directory to the script location, then click **Source**.

## Option B. Run from the R console

```r
source("mfuzz_pipeline_csv.R")
```

## Option C. Run from a terminal

```bash
Rscript mfuzz_pipeline_csv.R config.yaml
```

---

# 7. How `config.yaml` works

`config.yaml` is the central settings file for the pipeline.

You usually do **not** need to edit the R script.

You can control the following using `config.yaml`:

* input file path
* `gene_id` and `symbol` column names
* time-point sample groups
* mean or median aggregation
* optional log transform
* missing value filtering and filling
* low-variance filtering
* number of clusters
* color palettes
* PDF figure size
* PNG image resolution
* output folder name

---

# 8. Important configuration sections

## 8.1 Input file path

```yaml
data:
  input_file: "input.csv"
```

This tells the pipeline which file to read.

---

## 8.2 Gene annotation columns

```yaml
data:
  gene_id_col: "gene_id"
  symbol_col: "symbol"
```

These define which columns contain the gene identifier and gene symbol.

---

## 8.3 Time-point sample groups

```yaml
data:
  sample_groups:
    - name: "Ctrl"
      time: 0
      columns: ["HPctrl1", "HPctrl2", "HPctrl3", "HPctrl4"]

    - name: "1hr"
      time: 1
      columns: ["HP1hr1", "HP1hr2", "HP1hr3", "HP1hr4"]

    - name: "4hr"
      time: 4
      columns: ["HP4hr1", "HP4hr2", "HP4hr3", "HP4hr4"]
```

This is the **most important section**.

Each group needs:

* `name`: label shown in plots
* `time`: actual numeric time value
* `columns`: replicate columns for that time point

If you reuse this pipeline for another dataset, this is usually the first section you need to edit.

---

## 8.4 How replicates are combined

```yaml
preprocess:
  aggregate_fun: "mean"
```

Options:

* `mean`
* `median`

For normalized bulk RNA-seq data, `mean` is usually the best starting choice.

---

## 8.5 Optional log transform

```yaml
preprocess:
  apply_log2: false
  pseudo_count: 1
```

* Set `apply_log2: false` if your values are already log-scaled
* Set `apply_log2: true` if your values are normalized but not yet log-transformed

As a rough rule:

* values around `0–15` often suggest log-scale data
* values in the hundreds or thousands often suggest non-log normalized values

---

## 8.6 Missing value handling

```yaml
preprocess:
  filter_NA:
    enabled: true
    thres: 0.34

  fill_NA:
    enabled: true
    mode: "mean"
    k: 10
```

### `filter.NA`

Removes genes with too many missing values.

### `fill.NA`

Fills remaining missing values.

Supported modes:

* `mean`
* `median`
* `knn`
* `knnw`

---

## 8.7 Low-variance gene filtering

```yaml
preprocess:
  filter_std:
    enabled: true
    min_std: 0.5
    visu: false
```

This removes genes that change very little across time.

* larger `min_std` = stricter filtering
* smaller `min_std` = more genes kept, but possibly more noise

A good starting value is usually `0.5`.

---

## 8.8 Number of clusters

```yaml
clustering:
  cluster_number: 4
```

This is the `c` value for Mfuzz.

For a dataset with only 3 time points, starting with **4 clusters** is often a good practical choice.

---

## 8.9 Color palette settings

```yaml
plot:
  render_multiple_palettes: true
  palette_names_to_render: ["viridis", "magma", "inferno", "hot_cold"]
```

This renders the same plots using multiple palettes so you can compare them.

Supported palette names:

* `default`
* `fancy`
* `viridis`
* `magma`
* `inferno`
* `plasma`
* `cividis`
* `hot_cold`

If you want only one palette:

```yaml
plot:
  render_multiple_palettes: false
  active_palette: "viridis"
```

---

## 8.10 PNG resolution settings

```yaml
output:
  png:
    save_cluster_png: true
    width: 1800
    height: 1400
    res: 200
    bg: "white"
```

These control the cluster-wise PNG outputs.

* `save_cluster_png`: whether to save cluster-wise PNG files
* `width`, `height`: image size in pixels
* `res`: image resolution in dpi
* `bg`: background color

For publication or presentation figures, you may want:

```yaml
output:
  png:
    save_cluster_png: true
    width: 2400
    height: 1800
    res: 300
    bg: "white"
```

---

# 9. Output files

## 9.1 Intermediate processed matrices

### `01_group_matrix_used_for_mfuzz.csv`

Replicate-aggregated matrix used as the starting input to Mfuzz.

### `02_after_na_handling.csv`

Matrix after missing-value filtering and filling.

### `04_after_sd_filter.csv`

Matrix after low-variance gene filtering.

### `05_zscore_or_standardised_matrix.csv`

The standardized matrix actually used for clustering.

---

## 9.2 Clustering results

### `06_mestimate.csv`

The estimated fuzzifier `m`.

### `07_gene_cluster_assignment.csv`

The main cluster assignment table.

Columns:

* `gene_id`
* `symbol`
* `cluster`
* `max_membership`

### `11_membership_matrix.csv`

The full Mfuzz membership matrix.

### `12_cluster_centers.csv`

Cluster center profiles.

### `13_cluster_size_summary.csv`

Number of genes assigned to each cluster.

### `14_run_summary.csv`

A summary of the run, including the number of genes used and the cluster number.

---

## 9.3 Alpha-core files

### `08_alpha_core_cluster_1.csv`, `08_alpha_core_cluster_2.csv`, ...

Cluster-specific alpha-core gene lists.

### `08_alpha_core_cluster_total.csv`

A row-wise combined table of all alpha-core genes across all clusters.

Columns:

* `gene_id`
* `symbol`
* `MEM.SHIP`
* `cluster_id`

This is usually the most convenient file if you want one combined alpha-core summary.

---

## 9.4 Plot outputs

### Full PDF plots

* `09_mfuzz_plot_all_genes__<palette>.pdf`
* `10_mfuzz_plot_core_genes__<palette>.pdf`

Examples:

* `09_mfuzz_plot_all_genes__viridis.pdf`
* `10_mfuzz_plot_core_genes__magma.pdf`

### Cluster-wise PNG plots

* `09_mfuzz_plot_all_genes__<palette>__cluster_1.png`
* `09_mfuzz_plot_all_genes__<palette>__cluster_2.png`
* ...
* `10_mfuzz_plot_core_genes__<palette>__cluster_1.png`
* `10_mfuzz_plot_core_genes__<palette>__cluster_2.png`
* ...

These are useful for presentations, figure assembly, and cluster-by-cluster inspection.

---

# 10. Recommended starting settings

For a 3-time-point bulk RNA-seq dataset, a good starting setup is:

```yaml
preprocess:
  apply_log2: false

  filter_std:
    enabled: true
    min_std: 0.5

clustering:
  cluster_number: 4

plot:
  render_multiple_palettes: true
  palette_names_to_render: ["viridis", "magma", "inferno", "hot_cold"]

output:
  png:
    save_cluster_png: true
    width: 1800
    height: 1400
    res: 200
```

---

# 11. Troubleshooting

## Problem 1. “Required columns” error

Cause:

* column names in `config.yaml` do not match the input CSV header exactly

Fix:

* check the header row in `input.csv`
* update:

  * `gene_id_col`
  * `symbol_col`
  * `sample_groups.columns`

---

## Problem 2. Too few genes remain after filtering

Cause:

* `filter.NA` is too strict
* `filter_std.min_std` is too high

Fix:

* increase `filter_NA.thres`
* reduce `filter_std.min_std`

Example:

```yaml
preprocess:
  filter_std:
    enabled: true
    min_std: 0.3
```

---

## Problem 3. The PNG plots look blurry

Cause:

* resolution is too low

Fix:

* increase `output.png.res`
* also increase `width` and `height`

Example:

```yaml
output:
  png:
    save_cluster_png: true
    width: 2400
    height: 1800
    res: 300
```

---

## Problem 4. Warning about incomplete final line in `config.yaml`

Cause:

* the file does not end with a newline

Fix:

* open `config.yaml`
* go to the last line
* press **Enter**
* save the file

This warning is not usually fatal, but it is better to fix it.

---

## Problem 5. Palette-related errors

This pipeline uses `viridisLite` for:

* `viridis`
* `magma`
* `inferno`
* `plasma`
* `cividis`

If you still run into palette-related issues, try rendering a single palette first:

```yaml
plot:
  render_multiple_palettes: false
  active_palette: "viridis"
```

---

## Problem 6. The script does not seem to run when clicking Source

Make sure you did all of the following:

1. `config.yaml`, `mfuzz_pipeline_csv.R`, and `input.csv` are in the same folder
2. You opened `mfuzz_pipeline_csv.R` in RStudio
3. You clicked
   **Session → Set Working Directory → To Source File Location**
4. Then clicked **Source**

---

# 12. Notes and interpretation

* This pipeline assumes **normalized expression values** as input
* It is **not** a DESeq2 raw-count differential expression workflow
* Be careful when converting `0` values to `NA`

  * real missing values should be `NA`
  * real zero expression values should remain `0`
* If the number of time points is small, avoid using too many clusters
* Cluster plots are most meaningful when interpreted together with:

  * `cluster centers`
  * `membership values`
  * `alpha-core genes`

---

# 13. Reusing this pipeline for other datasets

This repository is intended to be reusable.

In most cases, to analyze a different dataset, you only need to change:

* `data.input_file`
* `data.gene_id_col`
* `data.symbol_col`
* `data.sample_groups`

For example, if your columns are:

* `ctrl_A`, `ctrl_B`, `ctrl_C`
* `t1_rep1`, `t1_rep2`, `t1_rep3`
* `t4_rep1`, `t4_rep2`

you can change:

```yaml
data:
  sample_groups:
    - name: "Ctrl"
      time: 0
      columns: ["ctrl_A", "ctrl_B", "ctrl_C"]

    - name: "1hr"
      time: 1
      columns: ["t1_rep1", "t1_rep2", "t1_rep3"]

    - name: "4hr"
      time: 4
      columns: ["t4_rep1", "t4_rep2"]
```

That allows the same script to work without editing the R code.

---

# 14. Future improvements

Possible future extensions for this repository:

* automatic heatmap generation
* top-gene summaries for each cluster
* GO / pathway enrichment integration
* transcript-level support
* interactive visualization
* batch processing of multiple datasets

---

# 15. License

Add the license that fits your project.

Example:

```text
MIT License
```

Or:

```text
This repository is intended for academic/research use.
```

```

If you want, I can also draft a matching `.gitignore` and a short GitHub repository description.
```
