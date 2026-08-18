# climate-wine-metabolomic

# A taste of climate change: Warming rewrites wine chemistry across time and geography

## Data and Code Availability
This repository contains the R scripts used to perform the statistical analyses and generate the figures for the manuscript. The raw and processed datasets required to run these scripts are provided as Supplementary Data files (`Data S1.xlsx`, `Data S2.xlsx`, `Data S3.xlsx`, and `Data S4.xlsx`).

To reproduce the analysis, download the Supplementary Data files and ensure they are located in your working directory alongside the respective R scripts.

---

### Figure 2
**Associated Data File:** `Data S1.xlsx`
* **`figure2A_climate.R`** -> Reads sheet **Fig 2A** (Historical climate trends)
* **`figure2B_C13C12.R`** -> Reads sheet **Fig 2B** (Stable carbon isotope ratios)
* **`figure2C_pH,anthocyanin&6forest.R`** -> Reads sheet **Fig 2C** (Physicochemical attributes)
* **`figure2D_e-nose,tongueheatmap.R`** -> Reads sheets **Fig 2D_e-nose** and **Fig 2D_e-tongue** (Sensory digitization heatmaps)
* **`figure2E_only PCA.R`** & **`figure2E_only distance(dot).R`** -> Reads sheets **Fig 2E_GCMS**, **Fig 2E_headspaceGCMS**, **Fig 2E_LCp**, and **Fig 2E_LCn** (Untargeted metabolomic variance and PC1 distance)
* **`figure2F_vplot.R`** -> Reads sheets **Fig 2F__histogram_edit** and **Fig 2F_metabolite_correlation** (Temporal chemical drivers and volcano plots)
* **`figure2G_bubble.R`** -> Reads sheet **Fig 2G** (Sensory-chemical integrative correlation)

---

### Figure 3
**Associated Data File:** `Data S2.xlsx`
* **`figure3A_tem,HI.R`** -> Reads sheets **Fig 3A_map** and **Fig 3A_tem,HI** (Regional climate profiles and Huglin Index)
* **`figure3B_PCA,PCvsHI.R`** -> Reads sheet **Fig 3B** (Global metabolic variance and PCA score plots)
* **`figure3C_significant metabolite.R`** & **`figure3C_vplot.R`** -> Reads sheet **Fig 3C** (Heat-responsive global biomarkers)
* **`figure3D_distance.R`** -> Reads sheets **Fig 3D_0** through **Fig 3D_5** (Temperature-dependent simulated regional trajectories)

---

### Figure 4
**Associated Data File:** `Data S3.xlsx`
* **`figure4B_metabolite.R`** & **`figure4B_vplot.R`** -> Reads sheet **Fig 4B** (Longitudinal metabolic volatility and forecasting)
* **`figure4C_16metabolite barplot.R`** -> Reads sheet **Fig 4C** (Precursor dose-response shifts)
* **`figure4D_heatmap(change dot plot).R`** -> Reads sheet **Fig 4D** (Downstream effect sizes/slopes)
*(Note: Panel 4E predictive cross-validation utilizes sheet **Fig 4E**)*

---

### Figure 5 & Supplementary Fig. S1
**Associated Data File:** `Data S4.xlsx`
* **`figureS1B_markers.R`** -> Reads sheets **FigS1B_Past**, **FigS1B_Present**, and **FigS1B_Future** (Cross-timeline thermal sensitivity swarm plots)
* **`figureS1CD.R`** -> Reads sheets **FigS1C** and **FigS1D** (Tripartite climatic/spatial validation matrix)
* **`figure5B_variables.R`** -> Reads sheets **Fig5B_past**, **Fig5B_present**, and **Fig5B_future** (Categorized biological class trajectories)
* **`figure5D_taste.R`** -> Reads sheet **Fig5D** (Targeted gustatory and sensory shift distribution)
