# ReMiTT
R codes for ReMiTT
Please do the analysis following the steps in analyze_ovarian1_git.R
Before that, please run codes in Tcells_lung_cancer.R to create seurat object for single cell data of the T cells and save it to the working directory.\\
Then, the first step is to download the spatial transcriptomics data and prepare it for analysis with the function "data_preparation"\\
Then, please check the data with plots and figures and define tumor regions by high expression or cancer markers\\
Then, please use function "get_all_paths" to get potential trails with increasing exhaustion score\\
Lastly, use function "path_selection" to select trails with high correlation of genes related to T cell stage\\
