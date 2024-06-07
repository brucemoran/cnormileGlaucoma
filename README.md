# cnormileGlaucoma R package

This R package contains functions used to analyse the data contained in 'extdata'.

The PDF and XLSX files produced show significant pathways and the genes comprising them.

To install the repo: `devtools::install_github("brucemoran/cnormileGlaucoma")`

To run analyses: `cnormileGlaucoma::run_analyses()`

## GSEA

The `{msigdbr}` package function `msigdbr` was used to with `species = "Homo sapiens"` and `category = "H"` to generate a pathway list object, `msigdb_pathlist`.

Data was read from XLSX files in the `extdata` directory. This data includes gene name and `avg_log2FC` columns. The `{fgsea}` function `fgsea` was run using `pathways = msigdb_pathlist`,
`stats = avg_log2FC` as per the XLSX files which were sorted in decreasing rank as required, `minSize = 15` and `maxSize = 500`. Results were sorted by descending normalised enrichment score (`NES`) value. For the `group_markers_2024.xlsx` analysis, `cluster` column was used to define the gene set and perform GSEA.

From these results, XLSX output was written for pathways with p-value < 0.05 for further interpretation. Plots per input XLSX were saved to PDF showing all analysed pathways, with indication of significant p-value pathways. The `xlsx` and `pdf` directories of the package hold these files. The `ALL_HALLMARK_PATHWAYS.xlsx` file in `xlsx` is the result of `write_all_hallmarks(TAG, OUTDIR)` function to write XLSX of Hallmark Pathways.
