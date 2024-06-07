#' Run fgsea analysis on XLSX
#' @param XLSX character filename
#' @param TAG character output naming scheme
#' @rdname run_analyses
#' @export

run_analyses <- function(){
  print("Running analyses...")
  OUTDIR <- "~/Downloads/cnormileGlaucoma"
  dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
  run_analysis(XLSX = "all_glaucoma_vs_all_control.xlsx",
               TAG = "All_Glaucoma_vs_All_Control",
               OUTDIR = OUTDIR)
  run_analysis(XLSX = "glaucoma_equal_sample_vs_all_control.xlsx",
              TAG = "Equal_Glaucoma_vs_All_Control",
              OUTDIR = OUTDIR)
  run_analysis(XLSX = "dividing_MSC_run.xlsx",
               TAG = "Dividing_MSC_Run",
               OUTDIR = OUTDIR)
  run_analysis(XLSX = "group4_COMBINED.xlsx",
               TAG = "Group4_Combined",
               OUTDIR = OUTDIR)
  run_group_analysis(XLSX = "group_markers_2024.xlsx",
               TAG = "Group_Markers",
               OUTDIR = OUTDIR)
}
