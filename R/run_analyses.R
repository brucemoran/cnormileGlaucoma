#' Run fgsea analysis on XLSX
#' @param XLSX character filename
#' @param TAG character output naming scheme
#' @rdname run_analyses
#' @export

run_analyses <- function(){
  run_analysis(XLSX = "all_glaucoma_vs_all_control.xlsx",
               TAG = "All_Glaucoma_vs_All_Control")
  run_analysis(XLSX = "glaucoma_equal_sample_vs_all_control.xlsx",
              TAG = "Equal_Glaucoma_vs_All_Control")
  run_analysis(XLSX = "dividing_MSC_run.xlsx",
               TAG = "Dividing_MSC_Run")
  run_analysis(XLSX = "group4_COMBINED.xlsx",
               TAG = "Group4_Combined")
}
