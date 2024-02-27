#' Run fgsea analysis on XLSX
#' @param XLSX character 'all' or 'equal' to specify which data to load
#' @param TAG character to tag output and name plots etc
#' @rdname run_analysis
#' @export

run_analysis <- function(XLSX, TAG){
  data_list <- read_xlsx_glaucoma(XLSX)
  msigdb_pathlist <- msigdb_pathways_to_list()
  tb_stats <- data_list$tb$avg_log2FC
  names(tb_stats) <- data_list$tb$X1
  fgseaRes <- fgsea::fgsea(pathways = msigdb_pathlist,
                           stats    = tb_stats,
                           minSize  = 15,
                           maxSize  = 500)
  fgseaResTidy <- dplyr::arrange(.data = tibble::as_tibble(fgseaRes),
                                 desc(NES))
  msigdb_pway_xlsx(fgseaRes, data_list$tb, padj = 0.05, TAG)
  ggplot2::ggplot(fgseaResTidy, ggplot2::aes(reorder(pathway, NES), NES)) +
    ggplot2::geom_col(ggplot2::aes(fill=padj<0.05)) +
    ggplot2::coord_flip() +
    ggplot2::labs(x="Pathway", y="Normalized Enrichment Score",
         title=paste0(TAG, " - Hallmark pathways NES from GSEA")) +
    ggplot2::theme_minimal()
  filenam <- paste0(TAG, ".", gsub(".xlsx", ".pathway.pdf", data_list$filename))
  ggplot2::ggsave(filename = filenam)
}

#' Read data from XLSX in extdata
#' @param XLSX character filename
#' @rdname read_xlsx_glaucoma
#' @export

read_xlsx_glaucoma <- function(XLSX) {
  RXLSX <- system.file(package = "cobrienGlaucoma", "extdata",
                       XLSX)

  xlsx_tb <- tibble::as_tibble(openxlsx::read.xlsx(RXLSX))
  return(list(tb = xlsx_tb,
              filename = basename(RXLSX)))
}

#' msigdb pathways in a nice list
#'
#' @param msigdb_species character one of msigdbr::msigdbr_show_species()
#' @param msigdb_cat character one of 'c("H", paste0("C", c(1:7)))',
#'        see: gsea-msigdb.org/gsea/msigdb/collections.jsp
#' @return msigdb_pathlist list object
#' @export

msigdb_pathways_to_list <- function(msigdb_species = "Homo sapiens", msigdb_cat = "H"){
  msigdb_pathway <- msigdbr::msigdbr(species = msigdb_species, category = msigdb_cat)
  ##create list
  msigdb_pathlist <- lapply(unique(msigdb_pathway$gs_name), function(f){
    fsig <- dplyr::filter(.data = msigdb_pathway, gs_name %in% !!f)
    return(as.vector(unlist(dplyr::select(.data = fsig, gene_symbol))))
  })
  names(msigdb_pathlist) <- unique(msigdb_pathway$gs_name)
  return(msigdb_pathlist)
}

#' msigdb pathways genes parsed from XLSX tb
#'
#' @param fgseaRes tibble of fgseaResults
#' @param tb tibble of XLSX results
#' @param padj numeric value for padj below which we return results
#' @param TAG character to tag output and name plots etc
#' @return NULL, writes XLSX containing pathways as sheets, and genes in those significant pathway
#' @export

msigdb_pway_xlsx <- function(fgseaRes, tb, padj, TAG){
  wb <- openxlsx::createWorkbook(TAG)
  olist <- lapply(seq_along(fgseaRes$padj), function(f){
    if(fgseaRes$padj[f] < padj){
      pway_tb <- dplyr::filter(.data = tb, X1 %in% unlist(fgseaRes[f,]$leadingEdge))
      ll <- strsplit(gsub("HALLMARK_", "", unlist(fgseaRes[f,"pathway"])), "")[[1]]
      if(length(ll) > 31){
        name <- paste(ll[1:31], collapse = "")
      } else {
        name <- gsub("HALLMARK_", "", unlist(fgseaRes[f,"pathway"]))
      }
      openxlsx::addWorksheet(wb, sheetName = name)
      openxlsx::writeData(wb = wb, sheet = name, x = pway_tb)
      openxlsx::saveWorkbook(wb = wb,
                             file = paste0(TAG, ".sig_", padj, ".fgsea_pathway_genes.xlsx"),
                             overwrite = TRUE)
    }
  })
}
