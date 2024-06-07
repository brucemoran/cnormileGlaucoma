#' Run fgsea analysis on XLSX
#' @param XLSX character 'all' or 'equal' to specify which data to load
#' @param TAG character to tag output and name plots etc
#' @param OUTDIR character where output is saved
#' @rdname run_analysis
#' @export

run_analysis <- function(XLSX, TAG, OUTDIR){
  print(paste0(" - ", TAG))
  data_list <- read_xlsx_glaucoma(XLSX)
  msigdb_pathlist <- msigdb_pathways_to_list()
  tb_stats <- data_list$tb$avg_log2FC
  colnames(data_list$tb)[1] <- "gene"
  names(tb_stats) <- data_list$tb$gene
  fgseaRes <- fgsea::fgsea(pathways = msigdb_pathlist,
                           stats    = tb_stats,
                           minSize  = 15,
                           maxSize  = 500)
  fgseaResTidy <- dplyr::arrange(.data = tibble::as_tibble(fgseaRes),
                                 desc(NES))
  msigdb_pway_xlsx(fgseaRes, data_list$tb, padj = 0.05, TAG, OUTDIR)
  ggplot2::ggplot(fgseaResTidy, ggplot2::aes(reorder(pathway, NES), NES)) +
    ggplot2::geom_col(ggplot2::aes(fill=padj<0.05)) +
    ggplot2::coord_flip() +
    ggplot2::labs(x="Pathway", y="Normalized Enrichment Score",
         title=paste0(TAG, " - Hallmark pathways NES from GSEA")) +
    ggplot2::theme_minimal()
  dir.create(paste0(OUTDIR, "/pdf/"), recursive = TRUE, showWarnings = FALSE)
  filenam <- paste0(OUTDIR, "/pdf/", TAG, ".", gsub(".xlsx", ".pathway.pdf", data_list$filename))
  ggplot2::ggsave(filename = filenam)
}

#' Run fgsea analysis on grouped XLSX
#' @param XLSX character 'all' or 'equal' to specify which data to load
#' @param TAG character to tag output and name plots etc
#' @param OUTDIR character where output is saved
#' @rdname run_group_analysis
#' @export

run_group_analysis <- function(XLSX, TAG, OUTDIR){
  print(paste0(" - ", TAG))
  data_list <- read_xlsx_glaucoma(XLSX)
  msigdb_pathlist <- msigdb_pathways_to_list()
  ##split on cluster
  res_list <- lapply(unique(unlist(data_list$tb$cluster)), function(ff){
    TAG <- paste0(TAG, "_Cluster_", ff)
    tb_ff <- dplyr::filter(.data = data_list$tb, cluster %in% ff)
    tb_stats <- tb_ff$avg_log2FC
    names(tb_stats) <- tb_ff$gene
    fgseaRes <- fgsea::fgsea(pathways = msigdb_pathlist,
                             stats    = tb_stats,
                             minSize  = 15,
                             maxSize  = 500)
    fgseaResTidy <- dplyr::arrange(.data = tibble::as_tibble(fgseaRes),
                                   desc(NES))
    msigdb_pway_xlsx(fgseaRes, tb_ff, padj = 0.05, TAG, OUTDIR)
    ggplot2::ggplot(fgseaResTidy, ggplot2::aes(reorder(pathway, NES), NES)) +
      ggplot2::geom_col(ggplot2::aes(fill = padj<0.05)) +
      ggplot2::coord_flip() +
      ggplot2::labs(x="Pathway", y="Normalized Enrichment Score",
                    title=paste0(TAG, " - Hallmark pathways NES from GSEA")) +
      ggplot2::theme_minimal()
    dir.create(paste0(OUTDIR, "/pdf/"), recursive = TRUE, showWarnings = FALSE)
    filenam <- paste0(OUTDIR, "/pdf/", TAG, ".", gsub(".xlsx", ".pathway.pdf", data_list$filename))
    ggplot2::ggsave(filename = filenam)
  })
}

#' Read data from XLSX in extdata
#' @param XLSX character filename
#' @rdname read_xlsx_glaucoma
#' @export

read_xlsx_glaucoma <- function(XLSX) {
  RXLSX <- system.file(package = "cnormileGlaucoma", "extdata",
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
#' @param OUTDIR character where output is saved
#' @return NULL, writes XLSX containing pathways as sheets, and genes in those significant pathway
#' @export

msigdb_pway_xlsx <- function(fgseaRes, tb, padj, TAG, OUTDIR){
  wb <- openxlsx::createWorkbook(TAG)
  olist <- lapply(seq_along(fgseaRes$padj), function(f){
    if(fgseaRes$padj[f] < padj){
      pway_tb <- dplyr::filter(.data = tb, gene %in% unlist(fgseaRes[f,]$leadingEdge))
      ll <- strsplit(gsub("HALLMARK_", "", unlist(fgseaRes[f,"pathway"])), "")[[1]]
      if(length(ll) > 31){
        name <- paste(ll[1:31], collapse = "")
      } else {
        name <- gsub("HALLMARK_", "", unlist(fgseaRes[f,"pathway"]))
      }
      openxlsx::addWorksheet(wb, sheetName = name)
      openxlsx::writeData(wb = wb, sheet = name, x = pway_tb)
      dir.create(paste0(OUTDIR, "/xlsx/"), recursive = TRUE, showWarnings = FALSE)
      openxlsx::saveWorkbook(wb = wb,
                             file = paste0(OUTDIR, "/xlsx/", TAG, ".sig_", padj, ".fgsea_pathway_genes.xlsx"),
                             overwrite = TRUE)
    }
  })
}

#' msigdb pathways genes to XLSX
#'
#' @param TAG character to tag output
#' @param OUTDIR character where output is saved
#' @return NULL, writes XLSX containing pathways as sheets, and genes in those significant pathway
#' @export

write_all_hallmarks <- function(TAG = "All_Hallmark_Pathways", OUTDIR){
  msigdb_pathlist <- msigdb_pathways_to_list()
  wb <- openxlsx::createWorkbook(TAG)
  lapply(names(msigdb_pathlist), function(name){
    ssn <- strsplit(name, "")[[1]]
    if(length(ssn) > 30){
      sname <- paste(strsplit(name, "")[[1]][1:30], collapse = "")
    } else {
      sname <- name
    }
    openxlsx::addWorksheet(wb, sheetName = sname)
    openxlsx::writeData(wb = wb, sheet = sname, x = msigdb_pathlist[[name]])
    openxlsx::saveWorkbook(wb = wb,
                           file = paste0(paste0(OUTDIR, "/xlsx/"), TAG, ".xlsx"),
                           overwrite = TRUE)
  })
}

