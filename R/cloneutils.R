
#' Compute clonealign fits across a range of highly variable genes
#'
#' @export
compute_cal_fits <- function(hvg_props,
                             sce,
                             cnv_data) {
  rvs <- rowVars(as.matrix(logcounts(sce)))
  ca_reps <- lapply(hvg_props, function(rvq) {
    genes_to_sample <- rvs > quantile(rvs, rvq)
    em <- clonealign(sce[genes_to_sample,], cnv_data[genes_to_sample,], rel_tol_em = 1e-8)
    em
  })
  ca_reps
}

#' Construct a precision recall curve using the output of compute_cal_fits
#'
#' @export
construct_pr_df <- function(cal_fits, base_fit, hvg_props) {
  calculate_precision_recall <- function(base_fit, new_fit, hvgp) {
    mat <- table(new_fit$clone, base_fit$clone)
    precision <- diag(mat) / rowSums(mat)
    recall <- (diag(mat) / colSums(mat))
    data_frame(precision = precision, recall = recall, clone = names(precision), hvg_prop = hvgp)
  }

  dfs <- lapply(seq_along(hvg_props), function(i) calculate_precision_recall(base_fit, cal_fits[[i]], hvg_props[i]))
  df_pr <- bind_rows(dfs) %>%
    gather(measure, value, -clone, -hvg_prop)
  df_pr
}

#' Plot the precision recall curve output by construct_pr_df
#'
#' @export
plot_pr <- function(df_pr) {
  ggplot(df_pr, aes(x = hvg_prop, y = value, color = clone, linetype = measure, shape = measure)) +
    geom_point() +
    geom_line() +
    ggsci::scale_color_d3() +
    labs(x = "Quantile highly variable genes", y = "Measure value")
}

#' Subset a SingleCellExperiment to have genes with a predicted dosage effect
#' 
#' The SingleCellExperiment should have \code{Symbol} in rowData
#' dosage df should at least have columns "gene" and "q.value"
#' 
#' The returned SingleCellExperiment contains only genes whose
#' q.value is less than 0.05
#' 
#' @export
dosage_subset_sce <- function(sce, dosage_df, threshold = 0.05) {
  genes_to_subset <- dplyr::filter(dosage_df, q.value < threshold)$gene
  sce[rowData(sce)$Symbol %in% genes_to_subset]
}

