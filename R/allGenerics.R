#### subsetting ----

#' Subsetting a ProMetIS \code{SummarizedExperiment}, \code{MultiAssayExperiment}, \code{ExpressionSet} or \code{MultiDataSet}
#'
#' Subsetting ProMetIS datasets according to genes, sex and tissues, and filtering out
#' features with too many NAs, too low variance, or too many imputed values (proteomics
#' datasets only)
#'
#' @param x An S4 object of class \code{SummarizedExperiment}, \code{MultiAssayExperiment}, \code{ExpressionSet} or \code{MultiDataSet}
#' @param set.c Character: name of the data set
#' @param genes.vc Character vector: with elements in 'LAT', 'MX2', and 'WT'; when
#' set to 'all', the 'c('WT', 'LAT', 'MX2')' vector will be used
#' @param sex.vc Character vector: with elements in 'M' and 'F'; when
#' set to 'all', the 'c('M', 'F')' vector will be used
#' @param tissues.vc Character vector: with elements in 'liver' and 'plasma'; when
#' set to 'all', the 'c('liver', 'plasma')' vector will be used
#' @param common_samples.l Logical: should the datasets be restricted to common samples?
#' @param na_thresh.n Numeric: maximal proportion of NAs for a feature to be kept
#' @param var_thresh.n Numteric: minimal variance for a feature to be kept
#' @param imputed_thresh.n Numeric: for proteomics datasets, the features with a
#' too high proportion of imputed values in all conditions to be compared will be
#' discarded 
#' @return \code{ExpressionSet} or \code{MultiDataSet} with the selected sets, samples,
#' and features (after filtering for the proportion of NAs, the minimum of variance,
#' and the proportion of imputed values for proteomics datasets)
#' @rdname subsetting
#' @export
#' @examples
#' # MultiAssayExperiment
#' prometis.mae <- phenomis::reading(ProMetIS::post_processed_dir.c(),
#'                                 report.c = "none")
#' prometis.mae <- prometis.mae[, , ProMetIS::sets.vc()]
#' plasma_mx2.mae <- ProMetIS::subsetting(prometis.mae,
#'                                        genes.vc = c("WT", "MX2"),
#'                                        tissues.vc = "plasma")
#' # SummarizedExperiment
#' set.c <- "proteomics_liver"
#' proteo.mae <- phenomis::reading(ProMetIS::post_processed_dir.c(),
#'                                 subsets.vc = set.c,
#'                                 report.c = "none")
#' proteo.se <- proteo.mae[[set.c]]
#' ## adding (common) sample metadata to the summarized experiment
#' sample_meta.DF <- SummarizedExperiment::colData(proteo.mae)
#' SummarizedExperiment::colData(proteo.se)[, "gene"] <- sample_meta.DF[colnames(proteo.se), "gene"]
#' SummarizedExperiment::colData(proteo.se)[, "sex"] <- sample_meta.DF[colnames(proteo.se), "sex"]
#' proteo_mx2.se <- ProMetIS::subsetting(proteo.se,
#'                                       set.c = set.c,
#'                                       genes.vc = c("WT", "MX2"))
#' # MultiDataSet
#' prometis.mds <- phenomis::reading(ProMetIS::post_processed_dir.c(),
#'                                 output.c = "set",
#'                                 report.c = "none")
#' prometis.mds <- prometis.mds[, ProMetIS::sets.vc()]
#' plasma_mx2.mds <- ProMetIS::subsetting(prometis.mds,
#'                                        genes.vc = c("WT", "MX2"),
#'                                        tissues.vc = "plasma")
#' # ExpressionSet
#' set.c <- "proteomics_liver"
#' proteo.mset <- phenomis::reading(ProMetIS::post_processed_dir.c(),
#'                                  subsets.vc = set.c,
#'                                  output.c = "set",
#'                                  report.c = "none")
#' proteo.set <- proteo.mset[[set.c]]
#' proteo_mx2.set <- ProMetIS::subsetting(proteo.set,
#'                                        set.c = set.c,
#'                                        genes.vc = c("WT", "MX2"))
setGeneric("subsetting",
           function(x,
                    set.c = NULL,
                    genes.vc = "all",
                    sex.vc = "all",
                    tissues.vc = "all",
                    common_samples.l = FALSE,
                    na_thresh.n = 0.2,
                    var_thresh.n = .Machine$double.eps,
                    imputed_thresh.n = 0.2)
             standardGeneric("subsetting"))


#### filter_overimputed ####

#' Filtering overimputation of proteomics datasets
#'
#' Filtering overimputation of proteomics datasets
#'
#' @param x An S4 object of class \code{SummarizedExperiment} or \code{ExpressionSet}
#' @param set.c Character: name of the data set
#' @param genes.vc Character vector: with elements in 'LAT', 'MX2', and 'WT'
#' @param sex.vc Character vector: with elements in 'M' and 'F'
#' @param imputed_thresh.n Numeric: for proteomics datasets, the features with a
#' too high proportion of imputed values in all conditions to be compared will be
#' discarded 
#' @return logical vector corresponding to the selected features
#' @rdname filter_overimputed
#' @examples
#' set.c <- "proteomics_plasma"
#' 
#' # SummarizedExperiment
#' 
#' proteo.mae <- phenomis::reading(ProMetIS::post_processed_dir.c(),
#'                                 subsets.vc = set.c)
#' proteo.se <- proteo.mae[[set.c]]
#' ## adding (common) sample metadata to the summarized experiment
#' sample_meta.DF <- SummarizedExperiment::colData(proteo.mae)
#' SummarizedExperiment::colData(proteo.se)[, "gene"] <- sample_meta.DF[colnames(proteo.se), "gene"]
#' SummarizedExperiment::colData(proteo.se)[, "sex"] <- sample_meta.DF[colnames(proteo.se), "sex"]
#' feat_select.vl <- ProMetIS:::filter_overimputed(proteo.se, set.c = set.c, genes.vc = "MX2")
#' message(round(sum(feat_select.vl) / length(feat_select.vl) * 100), "% of selected features")
#' 
#' # ExpressionSet
#' 
#' proteo.mset <- phenomis::reading(ProMetIS::post_processed_dir.c(),
#'                                  subsets.vc = set.c,
#'                                  output.c = "set")
#' proteo.eset <- proteo.mset[[set.c]]
#' feat_select.vl <- ProMetIS:::filter_overimputed(proteo.eset, set.c = set.c, genes.vc = "MX2")
setGeneric("filter_overimputed",
           function(x,
                    set.c,
                    genes.vc = "all",
                    sex.vc = "all",
                    imputed_thresh.n = 0.2)
             standardGeneric("filter_overimputed"))


#### imputation_info ####

#' Checking overimputation in proteomics datasets
#'
#' Checking overimputation in proteomics datasets
#'
#' @param x An S4 object of class \code{SummarizedExperiment} or \code{ExpressionSet}
#' @param set.c Character: name of the data set
#' @return integer matrix indicating the imputed values in the intensity matrix
#' @rdname imputation_info
#' @examples
#' set.c <- "proteomics_plasma"
#' 
#' # SummarizedExperiment
#' 
#' proteo.mae <- phenomis::reading(ProMetIS::post_processed_dir.c(),
#'                                 subsets.vc = set.c)
#' proteo.se <- proteo.mae[[set.c]]
#' imputed.mi <- ProMetIS:::imputation_info(proteo.se, set.c = set.c)
#' # computing the number of imputed values per feature and per genotype
#' genotype.fc <- factor(substr(colnames(imputed.mi), 1, 1), levels = c("L", "X", "W"))
#' geno_impute.mn <- t(apply(imputed.mi, 1, function(feat.vn) tapply(feat.vn, genotype.fc, sum)))
#' colnames(geno_impute.mn) <- c("LAT", "MX2", "WT")
#' colSums(geno_impute.mn)
#' geno_impute.df <- tidyr::gather(as.data.frame(geno_impute.mn), genotype, imputation, LAT:WT, factor_key = TRUE)
#' ggplot2::ggplot(geno_impute.df, ggplot2::aes(x = imputation, colour = genotype)) + ggplot2::geom_density()
#' 
#' # ExpressionSet
#' 
#' proteo.mset <- phenomis::reading(ProMetIS::post_processed_dir.c(),
#'                                  subsets.vc = set.c,
#'                                  output.c = "set")
#' proteo.eset <- proteo.mset[[set.c]]
#' imputed.mi <- ProMetIS:::imputation_info(proteo.eset, set.c = set.c)
setGeneric("imputation_info",
           function(x,
                    set.c)
             standardGeneric("imputation_info"))
