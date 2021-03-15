#### subsetting ####

#' Subsetting a ProMetIS \code{ExpressionSet} or \code{MultiDataSet}
#'
#' Subsetting ProMetIS datasets according to genes, sex and tissues, and filtering out
#' features with too many NAs, too low variance, or too many imputed values (proteomics
#' datasets only)
#'
#' @param x An S4 object of class \code{ExpressionSet} or \code{MultiDataSet}
#' @param set.c Character: name of the \code{ExpressionSet}
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
#' latmx2.mset <- phenomis::reading(ProMetIS::statistics_singleomics_dir.c(),
#'                                  report.c = "none")
#' plasma_lat.mset <- ProMetIS::subsetting(latmx2.mset,
#'                                         genes.vc = c("WT", "LAT"),
#'                                         tissues.vc = "plasma")
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