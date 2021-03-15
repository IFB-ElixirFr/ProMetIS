#' @export
genes.vc <- function() {
  c("LAT", "MX2")
}

#' @export
wtgenes.vc <- function() {
  c("WT", ProMetIS::genes.vc())
}

#' @export
palette.vc <- function(light.l = FALSE) {
  palette.df <- read.table(system.file("extdata/palette.tsv",
                                       package = "ProMetIS"),
                           header = TRUE,
                           row.names = 1,
                           quote = "",
                           sep = "\t",
                           stringsAsFactors = FALSE,
                           comment.char = "")
  
  if (light.l) {
    palette.vc <- palette.df[, "color_light"]
  } else
    palette.vc <- palette.df[, "color"]
  
  names(palette.vc) <- rownames(palette.df)
  
  palette.vc
  
}

#' @export
sets.vc <- function() {
  c("preclinical",
    
    "proteomics_liver",
    "proteomics_plasma",
    
    "metabolomics_liver_c18hypersil_pos",
    "metabolomics_liver_hilic_neg",

    "metabolomics_plasma_c18hypersil_pos",
    "metabolomics_plasma_hilic_neg",
    "metabolomics_plasma_c18acquity_pos",
    "metabolomics_plasma_c18acquity_neg")
}


#' @export
proteo_sets.vc <- function() {
  c("proteomics_liver",
    "proteomics_plasma")
}

#' @export
metabo_sets.vc <- function() {
  c("metabolomics_liver_c18hypersil_pos",
    "metabolomics_liver_hilic_neg",

    "metabolomics_plasma_c18hypersil_pos",
    "metabolomics_plasma_hilic_neg",
    "metabolomics_plasma_c18acquity_pos",
    "metabolomics_plasma_c18acquity_neg")
}

#' @export
liver_sets.vc <- function() {
  c("proteomics_liver",
    "metabolomics_liver_c18hypersil_pos",
    "metabolomics_liver_hilic_neg")
}

#' @export
plasma_sets.vc <- function() {
  c("proteomics_plasma",

    "metabolomics_plasma_c18hypersil_pos",
    "metabolomics_plasma_hilic_neg",
    "metabolomics_plasma_c18acquity_pos",
    "metabolomics_plasma_c18acquity_neg")
}

#' @export
sex.vc <- function() {
  c("M", "F")
}

#' @export
tissues.vc <- function() {
  c("liver", "plasma")
}

#' @export
data_dir.c <- function() {
  system.file("extdata", package = "ProMetIS")
}

#' @export
processed_dir.c <- function() {
  file.path(ProMetIS::data_dir.c(), "1_processed")
}

#' @export
post_processed_dir.c <- function() {
  file.path(ProMetIS::data_dir.c(), "2_post_processed")
}

#' @export
statistics_singleomics_dir.c <- function() {
  file.path(ProMetIS::data_dir.c(), "3_statistics_singleomics")
}

#' @export
statistics_integrative_dir.c <- function() {
  file.path(ProMetIS::data_dir.c(), "4_statistics_integrative")
}

#' @export
aggregated_dir.c <- function(gene.c) {
  stopifnot(gene.c %in% ProMetIS::genes.vc())
  file.path(ProMetIS::data_dir.c(), paste0("5", ifelse(gene.c == "LAT", "a", "b"), "_aggregated_", gene.c))
}

#' @export
supplementary_dir.c <- function(gene.c) {
  file.path(ProMetIS::data_dir.c(), "6_supplementary")
}