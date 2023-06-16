testthat::context("Datasets after post-processing")

post_proc.mset <- phenomis::reading(ProMetIS::post_processed_dir.c(), output.c = "set")


## values ----

testthat::test_that("values", {
  
  testthat::expect_equal(Biobase::exprs(post_proc.mset[["proteomics_liver"]])["Q8C196_Carbamoyl-phosphate synt.",
                                                                              "W627m"],
                         32.1570347762747, tol = 1e-13)
  
  testthat::expect_equal(Biobase::exprs(post_proc.mset[["proteomics_plasma"]])["P01027_Complement C3",
                                                                               "W617f"],
                         32.4223086299966, tol = 1e-13)
  
})


## dimensions ----

testthat::test_that("dimensions", {
  
  post_proc_dim.mn <- t(sapply(names(post_proc.mset),
                               function(set.c) dim(post_proc.mset[[set.c]])))
  
  # paste(paste(rownames(post_proc_dim.mn),
  #             apply(post_proc_dim.mn, 1, function(x) paste(x, collapse = "|")),
  #             sep = "|"),
  #       collapse = "', '")
  test_dim.vc <- c("metabolomics_liver_c18hypersil_pos|5665|42",
                   "metabolomics_liver_hilic_neg|2866|42",
                   "metabolomics_plasma_c18acquity_neg|1584|42",
                   "metabolomics_plasma_c18acquity_pos|6104|42",
                   "metabolomics_plasma_c18hypersil_pos|4788|42",
                   "metabolomics_plasma_hilic_neg|3131|42",
                   "preclinical|200|42",
                   "proteomics_liver|2187|42",
                   "proteomics_plasma|446|36")
  test_dim.mn <- sapply(test_dim.vc, function(x) unlist(strsplit(x, "|", fixed = TRUE))[c(2, 3)])
  mode(test_dim.mn) <- "integer"
  dimnames(test_dim.mn) <- list(c("Features", "Samples"),
                                sapply(test_dim.vc, function(x) unlist(strsplit(x, "|", fixed = TRUE))[1], USE.NAMES = FALSE))
  test_dim.mn <- t(test_dim.mn)
  testthat::expect_identical(post_proc_dim.mn,
                             test_dim.mn)
  
})


## impute_info ----

testthat::test_that("impute_info", {
  
  set.c <- "proteomics_plasma"
  
  # SummarizedExperiment
  proteo.mae <- phenomis::reading(ProMetIS::post_processed_dir.c(),
                                  subsets.vc = set.c)
  proteo.se <- proteo.mae[[set.c]]
  imputed.mi <- ProMetIS:::imputation_info(proteo.se, set.c = set.c)
  
  genotype.fc <- factor(substr(colnames(imputed.mi), 1, 1), levels = c("L", "X", "W"))
  geno_impute.mn <- t(apply(imputed.mi, 1, function(feat.vn) tapply(feat.vn, genotype.fc, sum)))
  colnames(geno_impute.mn) <- c("LAT", "MX2", "WT")
  
  testthat::expect_equivalent(colSums(geno_impute.mn),
                              c(225, 221, 351))
  
  # ExpressionSet
  
  proteo.mset <- phenomis::reading(ProMetIS::post_processed_dir.c(),
                                   subsets.vc = set.c,
                                   output.c = "set")
  proteo.eset <- proteo.mset[[set.c]]
  imputed.mi <- ProMetIS:::imputation_info(proteo.eset, set.c = set.c)
  # computing the number of imputed values per feature and per genotype
  genotype.fc <- factor(substr(colnames(imputed.mi), 1, 1), levels = c("L", "X", "W"))
  geno_impute.mn <- t(apply(imputed.mi, 1, function(feat.vn) tapply(feat.vn, genotype.fc, sum)))
  colnames(geno_impute.mn) <- c("LAT", "MX2", "WT")

  testthat::expect_equivalent(colSums(geno_impute.mn),
                              c(225, 221, 351))
  
})


## filter_overimputed ----

testthat::test_that("filter_overimputed", {
  
  set.c <- "proteomics_plasma"
  
  # SummarizedExperiment
  
  proteo.mae <- phenomis::reading(ProMetIS::post_processed_dir.c(),
                                  subsets.vc = set.c)
  proteo.se <- proteo.mae[[set.c]]
  sample_meta.DF <- SummarizedExperiment::colData(proteo.mae)
  SummarizedExperiment::colData(proteo.se)[, "gene"] <- sample_meta.DF[colnames(proteo.se), "gene"]
  SummarizedExperiment::colData(proteo.se)[, "sex"] <- sample_meta.DF[colnames(proteo.se), "sex"]
  feat_select.vl <- ProMetIS:::filter_overimputed(proteo.se, set.c = set.c, genes.vc = "LAT")
  
  testthat::expect_equal(sum(feat_select.vl),
                         416, tol = 1e-10)

  # ExpressionSet
  
  proteo.mset <- phenomis::reading(ProMetIS::post_processed_dir.c(),
                                   subsets.vc = set.c,
                                   output.c = "set")
  proteo.eset <- proteo.mset[[set.c]]
  feat_select.vl <- ProMetIS:::filter_overimputed(proteo.eset, set.c = set.c, genes.vc = "LAT")

  testthat::expect_equal(sum(feat_select.vl),
                         416, tol = 1e-10)
})

## subsetting ----

testthat::test_that("subsetting", {
  
  # SummarizedExperiment
  
  set.c <- "proteomics_liver"
  proteo.mae <- phenomis::reading(ProMetIS::post_processed_dir.c(),
                                  subsets.vc = set.c,
                                  report.c = "none")
  proteo.se <- proteo.mae[[set.c]]
  ## adding (common) sample metadata to the summarized experiment
  sample_meta.DF <- SummarizedExperiment::colData(proteo.mae)
  SummarizedExperiment::colData(proteo.se)[, "gene"] <- sample_meta.DF[colnames(proteo.se), "gene"]
  SummarizedExperiment::colData(proteo.se)[, "sex"] <- sample_meta.DF[colnames(proteo.se), "sex"]
  proteo_lat.se <- ProMetIS::subsetting(proteo.se,
                                        set.c = set.c,
                                        genes.vc = c("WT", "LAT"))
  
  testthat::expect_equal(dim(proteo.se), c(2187, 42))
  testthat::expect_equal(dim(proteo_lat.se), c(2098, 28))
  
  
})
