testthat::context("Datasets after post-processing")

post_proc.mset <- phenomis::reading(ProMetIS::post_processed_dir.c())

testthat::test_that("values", {
  
  testthat::expect_equal(Biobase::exprs(post_proc.mset[["proteomics_liver"]])["Q8C196_Carbamoyl-phosphate synt.",
                                                                              "W627m"],
                         32.1570347762747, tol = 1e-13)
  
  testthat::expect_equal(Biobase::exprs(post_proc.mset[["proteomics_plasma"]])["P01027_Complement C3",
                                                                               "W617f"],
                         32.4223086299966, tol = 1e-13)
  
})

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
                   "preclinical|236|42",
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

