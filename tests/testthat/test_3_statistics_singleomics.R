testthat::context("Datasets after single omics statistics")

stat_intra.mset <- phenomis::reading(ProMetIS::statistics_singleomics_dir.c(), output.c = "set")

## dimensions ----

testthat::test_that("dimensions", {
  
  stat_intra_dim.mn <- t(sapply(names(stat_intra.mset),
                                function(set.c) dim(stat_intra.mset[[set.c]])))
  
  # paste(paste(rownames(stat_intra_dim.mn),
  #             apply(stat_intra_dim.mn, 1, function(x) paste(x, collapse = "|")),
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
  testthat::expect_identical(stat_intra_dim.mn,
                             test_dim.mn)
  
})

## significant ----

testthat::test_that("significant", {
  
  stat_intra_signif.mn <- t(sapply(names(stat_intra.mset),
                                   function(set.c) {
                                     set_fda.df <- Biobase::fData(stat_intra.mset[[set.c]])
                                     if (set.c == "proteomics_liver") {
                                       return(apply(set_fda.df[, c("limmaM_WT.LAT_signif",
                                                                   "limmaF_WT.LAT_signif",
                                                                   "limma2ways_WT.MX2_signif")],
                                                    2,
                                                    function(y) sum(y, na.rm = TRUE)))
                                     } else {
                                       return(c(sum(set_fda.df[, c("limma2ways_WT.LAT_signif")], na.rm = TRUE),
                                                NA,
                                                sum(set_fda.df[, c("limma2ways_WT.MX2_signif")], na.rm = TRUE)))
                                     }
                                   }))
  colnames(stat_intra_signif.mn) <- c("LAT(_M)", "LAT_F", "MX2")
  
  # paste(paste(rownames(stat_intra_signif.mn),
  #             apply(stat_intra_signif.mn, 1, function(x) paste(x, collapse = "|")),
  #             sep = "|"),
  #       collapse = "', '")
  test_signif.vc <- c("metabolomics_liver_c18hypersil_pos|1608|NA|27",
                      "metabolomics_liver_hilic_neg|826|NA|9",
                      "metabolomics_plasma_c18acquity_neg|2|NA|1",
                      "metabolomics_plasma_c18acquity_pos|13|NA|86",
                      "metabolomics_plasma_c18hypersil_pos|3|NA|81",
                      "metabolomics_plasma_hilic_neg|585|NA|48",
                      "preclinical|2|NA|2",
                      "proteomics_liver|1|258|263",
                      "proteomics_plasma|7|NA|19")
  test_signif.mn <- sapply(test_signif.vc, function(x) unlist(strsplit(x, "|", fixed = TRUE))[2:4])
  suppressWarnings(mode(test_signif.mn) <- "integer")
  dimnames(test_signif.mn) <- list(c("LAT(_M)", "LAT_F", "MX2"),
                                   sapply(test_signif.vc, function(x) unlist(strsplit(x, "|", fixed = TRUE))[1], USE.NAMES = FALSE))
  test_signif.mn <- t(test_signif.mn)
  testthat::expect_identical(stat_intra_signif.mn,
                             test_signif.mn)
  
})

## proteomics_liver ----

testthat::test_that("proteomics_liver", {
  
  proteomics_liver.eset <- stat_intra.mset[["proteomics_liver"]]
  
  proteoliv_WL.eset <- ProMetIS::subsetting(proteomics_liver.eset,
                                          set.c = "proteomics_liver",
                                          genes.vc = c("WT", "LAT"))
  testthat::expect_equal(Biobase::dims(proteoliv_WL.eset )["Features", ],
                         2098)
  
  proteoliv_mWL.eset <- ProMetIS::subsetting(proteomics_liver.eset,
                                           set.c = "proteomics_liver",
                                           genes.vc = c("WT", "LAT"),
                                           sex.vc = "M")
  testthat::expect_equal(Biobase::dims(proteoliv_mWL.eset )["Features", ],
                         2120)
  
  proteoliv_fWL.eset <- ProMetIS::subsetting(proteomics_liver.eset,
                                           set.c = "proteomics_liver",
                                           genes.vc = c("WT", "LAT"),
                                           sex.vc = "F")
  testthat::expect_equal(Biobase::dims(proteoliv_fWL.eset )["Features", ],
                         2139)
  
  proteoliv_mfW.eset <- ProMetIS::subsetting(proteomics_liver.eset,
                                           set.c = "proteomics_liver",
                                           genes.vc = "WT")
  testthat::expect_equal(Biobase::dims(proteoliv_mfW.eset )["Features", ],
                         2126)
  
  proteoliv_mfL.eset <- ProMetIS::subsetting(proteomics_liver.eset,
                                           set.c = "proteomics_liver",
                                           genes.vc = "LAT")
  testthat::expect_equal(Biobase::dims(proteoliv_mfL.eset )["Features", ],
                         2145)
  
})

