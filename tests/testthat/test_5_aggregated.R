testthat::context("Final aggregated datasets")


aggregated.ls <- sapply(ProMetIS::genes.vc(),
                        function(gene.c) {
                          phenomis::reading(ProMetIS::aggregated_dir.c(gene.c))
                        })


testthat::test_that("dimensions", {
  
  aggregated_dim.ls <- lapply(aggregated.ls,
                              function(gene.mset) {
                                t(sapply(names(gene.mset),
                                         function(set.c) dim(gene.mset[[set.c]])))
                              })

  # gene.c <- "MX2"
  # paste(paste(rownames(aggregated_dim.ls[[gene.c]]),
  #             apply(aggregated_dim.ls[[gene.c]], 1, function(x) paste(x, collapse = "|")),
  #             sep = "|"),
  #       collapse = "', '")
  test_dim.ls <- list(LAT = c("metabolomics_liver_c18hypersil_pos|5665|28",
                              "metabolomics_liver_hilic_neg|2866|28",
                              "metabolomics_plasma_c18acquity_neg|1584|28",
                              "metabolomics_plasma_c18acquity_pos|6104|28",
                              "metabolomics_plasma_c18hypersil_pos|4788|28",
                              "metabolomics_plasma_hilic_neg|3131|28",
                              "preclinical|203|28",
                              "proteomics_liver|2098|28",
                              "proteomics_plasma|419|24"),
                      MX2 = c("metabolomics_liver_c18hypersil_pos|5665|29",
                              "metabolomics_liver_hilic_neg|2866|29",
                              "metabolomics_plasma_c18acquity_neg|1584|29",
                              "metabolomics_plasma_c18acquity_pos|6104|29",
                              "metabolomics_plasma_c18hypersil_pos|4788|29",
                              "metabolomics_plasma_hilic_neg|3131|29",
                              "preclinical|234|29",
                              "proteomics_liver|2090|29",
                              "proteomics_plasma|422|25"))
  test_dim.ls <- lapply(test_dim.ls,
                        function(test_dim.vc) {
                          test_dim.mn <- sapply(test_dim.vc,
                                                function(x)
                                                  unlist(strsplit(x, "|", fixed = TRUE))[c(2, 3)])
                          mode(test_dim.mn) <- "integer"
                          dimnames(test_dim.mn) <- list(c("Features", "Samples"),
                                                        sapply(test_dim.vc,
                                                               function(x)
                                                                 unlist(strsplit(x, "|", fixed = TRUE))[1],
                                                               USE.NAMES = FALSE))
                          t(test_dim.mn)
                        })
  testthat::expect_identical(aggregated_dim.ls,
                             test_dim.ls)
  
})

testthat::test_that("significant", {
  
  aggregated_signif.mn <- sapply(names(aggregated.ls),
                                 function(gene.c) {
                                   gene.mset <- aggregated.ls[[gene.c]]
                                   sapply(names(gene.mset),
                                            function(set.c) {
                                              set_fda.df <- Biobase::fData(gene.mset[[set.c]])
                                              sum(set_fda.df[, "WT.KO_signif"], na.rm = TRUE)
                                              # if (gene.c == "LAT" && set.c == "proteomics_liver") {
                                              #   return(apply(set_fda.df[, c("limmaM_WT.LAT_signif",
                                              #                               "limmaF_WT.LAT_signif")],
                                              #                2,
                                              #                function(y) sum(y, na.rm = TRUE)))
                                              # } else {
                                              #   return(c(sum(set_fda.df[, paste0("WT.", gene.c, "_signif")],
                                              #                na.rm = TRUE),
                                              #            NA))
                                              # }
                                            })
                                 })
  
  # gene.c <- "MX2"
  # paste(paste(rownames(aggregated_signif.ls[[gene.c]]),
  #             apply(aggregated_signif.ls[[gene.c]], 1, function(x) paste(x, collapse = "|")),
  #             sep = "|"),
  #       collapse = "', '")
  test_signif.ls <- list(LAT = c("metabolomics_liver_c18hypersil_pos|1608",
                                 "metabolomics_liver_hilic_neg|826",
                                 "metabolomics_plasma_c18acquity_neg|2",
                                 "metabolomics_plasma_c18acquity_pos|13",
                                 "metabolomics_plasma_c18hypersil_pos|3",
                                 "metabolomics_plasma_hilic_neg|585",
                                 "preclinical|2",
                                 "proteomics_liver|257",
                                 "proteomics_plasma|7"),
                         MX2 = c("metabolomics_liver_c18hypersil_pos|27",
                                 "metabolomics_liver_hilic_neg|9",
                                 "metabolomics_plasma_c18acquity_neg|1",
                                 "metabolomics_plasma_c18acquity_pos|86",
                                 "metabolomics_plasma_c18hypersil_pos|81",
                                 "metabolomics_plasma_hilic_neg|48",
                                 "preclinical|2",
                                 "proteomics_liver|263",
                                 "proteomics_plasma|19"))
  test_signif.mn <- sapply(test_signif.ls,
                           function(test_signif.vc) {
                             test_signif.vn <- sapply(test_signif.vc, function(x) unlist(strsplit(x, "|", fixed = TRUE))[2])
                             mode(test_signif.vn) <- "integer"
                             names(test_signif.vn) <- sapply(test_signif.vc,
                                                             function(x) unlist(strsplit(x, "|", fixed = TRUE))[1],
                                                             USE.NAMES = FALSE)
                             test_signif.vn
                           })

  testthat::expect_identical(aggregated_signif.mn,
                             test_signif.mn)
  
})

testthat::test_that("proteomics_liver", {
  
  proteomics_liver.eset <- aggregated.ls[["LAT"]][["proteomics_liver"]]
  
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
                         2084)
  
  proteoliv_fWL.eset <- ProMetIS::subsetting(proteomics_liver.eset,
                                             set.c = "proteomics_liver",
                                             genes.vc = c("WT", "LAT"),
                                             sex.vc = "F")
  testthat::expect_equal(Biobase::dims(proteoliv_fWL.eset )["Features", ],
                         2094)
  
  proteoliv_mfW.eset <- ProMetIS::subsetting(proteomics_liver.eset,
                                             set.c = "proteomics_liver",
                                             genes.vc = "WT")
  testthat::expect_equal(Biobase::dims(proteoliv_mfW.eset )["Features", ],
                         2078)
  
  proteoliv_mfL.eset <- ProMetIS::subsetting(proteomics_liver.eset,
                                             set.c = "proteomics_liver",
                                             genes.vc = "LAT")
  testthat::expect_equal(Biobase::dims(proteoliv_mfL.eset )["Features", ],
                         2093)
  
})
