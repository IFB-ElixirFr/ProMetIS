## Metabolomics ----

# Postprocessing of metabolomics datasets
#  Imputation
#  RT filtering
#  Blank filtering
#  Pool dilution correlation
#  Signal drift correction
#  poolCV <= 0.3
#  poolCV_over_sampleCV <= 1
#  Chemical redundancy (Monnerie et al., 1999)
.metabo_postprocessing <- function(metabo.mset,
                                   drift_correct.c = c("none", "pool", "sample", "prometis")[4],
                                   span.n = 1,
                                   .technical_validation.l = FALSE) { # removing poolCV filters and keeping the pools for the "Technical Validation" section
  
  for (set.c in names(metabo.mset)) {
    
    print(set.c)
    
    eset <- metabo.mset[[set.c]]
    
    # # Imputation (MTH Paris)
    # 
    # if (grepl("c18hyper", set.c)) {
    #   
    #   exprs.mn <- Biobase::exprs(eset)
    #   
    #   exprs.mn[is.na(exprs.mn)] <- 1
    #   
    #   Biobase::exprs(eset) <- exprs.mn
    #   
    # }
    
    # RT filtering (MTH Clermont)
    
    ## Converting RT in seconds
    rt.vn <- Biobase::fData(eset)[, "rt"]
    if (max(rt.vn) < 60)
      Biobase::fData(eset)[, "rt"] <- rt.vn * 60
    
    if (grepl("c18acqui", set.c))
      eset <- eset[which(Biobase::fData(eset)[, "rt"] >= 0.4 * 60 &
                           Biobase::fData(eset)[, "rt"] <= 22 * 60), ]
    
    # Blank filtering
    
    eset <- phenomis::inspecting(eset, figure.c = "none", report.c = "none")
    
    select_sample_specific.vl <- Biobase::fData(eset)[, "blankMean_over_sampleMean"] <= 0.33
    if (.technical_validation.l && grepl("(hyper|hilic)", set.c)) {
      select_sample_specific.vl <- select_sample_specific.vl |
        Biobase::fData(eset)[, "standard"] != ""
    }

    eset <- eset[which(select_sample_specific.vl), ]

    # Pool dilution (MTH Paris)
    
    if (grepl("(hyper|hilic)", set.c)) {
      
      select_poolDil_correl.vl <- Biobase::fData(eset)[, "poolDil_cor"] >= 0.7
      if (.technical_validation.l) {
        select_poolDil_correl.vl <- select_poolDil_correl.vl |
          Biobase::fData(eset)[, "standard"] != ""
      }
      
      eset <- eset[which(select_poolDil_correl.vl), ]
      
      ## Setting pool1 to pool for subsequent use in the 'pool CV < 30%' filter
      
      Biobase::pData(eset)[, "sampleType"] <- gsub("pool1", "pool",
                                                   Biobase::pData(eset)[, "sampleType"])
      
    }
    
    # Discarding blanks and pool dilutions
    
    eset <- eset[, which(Biobase::pData(eset)[, "sampleType"] %in% c("sample", "pool"))]
    
    # Signal drift correction
    
    ## ordering eset according to injection order
    
    eset <- eset[, order(Biobase::pData(eset)[, "injectionOrder"])]
    
    if (grepl("acqui", set.c)) { ## resetting injection order min to 1 (MTH Clermont)
      pdata.df <- Biobase::pData(eset)
      pdata.df[, "injectionOrder"] <- pdata.df[, "injectionOrder"] - min(pdata.df[, "injectionOrder"]) + 1
      Biobase::pData(eset) <- pdata.df
    }
    
    ## keeping only the last 'pool' before the first 'sample' (to limit extrapolation)
    
    pdata.df <- Biobase::pData(eset)
    first_sample.i <- which(pdata.df[, "sampleType"] == "sample")[1]
    last_pool_before_first_sample.i <- first_sample.i - 1 
    stopifnot(pdata.df[last_pool_before_first_sample.i, "sampleType"] == "pool")
    
    eset <- eset[, seq(last_pool_before_first_sample.i, nrow(pdata.df))]
    
    ## signal drift correction
    
    if (drift_correct.c == "prometis") {
      if (grepl("plasma.+(hyper|hilic)", set.c)) {
        eset <- phenomis::correcting(eset,
                                     reference.c = "pool",
                                     title.c = gsub("metabolomics_", "", set.c),
                                     span.n = span.n,
                                     figure.c = ifelse(.technical_validation.l, "interactive", "none"))
      } else if (grepl("acqui", set.c))
        eset <- phenomis::correcting(eset,
                                     reference.c = "sample",
                                     title.c = gsub("metabolomics_", "", set.c),
                                     span.n = span.n,
                                     figure.c = ifelse(.technical_validation.l, "interactive", "none"))
    } else if (drift_correct.c != "none")
      eset <- phenomis::correcting(eset,
                                   reference.c = drift_correct.c,
                                   title.c = gsub("metabolomics_", "", set.c),
                                   span.n = span.n,
                                   figure.c = ifelse(.technical_validation.l, "interactive", "none"))
    
    # NAs and variances
    
    eset <- phenomis::filtering(eset, class.c = ifelse(grepl("(hyper|hilic)", set.c), "KO", "Type"), max_na_prop.n = 0.2)
    
    if (!.technical_validation.l) {
      
      # pool_CV <= 0.3
      
      eset <- phenomis::inspecting(eset, report.c = "none", figure.c = "none")
      
      eset <- eset[which(Biobase::fData(eset)[, "pool_CV"] <= 0.3), ]
      
      # poolCV_over_sampleCV <= 1
      
      eset <- eset[which(Biobase::fData(eset)[, "poolCV_over_sampleCV"] <= 1), ]
      
      # Discarding pools
      
      eset <- eset[, which(Biobase::pData(eset)[, "sampleType"] != "pool")]
      
    }
    
    # Reducing chemical redundancy (Monnerie et al., 2019)
    
    eset <- phenomis::reducing(eset)
    
    select_non_redund.vl <- Biobase::fData(eset)[, "redund_is"] < 1
    if (.technical_validation.l && grepl("(hyper|hilic)", set.c)) {
      select_non_redund.vl <- select_non_redund.vl |
        Biobase::fData(eset)[, "standard"] != ""
    }
    
    eset <- eset[which(select_non_redund.vl), ]
    
    # Updating the eset object
    
    stopifnot(methods::validObject(eset))
    
    metabo.mset <- MultiDataSet::add_eset(metabo.mset,
                                          eset,
                                          dataset.type = set.c,
                                          GRanges = NA,
                                          overwrite = TRUE,
                                          warnings = FALSE)
    
  }
  
  return(metabo.mset)
  
}


.format_metabonames <- function(metabo.mset,
                                mice_id.df) {
  # mice_id.df, mice_id.vc, mice_num.vc) {
  
  for (set.c in ProMetIS::metabo_sets.vc()) {
    
    eset <- metabo.mset[[set.c]]
    
    # sample names formatting and ordering by mouse number
    
    if (grepl("(hyper|hilic)", set.c)) {
      samp_num.vc <- substr(Biobase::sampleNames(eset), 11, 13)
    } else if (grepl("acqui", set.c)) {
      samp_num.vc <- substr(Biobase::sampleNames(eset), 10, 12)
    } else
      stop("Unknown metabolomics dataset name.")
    
    stopifnot(identical(sort(samp_num.vc), sort(as.character(mice_id.df[, "mouse_id"]))))
    
    if (grepl("acqui", set.c)) {
      eset <- eset[, order(as.numeric(samp_num.vc))]
      samp_num.vc <- substr(Biobase::sampleNames(eset), 10, 12)
    }
    
    stopifnot(identical(sort(samp_num.vc), sort(as.character(mice_id.df[, "mouse_id"]))))
    
    samp_ord.vi <- order(samp_num.vc)
    eset <- eset[, samp_ord.vi]
    samp_num.vc <- samp_num.vc[samp_ord.vi]
    
    stopifnot(identical(samp_num.vc, as.character(mice_id.df[, "mouse_id"])))
    
    Biobase::pData(eset) <- cbind.data.frame(mice_id.df,
                                             initial_name = Biobase::sampleNames(eset),
                                             Biobase::pData(eset))
    
    Biobase::sampleNames(eset) <- rownames(mice_id.df)
    
    # variable metadata: adding the name of the chromatographic column

    fdata.df <- Biobase::fData(eset)

    fdata.df[, "chromato"] <- rep(unlist(strsplit(set.c, split = "_"))[3],
                                  nrow(fdata.df))
    
    if (grepl("acqui", set.c)) {
      for (col.c in c("name", "chebi_id", "annot_confidence")) {
        if (col.c %in% colnames(fdata.df)) {
          fdata.df[, col.c] <- as.character(fdata.df[, col.c])
          fdata.df[is.na(fdata.df[, col.c]), col.c] <- ""
        }
      }
    }

    Biobase::fData(eset) <- fdata.df

    stopifnot(methods::validObject(eset))
    
    metabo.mset <- MultiDataSet::add_eset(metabo.mset,
                                          eset,
                                          dataset.type = set.c,
                                          GRanges = NA,
                                          overwrite = TRUE,
                                          warnings = FALSE)
    
  }
  
  metabo.mset
  
}