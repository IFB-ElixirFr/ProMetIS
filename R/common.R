## All ----

#### subsetting (MultiDataSet) ####

#' @rdname subsetting
#' @export
setMethod("subsetting", signature(x = "MultiDataSet"),
          function(x,
                   set.c = NULL,
                   genes.vc = "all",
                   sex.vc = "all",
                   tissues.vc = "all",
                   common_samples.l = FALSE,
                   na_thresh.n = 0.2,
                   var_thresh.n = .Machine$double.eps,
                   imputed_thresh.n = 0.2) {
            
            if (length(genes.vc) == 1 && genes.vc == "all")
              genes.vc <- c("WT", ProMetIS::genes.vc())
            stopifnot(all(genes.vc %in% c("WT", ProMetIS::genes.vc())))
            
            if (length(sex.vc) == 1 && sex.vc == "all")
              sex.vc <- ProMetIS::sex.vc()
            stopifnot(all(sex.vc %in% ProMetIS::sex.vc()))
            
            stopifnot(length(tissues.vc) == 1)
            
            if (tissues.vc != "all") {
              
              stopifnot(tissues.vc %in% ProMetIS::tissues.vc())
              
              x <- x[, names(x)[names(x) %in% c("preclinical",
                                                grep(tissues.vc,
                                                     ProMetIS::sets.vc(),
                                                     value = TRUE))]]
              
            }
            
            if (common_samples.l) {
              x <- MultiDataSet::commonSamples(x)
              ref.eset <- x[[names(x)[1]]]
            } else if ("preclinical" %in% names(x)) {
              ref.eset <- x[["preclinical"]]
            } else {
              sample_names.ls <- Biobase::sampleNames(x)
              all_samples.vc <- Reduce("union", sample_names.ls)
              ref_eset.i <- which.max(unlist(lapply(Biobase::sampleNames(x),
                                                    function(sample_names.vc)
                                                      sum(sample_names.vc %in% all_samples.vc))))[1]
              ref.eset <- x[[names(x)[[ref_eset.i]]]]
            }
            
            samples.vc <- Biobase::sampleNames(ref.eset)
            samples_sel.vl <- Biobase::pData(ref.eset)[, "gene"] %in% genes.vc &
              Biobase::pData(ref.eset)[, "sex"] %in% sex.vc
            
            sub.mset <- x[samples.vc[samples_sel.vl], ]
            
            sub_sets.vc <- names(sub.mset)
            
            for (set.c in sub_sets.vc) {
              
              eset <- sub.mset[[set.c]]
              
              eset <- subsetting(x = eset,
                                 set.c = set.c,
                                 genes.vc = genes.vc,
                                 sex.vc = sex.vc,
                                 tissues.vc = NULL,
                                 common_samples.l = NULL,
                                 na_thresh.n = na_thresh.n,
                                 var_thresh.n = var_thresh.n,
                                 imputed_thresh.n = imputed_thresh.n)
              
              sub.mset <- MultiDataSet::add_eset(sub.mset,
                                                 eset,
                                                 dataset.type = set.c,
                                                 GRanges = NA,
                                                 overwrite = TRUE,
                                                 warnings = FALSE)
              
            }
            
            sub.mset <- sub.mset[, sub_sets.vc]
            
            invisible(sub.mset)
            
          })


#### subsetting (ExpressionSet) ####

#' @rdname subsetting
#' @export
setMethod("subsetting", signature(x = "ExpressionSet"),
          function(x,
                   set.c,
                   genes.vc = "all",
                   sex.vc = "all",
                   tissues.vc = NULL,
                   common_samples.l = NULL,
                   na_thresh.n = 0.2,
                   var_thresh.n = .Machine$double.eps,
                   imputed_thresh.n = 0.2) {
            
            if (length(genes.vc) == 1 && genes.vc == "all")
              genes.vc <- ProMetIS::wtgenes.vc()
            stopifnot(all(genes.vc %in% ProMetIS::wtgenes.vc()))
            
            if (length(sex.vc) == 1 && sex.vc == "all")
              sex.vc <- ProMetIS::sex.vc()
            stopifnot(all(sex.vc %in% ProMetIS::sex.vc()))
            
            samples.vc <- Biobase::sampleNames(x)
            
            samples_sel.vl <- Biobase::pData(x)[, "gene"] %in% genes.vc &
              Biobase::pData(x)[, "sex"] %in% sex.vc
            
            x <- x[, samples.vc[samples_sel.vl]]
            
            filter.vi <- c(nas_zerovar = NA_integer_,
                           overimputed = NA_integer_)
            
            # NAs <= 20% and variance >= 1e-5
            if (length(genes.vc) > 1) {
              class.c <- "gene"
            } else if (length(sex.vc) > 1) {
              class.c <- "sex"
            } else
              stop("Only a single gene and sex selected.")
            filtered.eset <- phenomis::filtering(x,
                                                 class.c = class.c,
                                                 max_na_prop.n = na_thresh.n,
                                                 min_variance.n = var_thresh.n)
            na_zerovar_sel.vl <- Biobase::featureNames(x) %in% Biobase::featureNames(filtered.eset)
            
            # na_zerovar_sel.vl <- ProMetIS:::.filter_na_zerovar(t(Biobase::exprs(x)),
            #                                                    na_thresh.n = na_thresh.n,
            #                                                    var_thresh.n = var_thresh.n)
            
            filter.vi["nas_zerovar"] <- sum(!na_zerovar_sel.vl)
            
            # proteomics: observations >= 80% in at least one condition
            if (grepl("proteomics", set.c)) {
              overimputed_sel.vl <- ProMetIS:::.filter_overimputed(eset = x,
                                                                   set.c = set.c,
                                                                   genes.vc = genes.vc,
                                                                   sex.vc = sex.vc,
                                                                   imputed_thresh.n = imputed_thresh.n)
            } else
              overimputed_sel.vl <- rep(TRUE, dim(x)["Features"])
            
            filter.vi["overimputed"] <- sum(!overimputed_sel.vl)
            
            # intersection of both conditions
            feat_sel.vl <- na_zerovar_sel.vl & overimputed_sel.vl
            
            stopifnot(length(feat_sel.vl) == dim(x)["Features"] &&
                        !any(is.na(feat_sel.vl)))
            
            x <- x[feat_sel.vl, ]
            
            # if (sum(!feat_sel.vl))
              message("Nb of discard. feat. in '", set.c, "': ",
                      paste(paste0(names(filter.vi), ": ", filter.vi), collapse = ", "))
            
            invisible(x)
            
          })


# .filter_na_zerovar <- function(input.mn,
#                                class.c = "",
#                                na_thresh.n = 0.2,
#                                var_thresh.n = .Machine$double.eps) {
#   
#   # removing variables with > 20% NA (including 100% NA in females)
#   feat_na.vn <- apply(input.mn, 2, function(feat.vn)
#     sum(is.na(feat.vn)) / length(feat.vn))
#   feat_notna.vl <- feat_na.vn <= na_thresh.n
#   # sum(feat_notna.vl)
#   
#   # removing variables with variance < 1e-5
#   feat_var.vn <- apply(input.mn, 2, function(feat.vn)
#     var(feat.vn, na.rm = TRUE))
#   feat_notzerovar.vl <- !is.na(feat_var.vn) &
#     (feat_var.vn >= var_thresh.n)
#   # sum(feat_notzerovar.vl)
#   
#   feat_sel.vl <- feat_notna.vl & feat_notzerovar.vl
#   
#   stopifnot(length(feat_sel.vl) == ncol(input.mn) &&
#               !any(is.na(feat_sel.vl)))
#   
#   feat_sel.vl
#   
# }

.imputation_info <- function(eset,
                             set.c) {
  
  prot_pda.df <- Biobase::pData(eset)
  prot_fda.df <- Biobase::fData(eset)
  
  load(system.file("extdata/2_post_processed/metadata_supp.rdata", package = "ProMetIS"))
  
  supp_pda.df <- metadata_supp.ls[[set.c]][["pdata"]]
  supp_fda.df <- metadata_supp.ls[[set.c]][["fdata"]]
  
  stopifnot(all(rownames(prot_pda.df) %in% rownames(supp_pda.df)))
  stopifnot(all(rownames(prot_fda.df) %in% rownames(supp_fda.df)))
  
  prot_pda.df <- cbind.data.frame(prot_pda.df, supp_pda.df[rownames(prot_pda.df), , drop = FALSE])
  prot_fda.df <- cbind.data.frame(prot_fda.df, supp_fda.df[rownames(prot_fda.df), , drop = FALSE])
  
  
  ## checking that the sample names are ordered by increasing ID
  prot_samp.vi <- as.integer(substr(rownames(prot_pda.df), 2, 4))
  stopifnot(identical(prot_samp.vi, sort(prot_samp.vi)))
  
  ## getting imputation info
  value_origin.vl <- vapply(colnames(prot_fda.df), function(colname.c) {
    colname_split.vc <- unlist(strsplit(colname.c, split = "_"))
    grepl("^OriginOfValueabundance", colname.c) &
      colname_split.vc[length(colname_split.vc)] %in% rownames(prot_pda.df)
  }, FUN.VALUE = logical(1))
  value_origin.df <- prot_fda.df[, value_origin.vl]
  colnames(value_origin.df) <- gsub("_run90methode30K",
                                    "",
                                    gsub("_mgf", "",
                                         gsub("OriginOfValueabundance_", "",
                                              colnames(value_origin.df))))
  
  ## re-ordering imputation info to match sample names
  value_origin_samp.vc <- vapply(colnames(value_origin.df), function(colname.c) {
    colname_split.vc <- unlist(strsplit(colname.c, split = "_"))
    colname_split.vc[length(colname_split.vc)]}, FUN.VALUE = character(1))
  temp <- value_origin_samp.vc
  value_origin_samp.vc <- names(value_origin_samp.vc)
  names(value_origin_samp.vc) <- temp
  value_origin.df <- value_origin.df[, value_origin_samp.vc[rownames(prot_pda.df)]]
  
  stopifnot(identical(names(value_origin_samp.vc[rownames(prot_pda.df)]), Biobase::sampleNames(eset)))
  colnames(value_origin.df) <- Biobase::sampleNames(eset)
  
  imputed.mi <- apply(value_origin.df, 2, DAPAR_is.MV)
  mode(imputed.mi) <- "integer"
  
  stopifnot(!any(is.na(c(imputed.mi))))
  
  return(imputed.mi)
  
}

.filter_overimputed <- function(eset,
                                set.c,
                                genes.vc,
                                sex.vc,
                                imputed_thresh.n) {

  imputed.mi <- .imputation_info(eset = eset, set.c = set.c)
  
  stopifnot(identical(Biobase::sampleNames(eset), colnames(imputed.mi)))
  
  if (length(genes.vc) >= 2) {
    # general case: no restriction about sex
    # or e.g. LAT vs WT on males (or females) only
    
    factor.fc <- factor(Biobase::pData(eset)[, "gene"],
                        levels = ProMetIS::wtgenes.vc()[ProMetIS::wtgenes.vc() %in% genes.vc])
    
  } else if (length(genes.vc) == 1 && length(sex.vc) == 2) {
    # males vs females on LAT (or WT) only
    
    factor.fc <- factor(Biobase::pData(eset)[, "sex"],
                        levels = ProMetIS::sex.vc())
    
  } else
    stop("The corresponding imputation metric for this combination of genotype(s) and gender(s) is not currently available.",
         call. = FALSE)
  
  imputed_factor.mi <- t(apply(imputed.mi, 1, function(var.vn) {
    tapply(var.vn, factor.fc, sum)
  }))
  imputed_prop.mn <- sweep(imputed_factor.mi, 2, table(factor.fc), "/")
  
  imputed.ml <- imputed_prop.mn >= imputed_thresh.n
  imputed.vn <- rowSums(imputed.ml, na.rm = TRUE)
  
  feat_sel.vl <- imputed.vn <= 1
  
  stopifnot(length(feat_sel.vl) == dim(eset)["Features"] &&
              !any(is.na(feat_sel.vl)))
  
  feat_sel.vl
  
}


metadata_select <- function(mset,
                            step.c) { # e.g. step.c = "2_post_processed"
  
  if(is.null(names(mset))) { # preclinical expression set
    preclinical.eset <- mset
    mset <- MultiDataSet::createMultiDataSet()
    mset <- MultiDataSet::add_eset(mset,
                                   preclinical.eset,
                                   dataset.type = "preclinical",
                                   GRanges = NA,
                                   overwrite = TRUE,
                                   warnings = FALSE)
  }
  
  mset_names.vc <- names(mset)
  
  metadata_supp_file.c <- paste0("../inst/extdata/", step.c, "/metadata_supp.rdata")
  
  if(file.exists(metadata_supp_file.c)) {
    load(metadata_supp_file.c)
  } else {
    metadata_supp.ls <- vector(mode = "list", length = length(ProMetIS::sets.vc()))
    names(metadata_supp.ls) <- ProMetIS::sets.vc()
  }
  
  for (set.c in mset_names.vc) {
    
    eset <- mset[[set.c]]
    
    # sample metadata
    
    pdata.df <- Biobase::pData(eset)
    
    samplemeta.vc <- .sample_metadata_select(set.c)
    
    samplemeta.vc <- samplemeta.vc[samplemeta.vc %in% colnames(pdata.df)]
    
    samplemeta_supp.vc <- setdiff(colnames(pdata.df), samplemeta.vc)
    
    if(length(samplemeta_supp.vc)) {
      pdata_supp.df <- pdata.df[, samplemeta_supp.vc, drop = FALSE]
    } else
      pdata_supp.df <- data.frame()
    
    Biobase::pData(eset) <- pdata.df[, samplemeta.vc]
    
    # variable metadata
    
    fdata.df <- Biobase::fData(eset)
    
    variablemeta.vc <- .variable_metadata_select(fdata.df = fdata.df, set.c = set.c)
    
    variablemeta.vc <- variablemeta.vc[variablemeta.vc %in% colnames(fdata.df)]
    
    variablemeta_supp.vc <- setdiff(colnames(fdata.df), variablemeta.vc)
    
    if(length(variablemeta_supp.vc)) {
      fdata_supp.df <- fdata.df[, variablemeta_supp.vc, drop = FALSE]
    } else
      fdata_supp.df <- data.frame()
    
    if (grepl("metabolomics", set.c)) {
      annot_level.vc <- fdata.df[, "annot_level"]
      annot_level.vc <- as.character(annot_level.vc)
      annot_level_na.vi <- which(is.na(annot_level.vc) | annot_level.vc == "NA")
      annot_level.vc[annot_level_na.vi] <- ""
      fdata.df[, "annot_level"] <- annot_level.vc
    }
    
    Biobase::fData(eset) <- fdata.df[, variablemeta.vc]
    
    stopifnot(methods::validObject(eset))
    
    mset <- MultiDataSet::add_eset(mset,
                                   eset,
                                   dataset.type = set.c,
                                   GRanges = NA,
                                   overwrite = TRUE,
                                   warnings = FALSE)
    
    metadata_supp.ls[[set.c]] <- list(pdata = pdata_supp.df,
                                      fdata = fdata_supp.df)
 
    
  }
  
  mset <- mset[, mset_names.vc]
  
  save(metadata_supp.ls, file = metadata_supp_file.c)
  
  message("Supplementary metadata written in:\n", metadata_supp_file.c)
  
  return(invisible(mset))
  
}

.sample_metadata_select <- function(set.c) {
  
  first.vc <- c("gene",
                "mouse_id",
                "sex")
  
  if(set.c == "preclinical")
    first.vc <- c(first.vc,
                  c("gene_name",
                    "mgi_id",
                    "variant",
                    "impc_id",
                    "impc_genotype",
                    "impc_phenotype"))
  
  return(first.vc)
  
}

.variable_metadata_select <- function(fdata.df, set.c) {

  # post-processing
  
  if (set.c == "preclinical")
    varmeta.vc <- c("category",
                    "measurement",
                    "category_full",
                    "transformation")
  
  if (grepl("metabolomics", set.c))
    varmeta.vc <- c("chromato",
                    "name",
                    "chebi_id",
                    "annot_level",
                    "annot_confidence")
  
  if (grepl("proteomics", set.c))
    varmeta.vc <- c("accession",
                     "description",
                     "uniprot_id")

  
  # hypothesis testing
  
  limma_col.vc <- grep("limma", colnames(fdata.df), value = TRUE)
  if (length(limma_col.vc))
    varmeta.vc <- c(varmeta.vc,
                    limma_col.vc)
  
  # VIP
  
  vip_col.vc <- grep("OPLSDA_VIP-pred", colnames(fdata.df), value = TRUE)
  if (length(vip_col.vc))
    varmeta.vc <- c(varmeta.vc,
                    vip_col.vc)
  
  # Feature selection
  
  biosign_col.vc <- grep("biosign_", colnames(fdata.df), value = TRUE)
  if (length(biosign_col.vc))
    varmeta.vc <- c(varmeta.vc,
                    biosign_col.vc)
  
  # mixOmics
  
  mixomics_col.vc <- grep("mixomics_", colnames(fdata.df), value = TRUE)
  if (length(mixomics_col.vc))
    varmeta.vc <- c(varmeta.vc,
                    mixomics_col.vc)
  
  # complementary metabolomics annotation
  
  if (grepl("metabolomics", set.c)) {

    if (grepl("(hyper|hilic)", set.c))
      varmeta.vc <- c(varmeta.vc,
                      "kegg_id",
                      "kegg_pathway_family",
                      "kegg_pathways",
                      "kegg_subpathways",
                      "hmdb_id",
                      "pubchem_id",
                      "formula",
                      "monoisotopic_mass",
                      "inchikey",
                      "inchi")
    
    varmeta.vc <- c(varmeta.vc,
                    "MT",
                    "mz",
                    "rt",
                    "isotopes",
                    "adduct",
                    "pcgroup",
                    "redund_group",
                    "redund_iso_add_frag")
    
  }
    
  return(varmeta.vc)
  
}


#' @export
abbrev_mset <- function(mset, full_to_short.l = TRUE) {
  
  mset_names.vc <- names(mset)
  
  abbrev.mset <- MultiDataSet::createMultiDataSet()
  
  for (set.c in mset_names.vc) {
    
    eset <- mset[[set.c]]
    
    abbrev.c <- ProMetIS::sets_abbrev.vc(full_to_short.l = full_to_short.l)[set.c]
    
    exp_data <- Biobase::experimentData(eset)
    
    exp_data@title <- abbrev.c
    
    Biobase::experimentData(eset) <- exp_data
    
    stopifnot(methods::validObject(eset))
    
    abbrev.mset <- MultiDataSet::add_eset(abbrev.mset,
                                          eset,
                                          dataset.type = abbrev.c,
                                          GRanges = NA,
                                          overwrite = TRUE,
                                          warnings = FALSE)
    
  }
  
  abbrev.mset <- abbrev.mset[, ProMetIS::sets_abbrev.vc(full_to_short.l = full_to_short.l)[mset_names.vc]]
  
  invisible(abbrev.mset)
  
}

.sets_abbrev.vc <- function(full_to_short.l = TRUE) {
  
  abbrev.vc <- gsub("metabolomics", "met",
                    gsub("proteomics", "prot",
                         gsub("liver", "liv",
                              gsub("plasma", "plas",
                                   gsub("hyper", "h",
                                        gsub("acqui", "a",
                                             gsub("hilic", "hil",
                                                  gsub("-pos", "-p",
                                                       gsub("-neg", "-n",
                                                            ProMetIS::sets.vc(), fixed = TRUE),
                                                       fixed = TRUE))))))))
  full_to_short.vc <- abbrev.vc
  names(full_to_short.vc) <- ProMetIS::sets.vc()
  
  short_to_full.vc <- ProMetIS::sets.vc()
  names(short_to_full.vc) <- abbrev.vc
  
  if (full_to_short.l) {
    return(full_to_short.vc)
  } else {
    return(short_to_full.vc)
  }
  
}

## Proteomics ----

# Functions from the DAPAR R/Bioconductor package
DAPAR_is.OfType <- function(data, type) 
{
  return(type == data)
}
DAPAR_is.MV <- function(data)
{
  POV = DAPAR_is.OfType(data, "POV")
  MEC = DAPAR_is.OfType(data, "MEC")
  isNA = is.na(data)
  df <- POV | MEC | isNA
  return(df)
}

## Metabolomics ----

.standards <- function() {
  
  standard_file.c <- "//10.0.238.33/Data/Phenostore/data/ProMetIS/cs1_phenomin/0_raw/metabolomics_standards.xlsx"
  standard.df <- as.data.frame(readxl::read_excel(standard_file.c,
                                               sheet = 1,
                                               skip = 4,
                                               n_max = 17),
                            stringsAsFactors = FALSE)
  colnames(standard.df) <- gsub("Maximum shift (mn)", "max_shift_mn",
                             gsub("Mass accuracy (ppm)", "mass_accur_ppm",
                                  gsub("CV* (%)", "cv_percent",
                                       gsub("Rt  (min)", "rt_min",
                                            colnames(standard.df),
                                            fixed = TRUE),
                                       fixed = TRUE),
                                  fixed = TRUE),
                             fixed = TRUE)
  standard_colname_comp.df <- as.data.frame(readxl::read_excel(standard_file.c,
                                                            sheet = 1,
                                                            skip = 3,
                                                            n_max = 2),
                                         stringsAsFactors = FALSE)
  standard_colname_comp.vc <- gsub("Compound", "compound",
                                gsub("Elemental formula", "formula",
                                     gsub("Exact Mass", "mass_exact",
                                          gsub("EI+", "EI_pos",
                                               gsub("EI-", "EI_neg",
                                                    gsub("Concentration", "concentration",
                                                         gsub("Retention time", "rt",
                                                              gsub("Mass measurement", "mass_measured",
                                                                   colnames(standard_colname_comp.df),
                                                                   fixed = TRUE),
                                                              fixed = TRUE),
                                                         fixed = TRUE),
                                                    fixed = TRUE),
                                               fixed = TRUE),
                                          fixed = TRUE),
                                     fixed = TRUE),
                                fixed = TRUE)
  colnames(standard.df) <- paste0(vapply(standard_colname_comp.vc,
                                      function(x) {
                                        unlist(strsplit(x, split = "...", fixed = TRUE))[1]
                                      }, FUN.VALUE = character(1), USE.NAMES = FALSE),
                               vapply(colnames(standard.df),
                                      function(x) {
                                        name.c <- unlist(strsplit(x, split = "...", fixed = TRUE))[1]
                                        if (name.c != "")
                                          name.c <- paste0("@", name.c)
                                        if (name.c %in% c("@max_shift_mn", "@cv_percent"))
                                          name.c <- paste0("rt", name.c)
                                        return(name.c)
                                      }, FUN.VALUE = character(1), USE.NAMES = FALSE))
  colnames(standard.df)[1] <- "standard"
  
  standard.df[, "ei_pos_mz"] <- as.numeric(vapply(standard.df[, "EI_pos"],
                                               function(x) unlist(strsplit(x, split = " "))[1],
                                               FUN.VALUE = character(1)))
  standard.df[, "ei_pos_c8_rt"] <- as.numeric(vapply(standard.df[, 8],
                                                  function(x) unlist(strsplit(x, split = " "))[1],
                                                  FUN.VALUE = character(1))) * 60
  standard.df[, "ei_neg_c8_rt"] <- as.numeric(vapply(standard.df[, 12],
                                                  function(x) unlist(strsplit(x, split = " "))[1],
                                                  FUN.VALUE = character(1))) * 60
  standard.df[, "ei_neg_mz"] <- as.numeric(vapply(standard.df[, "EI_neg"],
                                               function(x) unlist(strsplit(x, split = " "))[1],
                                               FUN.VALUE = character(1)))
  standard.df[, "ei_pos_hilic_rt"] <- as.numeric(vapply(standard.df[, 24],
                                                     function(x) unlist(strsplit(x, split = " "))[1],
                                                     FUN.VALUE = character(1))) * 60
  standard.df[, "ei_neg_hilic_rt"] <- as.numeric(vapply(standard.df[, 28],
                                                     function(x) unlist(strsplit(x, split = " "))[1],
                                                     FUN.VALUE = character(1))) * 60
  
  standard.df <- standard.df[, c("standard",
                           "compound",
                           "formula",
                           "concentration",
                           "mass_exact",
                           "ei_pos_mz",
                           "ei_pos_c8_rt",
                           "ei_pos_hilic_rt",
                           "ei_neg_mz",
                           "ei_neg_c8_rt",
                           "ei_neg_hilic_rt")]
  
  colnames(standard.df)[5:11] <- c("mass",
                                "mz_pos",
                                "rt_pos_c8",
                                "rt_pos_hilic",
                                "mz_neg",
                                "rt_neg_c8",
                                "rt_neg_hilic")
  
  rownames(standard.df) <- standard.df[, "compound"]
  
  standard.df <- standard.df[order(standard.df[, "standard"], standard.df[, "compound"]), ]
  
  standard.df[, "color"] <- rev(rainbow(nrow(standard.df), end = 4/6))
  
  rownames(standard.df) <- gsub("AMPA 2-amino-3-(3-hydroxy-5-methyl-isoxazol-4-yl)propanoic acid",
                             "AMPA",
                             gsub("MCPA 2-methyl-4-chlorophenoxyacetic acid", "MCPA",
                                  rownames(standard.df), fixed = TRUE), fixed = TRUE)
  
  standard.df
  
}
