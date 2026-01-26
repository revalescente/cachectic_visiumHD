library(purrr)
library(SpatialExperiment)
library(scuttle)
library(SpaceTrooper)

# custom spatialQC metrics ----

vHD_spatialPerCellQC <- function(spe, micronConvFact=0.316, rmZeros=TRUE, ...) {
  stopifnot(is(object=spe, "SpatialExperiment"))
  
  is.mito <- BiocGenerics::grep("^mt-", BiocGenerics::rownames(spe), ignore.case = TRUE)
  spe <- scuttle::addPerCellQC(spe, subsets = list(Mito = is.mito), ...)
  
  if(!all(SpatialExperiment::spatialCoordsNames(spe) %in% names(colData(spe)))) {
    colData(spe) <- cbind.DataFrame(colData(spe), spatialCoords(spe))
  }
  
  spnc <- spatialCoords(spe) * micronConvFact
  colnames(spnc) <-  paste0(SpatialExperiment::spatialCoordsNames(spe), "_um")
  colData(spe) <- cbind.DataFrame(colData(spe), spnc)
  spe$Area_um <- spe$area * (micronConvFact^2)
  spe$AspectRatio <- spe$major_axis_length/spe$minor_axis_length
  if ("AspectRatio" %in% colnames(colData(spe))) {
    spe$log2AspectRatio <- log2(spe$AspectRatio) # not cosmx
  } else { warning("Missing aspect ratio in colData") }
  
  spe$CountArea <- spe$sum/spe$Area_um
  spe$log2CountArea <- log2(spe$CountArea)
  # spe$SignalDensity <- spe$sum/spe$Area_um # these are the new names of the CountArea
  # spe$log2SignalDensity <- log2(spe$SignalDensity)
  if (rmZeros) {
    if (sum(spe$sum==0) > 0) {
      message("Removing ", dim(spe[,spe$sum==0])[2],
              " cells with 0 counts!")
      spe <- spe[,!spe$sum==0]
    }
  }
  return(spe)
}

# Helper function of SpaceTrooper

.prepQCContext <- function(spe, metricList=c("log2SignalDensity", "Area_um",
                                             "log2AspectRatio"), verbose=FALSE) {
  
  stopifnot(is(spe,"SpatialExperiment"))
  
  # remove zero-count cells once
  zerocells <- spe$total==0
  if (sum(zerocells) > 0) {
    warning(paste0(sum(zerocells),
                   " cells with 0 counts were found. These cells will be removed."))
    spe <- spe[, sum(zerocells)]
  }
  spe1 <- computeOutliersQCScore(spe, metricList)
  spe1 <- checkOutliers(spe1, verbose)
  
  out_var <- metadata(spe1)$formula_variables
  df <- as.data.frame(colData(spe1))
  tech <- metadata(spe1)$technology
  
  return(list(df=df, out_var=out_var, tech=tech))
}

.checkSkw <- function(cd, metricList=c("log2SignalDensity", "Area_um",
                                       "log2AspectRatio")) {
  stopifnot(is(cd, "DataFrame"))
  stopifnot(all(metricList %in% names(cd)))
  method <- c()
  for (i in metricList) {
    skw <- e1071::skewness(cd[[i]], na.rm = TRUE)
    method[i] <- ifelse((skw>-1 & skw<1), "sc", "mc")
  }
  for (submeth in c("log2SignalDensity", "log2Ctrl_total_ratio"))
  {
    if (!submeth %in% names(method)) next
    idx <- switch(submeth,
                  "log2SignalDensity"        = cd[["total"]] > 0,
                  "log2Ctrl_total_ratio" = cd[["ctrl_total_ratio"]] != 0,
                  rep(TRUE, nrow(cd))
    )
    if (!any(idx, na.rm = TRUE)) next
    skw <- e1071::skewness(cd[[submeth]][idx], na.rm=TRUE)
    newm <- ifelse((skw>-1 & skw<1), "sc", "mc")
    if (!is.na(newm)) method[[submeth]] <- newm
  }
  return(method)
}

computeOutliersQCScore <- function(spe, metricList=c("log2SignalDensity","Area_um",
                                                     "log2AspectRatio")) {
  
  stopifnot(is(spe, "SpatialExperiment"))
  
  cd <- colData(spe)
  method <- .checkSkw(cd, metricList)
  
  spec <- list(
    log2SignalDensity = list(
      idx  = spe$total > 0,
      zero = spe$total == 0,
      tweak_lower = TRUE
    ),
    log2Ctrl_total_ratio = list(
      idx  = spe$ctrl_total_ratio != 0,
      zero = spe$ctrl_total_ratio == 0,
      tweak_lower = FALSE
    )
  )
  if(any(is.na(method)))
  {
    idx <- which(is.na(method))
    for (i in idx) {
      message("Metric ", names(method)[i], " produced NA skewness values.",
              "It will be removed from the QC score formula.")
      method <- method[-i]
    }
  }
  for (var in intersect(names(method), names(spec))) {
    s <- spec[[var]]
    spe_temp <- computeSpatialOutlier(spe[, s$idx], computeBy=var,
                                      method=method[[var]])
    out_var <- paste0(var, "_outlier_", method[[var]])
    low_thr  <- getFencesOutlier(spe_temp, out_var, "lower")
    high_thr <- getFencesOutlier(spe_temp, out_var, "higher")
    if (isTRUE(s$tweak_lower)) {
      if (low_thr < min(spe_temp[[var]], na.rm = TRUE)) {
        low_thr <- stats::quantile(spe[[var]], probs=0.01)
      }
    }
    tgt <- paste0(var, "_outlier_train")
    spe[[tgt]] <- dplyr::case_when(
      s$zero ~ "NO",
      spe[[var]] < low_thr  ~ "LOW",
      spe[[var]] > high_thr ~ "HIGH",
      TRUE ~ "NO"
    )
    spe[[tgt]] <- scuttle::outlier.filter(spe[[tgt]])#is this really needed?
    thr <- getFencesOutlier(spe_temp, out_var)
    if (isTRUE(s$tweak_lower)) thr[1] <- low_thr
    attr(spe[[tgt]], "thresholds") <- thr
  }
  
  submethod <- method[!names(method) %in%
                        c("log2SignalDensity", "log2Ctrl_total_ratio")]
  
  for(j in names(submethod)){
    spe <- computeSpatialOutlier(spe, computeBy=j, method=submethod[j])
  }
  
  out_var <- paste0(names(method), "_outlier_", method)
  names(out_var) <- names(method)
  # gives warning if one of the variables is missing, but still works!
  # out_var[names(out_var) %in% c("log2SignalDensity", "log2Ctrl_total_ratio")] <-
  #     c("log2SignalDensity_outlier_train", "log2Ctrl_total_ratio_outlier_train")
  mapping <- c(
    log2SignalDensity = "log2SignalDensity_outlier_train",
    log2Ctrl_total_ratio = "log2Ctrl_total_ratio_outlier_train"
  )
  present <- intersect(names(out_var), names(mapping))
  if (length(present) > 0L) {
    for (nm in present) out_var[nm] <- mapping[nm]
  }
  metadata(spe)$formula_variables <- out_var
  return(spe)
}

checkOutliers <- function(spe, verbose=FALSE) {
  warnstopmsg <- function(var, warnstop=c("w","s")) {
    warnstop <- match.arg(warnstop)
    m1 <- paste0("Not enough outlier cells for ", var, ".\n")
    m2 <- switch(warnstop,
                 s="QC score computation cannot be performed",
                 w="This variable will not be used in the final formula")
    return(paste0(m1, m2))
  }
  out_var <- metadata(spe)$formula_variables
  cd <- colData(spe)
  if (verbose) {
    for (i in names(out_var)) {
      message("Outliers found for ", i, ":")
      message(paste(table(cd[[out_var[i]]]), collapse = " "))
    }
  }
  stopifnot(
    "log2SignalDensity is not included in the QC score formula.\n
        QC score cannot be computed"=
      "log2SignalDensity" %in% names(out_var)
  )
  cfg <- list(
    log2SignalDensity=list(pattern="log2SignalDensity_outlier",
                           remove="log2SignalDensity_outlier_train", label="LOW", act=stop,
                           code="s"),
    Area_um=list(pattern="Area_um_outlier", remove="Area_um_outlier",
                 label="HIGH", act=warning, code="w"),
    log2Ctrl_total_ratio=list(pattern="log2Ctrl_total_ratio_outlier",
                              remove="log2Ctrl_total_ratio_outlier_train", label="HIGH",
                              act=warning, code="w")
  )
  for (v in intersect(names(cfg), names(out_var))) {
    r <- cfg[[v]]
    gi <- grep(r$pattern, out_var)
    if (length(gi)==0L) next
    col <- out_var[gi]
    labs <- factor(cd[[col]], levels=c("LOW","HIGH","NO"))
    tab <- table(labs)
    cnt <- tab[r$label]
    if (is.null(cnt) || is.na(cnt)) cnt <- 0L
    if (cnt < ncol(spe)*0.001) {
      r$act(warnstopmsg(v, r$code))
      out_var <- out_var[-grep(r$remove, out_var)]
    }
  }
  var <- "log2AspectRatio"
  pat <- paste0(var, "_outlier")
  is_cosmx <- metadata(spe)$technology %in%
    c("Nanostring_CosMx", "Nanostring_CosMx_Protein")
  if (is_cosmx && (var %in% names(out_var))) {
    idx <- grep(pat, out_var)
    if (length(idx)) {
      col <- out_var[idx]
      labs <- factor(cd[[col]], levels=c("LOW","HIGH","NO"))
      tab <- table(labs)
      nmin <- ncol(spe) * 0.001
      
      low  <- tab[["LOW"]];  if (is.na(low))  low  <- 0L
      high <- tab[["HIGH"]]; if (is.na(high)) high <- 0L
      
      if (low < nmin && high < nmin) {
        warning(warnstopmsg(var, "w"))
        out_var <- out_var[-grep(pat, out_var)]
      }
    }
  } else {
    out_var <- out_var[-grep(pat, out_var)]
  }
  metadata(spe)$formula_variables <- out_var
  return(spe)
}



# function to compute the Quality Control scores for each nucleus

vHD_computeQCScore <- function (spe, bestLambda = NULL, verbose = FALSE) 
{
  stopifnot(is(spe, "SpatialExperiment"))
  if (dim(spe[, spe$total == 0])[2] != 0) {
    warning(paste0(dim(spe[, spe$total == 0])[2], " cells with 0 counts were found. These cells will be removed."))
    spe <- spe[, spe$total > 0]
  }
  metricList <- c("log2CountArea", "Area_um", "log2AspectRatio")
  stopifnot(`Not all required metrics in the colData.\nPlease run spatialPerCellQC first.` = all(metricList %in% 
                                                                                                   names(colData(spe))))
  ctx <- .prepQCContext(spe, metricList, verbose)
  df <- ctx$df
  out_var <- ctx$out_var
  tech <- ctx$tech
  train_df <- computeTrainDF(df, out_var, tech, verbose)
  model_formula <- getModelFormula(out_var, verbose)
  model_matrix <- model.matrix(as.formula(model_formula), 
                               data = train_df)
  model <- trainModel(model_matrix, train_df)
  if (is.null(bestLambda)) {
    bestLambda <- computeLambda(train_df, model_formula)
  }
  if (verbose) {
    message("Model coefficients for every term used in the formula:")
    message(round(predict(model, s = bestLambda, type = "coefficients"), 
                  2))
  }
  full_matrix <- model.matrix(as.formula(model_formula), data = df)
  cd <- colData(spe)
  cd$QC_score <- as.vector(predict(model, s = bestLambda, 
                                   newx = full_matrix, type = "response"))
  spe$QC_score <- cd$QC_score
  return(spe)
}