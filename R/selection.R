#' Summary of a multi group AME object
#'
#' Summary method for a multi group AME object. This is adapted from
#' Peter Hoff's work for ame class
#'
#' @param object the output from an ameMG function
#' @param ... unused
#' @return a summary of parameter estimates and confidence intervals for a multi-group AME
#'
#' @author Daniel Griffin
#' @export
summary.ame_mg <- function(object, ...) {
  fit <- object
  tmp <- cbind(
    apply(fit$BETA, 2, mean),
    apply(fit$BETA, 2, sd) ,
    apply(fit$BETA, 2, mean) + apply(fit$BETA, 2, sd) * qnorm(.025),
    apply(fit$BETA, 2, mean) + apply(fit$BETA, 2, sd) * qnorm(.025, lower.tail = FALSE),
    apply(fit$BETA, 2, mean) / apply(fit$BETA, 2, sd) ,
    2 * (1 - pnorm(abs(
      apply(fit$BETA, 2, mean) / apply(fit$BETA, 2, sd)
    )))
  )
  colnames(tmp) <-
    c("pmean", "psd", "CI.low95", "CI.up95", "z-stat", "p-val")
  cat("\nRegression coefficients:\n")
  print(round(tmp, 3))


  if (fit$FitOptions[1]) {
    cat("\n______________________________________________________\n")
    tmpgof <- apply(fit$GOF, 2, mean)
  }
  if (fit$FitOptions[1]) {
    cat("\nModel Fit:\n")
    cat(
      "\n - R^2: ",
      round(fit$R2, 3),
      "\nIMPORTANT! This is bisaed for models with additive or multiplicitive random effects.\n"
    )
  }
  if (fit$FitOptions[1] && fit$FitOptions[2])
    cat(
      "\n - Psudo R^2: ",
      round(fit$R2_nm, 3),
      "\nThis accounts for the random effects and compars the full model to the model with the null model that includes the random effects.\n\n"
    )

  if (fit$FitOptions[1] && fit$FitOptions[3]) {
    cat("\n - Goodness of fit:\n")
    cat("First Order Structural Fit:\n")
    print(round(tmpgof[1:2], 3))
    cat("\n- - - - - - - - - -\n")

    cat("Second Order Structural Fit:\n")
    print(round(tmpgof[3], 3))
    cat("\n- - - - - - - - - -\n")

    cat("Third Order Structural Fit:\n")
    print(round(tmpgof[4:5], 3))
    cat("\n- - - - - - - - - -\n")

    cat(
      "SSFI = ",
      round(min(tmpgof), 3),
      "\nThe Structural Selection Fit Index (SSFI) is the minimun fit index.\n"
    )
  }

  tmp <- cbind(apply(fit$VC, 2, mean), apply(fit$VC, 2, sd))
  colnames(tmp) <- c("pmean", "psd")
  cat("\n______________________________________________________\n")

  cat("\nVariance parameters:\n")
  print(round(tmp, 3))

  cat(
    "\nLegend:\n",
    " va  - Ego Random Effect Variance\n",
    " cab - Ego-Alter Random Effect Covariance\n",
    " vb  - Alter Random Effect Variance\n",
    " rho - Dyadic Residual Correlation (a->b compared to b->a)\n",
    " ve  - Residual Variance"
  )
}



#' Get Fit indices for ame models
#'
#' This function takes an ame model output and
#' and provides the porportion of simulations for which
#' the observed value is more distant from the model mean than
#' the simulated vlaues.
#'
#' @param mdl output from ame models of class ame
#' @return a vector of porportions
getAMEFit <- function(mdl) {
  observed = mdl$GOF[1,]
  postGOF = mdl$GOF[-1,]
  postMeans = colMeans(postGOF)
  observed = abs(observed - postMeans)
  postGOF = abs(postGOF - matrix(
    postMeans,
    nrow = nrow(postGOF),
    ncol = ncol(postGOF),
    byrow = TRUE
  ))
  return(colSums(postGOF > observed) / nrow(postGOF))
}

#' Get a Prior for a posterior
#'
#' This function the posterior distribution from an ame model
#' and uses it to generate the prior for the next simulation.
#'
#' @param mdl output from ame models of class ame
#' @param lastPrior previous prior used for tracking degrees of freedom
#' @param strong booliean value indicating weather to treat this as a
#' stribg prior by reducing the variance in beta and increasing the degrees
#' of freedom
#' @param strongFactor value indicating the factor used to augment prior to make it
#' a strong prior.
#' @return a list of hyperparamters used as a prior by the ame model
getPrior <-
  function(mdl,
           strong = F,
           strongFactor = 1000000) {

    prior = list()

    ### Beta Values
    #Precision Matrix for Regression Coefficients
    prior$iV0 = solve(var(mdl$BETA))
    #Mean Regression Coefficients
    prior$m0 = colMeans(mdl$BETA)

    ### Additive Effects ###
    #Corilation table for additive effects
    VC_means = colMeans(mdl$VC)
    prior$Sab0 = matrix(c(VC_means[1], VC_means[2], VC_means[2], VC_means[3]), nrow = 2)
    #Degrees of Freedom table for additive effects (i.e., prior$eta0) set by defaults.

    ### Dyad Effects ###
    #Dyadic variance
    prior$s20 = mean(mdl$VC[, "ve"], na.rm = T)
    #DF for dyadic variance (i.e., prior$nu0) set by defaults

    ### Multiplicative Effects ###
    #Corilation table for multiplicative effects
    prior$Suv0 = var(cbind(mdl$U, mdl$V))
    #Degrees of Freedom table for multiplicative effects (i.e., prior$kappa0) set by defaults.


    # Strong Prior
    if (strong) {
      # Set DF very high so that there is little update
      prior$nu0    = strongFactor
      prior$kappa0 = round(sqrt(strongFactor))
      prior$eta0   = round(sqrt(strongFactor))

      # Set Beta distribution to have very little variance
      prior$iV0 = solve(var(mdl$BETA) / strongFactor)
    }
    return(prior)
  }

#' Multi Group AME model
#'
#' This function runs the ame model for a multi group setting
#' it also calculates fit statistics and variance explained by
#' the model.
#'
#' @param data a data frame in person parwise structure
#' @param Y a dyadic outcome variable of interest
#' @param Xego an individual level anticedent for the ego node
#' @param Xalter an individual level anticedent for the alter node
#' @param family set to "nrm" for weighted dyadic outcom and "bin" for unweighted outcome values.
#' See ame documentation for more information.
#' @param modelFitTest weather to run fit test models
#' @param nullModel whether to run a null_model with only random additive
#' and multiplicitive effects to calculate unique varinace explained by fixed effects.
#' @param verboseOutput print ame model out put or not
#' @param makePlot generate plots as models are running or not
#' @return a list of the class ame_mg which describes the posterior distribution
#' @export
mgame = function(data,
                 Y,
                 group,
                 Xdyad = NULL,
                 Xego = NULL,
                 Xalter = NULL,
                 family,
                 R = 0,
                 rvar = !(family == "rrl") ,
                 cvar = TRUE,
                 dcor = !symmetric,
                 nvar = TRUE,
                 intercept = !is.element(family, c("rrl", "ord")),
                 symmetric = FALSE,
                 seed = 1,
                 nscan = 5000,
                 burn = 500,
                 odens = 25,
                 odmax_grand = NULL,
                 modelFitTest = TRUE,
                 nullModel = TRUE,
                 gof = TRUE,
                 prior = list(),
                 verboseOutput = FALSE,
                 makePlot = FALSE) {
  cat("Running multi-group additive multiplicative effects network model.\n
      This may take a few minutes...")
  # R Code: Structure Data by groups

  # Get the list of all groups in the data
  groups = unique(data[[group]])
  pb <- txtProgressBar(max = length(groups)*(1 + modelFitTest + 2 * modelFitTest * nullModel),style = 3)
  progressIndex = 0
  w = getOption("width")
  # Initialize the places DV and IV's will be stored
  DV = list()
  X_ego = list()
  X_alter = list()
  X_dyad = list()

  cat("\r",rep(" ", w),sep = "")
  cat("\rLoading / Formating Data...\n")
  # For loop repeats following code for each team
  for (i in 1:length(groups)) {
    # Filter data for the ith team
    g = groups[i]

    defaultW <- getOption("warn")
    options(warn = -1)
    edges = dplyr::filter_(data, lazyeval::interp(quote(x == y), x = as.name(group), y =
                                                  g))
    options(warn = defaultW)
    id_list = unique(unlist(edges[,1:2]))
    group_size = length(id_list)
    id_df = data.frame(id_list)

    ### Make team graph object
    dyadMins = rep(NA, length(Xdyad) + 1)
    if (!is.null(Xdyad)) {
      for (j in 1:length(Xdyad)) {
        col = Xdyad[j]
        dyadMins[j] = min(edges[, col])
        edges[, col] = edges[, col] - min(edges[, col]) + 1
      }
    }
    dyadMins[length(Xdyad) + 1] = min(edges[, Y])
    edges[, Y] = edges[, Y] - min(edges[, Y]) + 1


    G = igraph::graph.data.frame(edges)

    ### Get DV matrix
    M = as.matrix(igraph::get.adjacency(G, attr = Y))
    M[M == 0] <- NA
    M = M - 1 + dyadMins[j]
    DV[[i]] = M

    # ### Ego and Alter IV'ss
    X_ego[[i]] = NULL
    if (!is.null(Xego)) {
      names(id_df) = names(edges)[1]
      temp_edges = dplyr::full_join(edges,id_df, by = names(id_df))

      for (j in 1:length(Xego)) {
        defaultW <- getOption("warn")
        options(warn = -1)
        nodes = dplyr::summarise_(dplyr::group_by_(temp_edges, names(edges)[1]),
                                  paste0("mean(", as.name(Xego[j]), ", na.rm = T)"))
        options(warn = defaultW)
        nodes = dplyr::arrange(nodes,by_group =TRUE)
        if (j > 1)
          X_ego[[i]] = cbind(X_ego[[i]] , nodes[, 2])
        else
          X_ego[[i]] = nodes[, 2]
      }
      X_ego[[i]] = array(unlist(X_ego[[i]]), dim = c(group_size, length(Xego)))
    }

    X_alter[[i]] = NULL
    if (!is.null(Xalter)) {
      names(id_df) = names(edges)[2]
      temp_edges = dplyr::full_join(edges,id_df, by = names(id_df))

      for (j in 1:length(Xalter)) {
        defaultW <- getOption("warn")
        options(warn = -1)
        nodes = dplyr::summarise_(dplyr::group_by_(temp_edges, names(edges)[2]),
                                  paste0("mean(", as.name(Xalter[j]), ", na.rm = T)"))
        options(warn = defaultW)
        nodes = dplyr::arrange(nodes,by_group =TRUE)
        if (j > 1)
          X_alter[[i]] = cbind(X_alter[[i]] , nodes[, 2])
        else
          X_alter[[i]] = nodes[, 2]
      }
      X_alter[[i]] = array(unlist(X_alter[[i]]), dim = c(group_size, length(Xalter)))
    }

    ### Dyadic IV's ###
    X_dyad[[i]] =  NULL
    if (!is.null(Xdyad)) {
      for (j in 1:length(Xdyad)) {
        M = as.matrix(igraph::get.adjacency(G, attr = Xdyad[j]))
        M[M == 0] <- NA
        M = M - 1 + dyadMins[j]
        M = M[order(id_list), order(id_list)]
        if (j > 1)
          X_dyad[[i]] = cbind(X_dyad[[i]], M)
        else
          X_dyad[[i]] = M
      }
      X_dyad[[i]] = array(unlist(X_dyad[[i]]), dim = c(group_size, group_size, length(Xdyad)))
    }
  }
  # R Code: Run models using posterior for each group as prior for next group

  # Clear priors
  null_prior = prior

  # Run code for each team
  cat("Running Models...\n\n")
  for (i in 1:length(groups)) {
    # Run Model:
    cat("\r",rep(" ", w),sep = "")
    cat("\r -  Running main ame model for Group",i,"\n")
    setTxtProgressBar(pb, progressIndex+.1)
    if (is.null(odmax_grand))
      ame_model = amen::ame(
        Y = DV[[i]],
        Xrow  = if (is.null(Xego))   NULL else X_ego[[i]],
        Xcol  = if (is.null(Xalter)) NULL else X_alter[[i]],
        Xdyad = if (is.null(Xdyad))  NULL else X_dyad[[i]],
        R = R,
        prior = prior,
        family = family,
        rvar = rvar,
        cvar = cvar,
        dcor = dcor,
        nvar = nvar,
        intercept = intercept,
        symmetric = symmetric,
        seed = seed,
        nscan = nscan,
        burn = burn,
        odens = odens,
        plot = makePlot,
        print = verboseOutput,
        gof = gof
      )
    else
      ame_model = amen::ame(
        Y = DV[[i]],
        Xrow  = if (is.null(Xego))   NULL else X_ego[[i]],
        Xcol  = if (is.null(Xalter)) NULL else X_alter[[i]],
        Xdyad = if (is.null(Xdyad))  NULL else X_dyad[[i]],
        R = R,
        prior = prior,
        family = family,
        rvar = rvar,
        cvar = cvar,
        dcor = dcor,
        nvar = nvar,
        intercept = intercept,
        symmetric = symmetric,
        odmax = odmax,
        seed = seed,
        nscan = nscan,
        burn = burn,
        odens = odens,
        plot = makePlot,
        print = verboseOutput,
        gof = gof
      )
    progressIndex = progressIndex + 1
    prior = getPrior(ame_model)

    if (nullModel && modelFitTest) {
      cat("\r",rep(" ", w),sep = "")
      cat("\r -  Running null ame model for comparison for Group:",i,"\n")
      setTxtProgressBar(pb, progressIndex)
      
      # Run null (Excluding all X's) model used to calculate R^2
      if (is.null(odmax_grand))
        null_model = amen::ame(
          Y = DV[[i]],
          R = R,
          prior = null_prior,
          family = family,
          rvar = rvar,
          cvar = cvar,
          dcor = dcor,
          nvar = nvar,
          intercept = intercept,
          symmetric = symmetric,
          seed = seed,
          nscan = nscan,
          burn = burn,
          odens = odens,
          plot = makePlot,
          print = verboseOutput,
          gof = gof
        )
      else
        null_model = amen::ame(
          Y = DV[[i]],
          R = R,
          prior = null_prior,
          family = family,
          rvar = rvar,
          cvar = cvar,
          dcor = dcor,
          nvar = nvar,
          intercept = intercept,
          symmetric = symmetric,
          odmax = odmax,
          seed = seed,
          nscan = nscan,
          burn = burn,
          odens = odens,
          plot = makePlot,
          print = verboseOutput,
          gof = gof
        )
      progressIndex = progressIndex + 1
      null_prior = getPrior(null_model)
    }

    # Get prior from output posterior distributions
  }

  # R Code: Assess Fit and R^2 of Model
  residuals = c()
  avgResiduals = c()
  res_null = c()
  GOFPost = matrix(NA, nrow = length(groups), ncol = 5)

  # Get a strong prior set for the model to strongly prefer so that model will not update
  # when ran to collect information for fit and R^2 values.
  strongPrior = getPrior(ame_model, strong = TRUE)
  if (nullModel &&
      modelFitTest)
    strongPrior_null = getPrior(null_model, strong = TRUE)

  # Loop through all groups again, but always using the same strong prior so that the model
  # give us residuals and fit information but doesn't update
  cat("\r",rep(" ", w),sep = "")
  cat("\nCalculating  Model Fit and Structural Fit...\n")
  if (modelFitTest) {
    for (i in 1:length(groups)) {
      #Run model with input for each team
      cat("\r",rep(" ", w),sep = "")
      cat("\r -  Main fit for group:",i,"\n")
      setTxtProgressBar(pb, progressIndex)
      
      if (is.null(odmax_grand))
        fit_model = amen::ame(
          DV[[i]],
          Xrow  = if (is.null(Xego))   NULL else X_ego[[i]],
          Xcol  = if (is.null(Xalter)) NULL else X_alter[[i]],
          Xdyad = if (is.null(Xdyad))  NULL else X_dyad[[i]],
          prior = strongPrior,
          R = R,
          burn = 0,
          family = family,
          rvar = rvar,
          cvar = cvar,
          dcor = dcor,
          nvar = nvar,
          intercept = intercept,
          symmetric = symmetric,
          seed = seed,
          nscan = nscan,
          odens = odens,
          plot = makePlot,
          print = verboseOutput,
          gof = gof
        )
      else
        fit_model = amen::ame(
          DV[[i]],
          Xrow  = if (is.null(Xego))   NULL else X_ego[[i]],
          Xcol  = if (is.null(Xalter)) NULL else X_alter[[i]],
          Xdyad = if (is.null(Xdyad))  NULL else X_dyad[[i]],
          prior = strongPrior,
          R = R,
          burn = 0,
          family = family,
          rvar = rvar,
          cvar = cvar,
          dcor = dcor,
          nvar = nvar,
          intercept = intercept,
          symmetric = symmetric,
          odmax = odmax,
          seed = seed,
          nscan = nscan,
          odens = odens,
          plot = makePlot,
          print = verboseOutput,
          gof = gof
        )
      progressIndex = progressIndex + 1
      
      # Record goodness of fit information
      GOFPost[i,] = getAMEFit(fit_model)

      # Record Residuals
      residuals = append(residuals, c((fit_model$YPM - DV[[i]]) ^ 2))
      avgResiduals = append(avgResiduals, c((DV[[i]] - mean(DV[[i]], na.rm = T)) ^
                                              2))

      if (nullModel) {
        cat("\r",rep(" ", w),sep = "")
        cat("\r -  Null model fit for group:",i,"\n")
        setTxtProgressBar(pb, progressIndex)
        
        #Run null model with input for each team
        if (is.null(odmax_grand))
          fit_null = amen::ame(
            DV[[i]],
            prior = strongPrior_null,
            R = R,
            family = family,
            burn = 0,
            rvar = rvar,
            cvar = cvar,
            dcor = dcor,
            nvar = nvar,
            intercept = intercept,
            symmetric = symmetric,
            seed = seed,
            nscan = nscan,
            odens = odens,
            plot = makePlot,
            print = verboseOutput,
            gof = gof
          )
        else
          fit_null = amen::ame(
            DV[[i]],
            prior = strongPrior_null,
            R = R,
            family = family,
            burn = 0,
            rvar = rvar,
            cvar = cvar,
            dcor = dcor,
            nvar = nvar,
            intercept = intercept,
            symmetric = symmetric,
            odmax = odmax,
            seed = seed,
            nscan = nscan,
            odens = odens,
            plot = makePlot,
            print = verboseOutput,
            gof = gof
          )

        progressIndex = progressIndex + 1
        res_null = append(res_null, c((fit_null$YPM - DV[[i]]) ^ 2))
      }
    }
    cat("\r",rep(" ", w),sep = "")
    cat("\rFinalizing Results...\n")
    setTxtProgressBar(pb, progressIndex)
    
    if (modelFitTest) {
      # Calculate Goodness of fit
      ame_model$GOF = GOFPost

      # Calculate basic R^2: variance explained by the full model over the intercept only model
      r2 = 1 - sum(residuals, na.rm = TRUE) / sum(avgResiduals, na.rm = TRUE)
      ame_model$R2 = r2
    }
    if (nullModel) {
      # Calculate fixed effects R^2: variance explained by full model over the random effects only model
      r2_nm = 1 - sum(residuals, na.rm = TRUE) / sum(res_null, na.rm = TRUE)
      ame_model$R2_nm = r2_nm
    }
  }
  # R Code: Get model regression parameters and variance components
  colnames(ame_model$BETA) = c(
    "(constant)",
    if(!is.null(Xego)) paste0(Xego, "_(ego)"),
    if(!is.null(Xalter)) paste0(Xalter, "_(alter)"),
    if(!is.null(Xdyad))paste0(Xdyad, "_(dyad)")
  )
  colnames(ame_model$GOF) = c("Ego",
                              "Alter",
                              "Reciprocity",
                              "Cycle Closure",
                              "Transitivity")

  ame_model$FitOptions = c(modelFitTest, nullModel, gof)

  class(ame_model) <- "ame_mg"

  return(ame_model)
}

