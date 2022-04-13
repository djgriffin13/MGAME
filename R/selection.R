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
summary.ame_mg <- function(object, ...)
{
  fit<-object
  require(amen)
  tmp<-cbind(apply(fit$BETA,2,mean), apply(fit$BETA,2,sd) ,
             apply(fit$BETA,2,mean)/apply(fit$BETA,2,sd) ,
             2*(1-pnorm( abs(apply(fit$BETA,2,mean)/apply(fit$BETA,2,sd)))))
  colnames(tmp)<-c("pmean","psd","z-stat","p-val")
  cat("\nRegression coefficients:\n")
  print(round(tmp,3))

  cat("\n_________________________________________________\n")

  if(!is.null(fit$GOF)) tmpgof<-apply(fit$GOF,2,mean)
  if(!is.null(fit$GOF) || !.is.null(fit$R2) || !.is.null(fit$R2_nm)) cat("\nModel Fit:\n")
  if(!is.null(fit$R2)) cat("\nR^2 for model is: ", round(fit$R2,3))
  if(!is.null(fit$R2_nm)) cat("\nR^2 compared againset random effects only: ", round(fit$R2_nm,3), "\n\n")

  if(!is.null(fit$GOF)) print(round(tmpgof,3))
  if(!is.null(fit$GOF))cat("\nMinimum index fit:",round(min(tmpgof),3),"\n")


  tmp<-cbind(apply(fit$VC,2,mean), apply(fit$VC,2,sd) )
  colnames(tmp)<-c("pmean","psd")
  cat("\n_________________________________________________\n")

  cat("\nVariance parameters:\n")
  print(round(tmp,3))
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
  observed = mdl$GOF[1, ]
  postGOF = mdl$GOF[-1, ]
  postMeans = colMeans(postGOF)
  observed = abs(observed - postMeans)
  postGOF = abs(postGOF - postMeans)
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
           lastPrior = list(),
           strong = F,
           strongFactor = 1000000) {

    if (length(lastPrior) == 0) {
      lastPrior = list()
      lastPrior$nu0 = 0
      lastPrior$kappa0 = 0
      lastPrior$eta0 = 0
    }

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
    #Degrees of Freedom table for multiplicative effects
    prior$eta0 = length(mdl$APM) + lastPrior$eta0


    ### Dyad Effects ###
    #Dyadic variance
    prior$s20 = mean(mdl$VC[, "ve"], na.rm = T)
    #DF for dyadic variance
    prior$nu0 = nrow(mdl$VC) + lastPrior$nu0

    ### Multiplicative Effects ###
    #Corilation table for multiplicative effects
    prior$Suv0 = var(cbind(mdl$U, mdl$V))
    #Degrees of Freedom table for multiplicative effects
    prior$kappa0 = nrow(mdl$U) + lastPrior$kappa0


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
ameMG = function(data,
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
                 odmax = rep(max(apply(Y > 0, 1, sum, na.rm = TRUE)), nrow(Y)),
                 seed = 1,
                 nscan = 5000,
                 burn = 500,
                 odens = 25,
                 modelFitTest = TRUE,
                 nullModel = TRUE,
                 gof = TRUE,
                 prior = list(),
                 verboseOutput = TRUE,
                 makePlot = FALSE) {
  # R Code: Structure Data by groups

  # Get the list of all groups in the data
  groups = unique(data[[group]])

  # Initialize the places DV and IV's will be stored
  DV = list()
  X_ego = list()
  X_alter = list()
  X_dyad = list()

  # For loop repeats following code for each team
  for (i in 1:length(groups)) {
    # Filter data for the ith team
    g = groups[i]

    edges = dplyr::filter_(data, lazyeval::interp(quote(x == y), x = as.name(group), y =
                                                    g))

    group_size = length(unique(unlist(c(edges[, 1], edges[, 2]))))


    # ### Ego and Alter IV'ss
    X_ego[[i]] = NULL
    for (j in 1:length(Xego)) {
      nodes = dplyr::summarise_(dplyr::group_by_(edges, names(edges)[1]),
                                paste0("mean(", as.name(Xego[j]), ", na.rm = T)"))
      if (j > 1)
        X_ego[[i]] = cbind(X_ego[[i]] , nodes[, 2])
      else
        X_ego[[i]] = nodes[, 2]
    }
    X_ego[[i]] = array(unlist(X_ego[[i]]), dim = c(group_size, length(Xego)))


    X_alter[[i]] = NULL
    for (j in 1:length(Xego)) {
      nodes = dplyr::summarise_(dplyr::group_by_(edges, names(edges)[1]),
                                paste0("mean(", as.name(Xalter[j]), ", na.rm = T)"))
      if (j > 1)
        X_alter[[i]] = cbind(X_alter[[i]] , nodes[, 2])
      else
        X_alter[[i]] = nodes[, 2]
    }
    X_alter[[i]] = array(unlist(X_alter[[i]]), dim = c(group_size, length(Xalter)))

    ### Make team graph object
    G = igraph::graph.data.frame(edges)

    ### Get DV matrix
    DV[[i]] = as.matrix(igraph::get.adjacency(G, attr = Y))

    ### Dyadic IV's ###
    X_dyad[[i]] =  NULL
    for (j in 1:length(Xdyad)) {
      if (j > 1)
        X_dyad[[i]] = cbind(X_dyad[[i]], as.matrix(igraph::get.adjacency(G, attr = Xdyad[i])))
      else
        X_dyad[[i]] = as.matrix(igraph::get.adjacency(G, attr = Xdyad[j]))
    }
    X_dyad[[i]] = array(unlist(X_dyad[[i]]), dim = c(group_size, group_size, length(Xdyad)))
  }

  # R Code: Run models using posterior for each group as prior for next group

  # Clear priors
  null_prior = prior

  # Run code for each team
  for (i in 1:length(groups)) {
    # Run Model:
    ame_model = amen::ame(
      Y = DV[[i]],
      Xrow = X_ego[[i]],
      Xcol = X_alter[[i]],
      Xdyad = X_dyad[[i]],
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
    prior = getPrior(ame_model, prior)

    if (nullModel && modelFitTest) {
      # Run null (Excluding all X's) model used to calculate R^2
      null_model = amen::ame(
        Y = DV[[i]],
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
      null_prior = getPrior(null_model, null_prior)
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
  if (modelFitTest) {
    for (i in 1:length(groups)) {
      #Run model with input for each team
      fit_model = amen::ame(
        DV[[i]],
        Xrow = X_ego[[i]],
        Xcol = X_alter[[i]],
        Xdyad = X_dyad[[i]],
        prior = strongPrior,
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

      # Record goodness of fit information
      GOFPost[i, ] = getAMEFit(fit_model)
      colnames(GOFPost) = colnames(ame_model$GOF)

      # Record Residuals
      residuals = append(residuals, c((fit_model$YPM - DV[[i]]) ^ 2))
      avgResiduals = append(avgResiduals, c((DV[[i]] - mean(DV[[i]], na.rm = T)) ^
                                              2))

      if (nullModel) {
        #Run null model with input for each team
        fit_null = amen::ame(
          DV[[i]],
          prior = strongPrior_null,
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
        res_null = append(res_null, c((fit_null$YPM - DV[[i]]) ^ 2))
      }
    }

    if (modelFitTest) {
      # Calculate Goodness of fit
      ame_model$GOF = GOFPost
      ame_model$GOF_OBT = colMeans(GOFPost)

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
  colnames(ame_model$BETA) = c("(constant)", Xego, Xalter, Xdyad)

  class(ame_model) <- "ame_mg"

  return(ame_model)
}

