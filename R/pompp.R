#' @include initial-class.R prior-class.R model-class.R fit-class.R
NULL

#' Fit presence-only data using the Presence-Only for Marked Point Process model
#'
#' The model uses a data augmentation scheme to avoid performing approximations
#' on the likelihood function.
#' @param object Either a \code{pompp_model} or \code{pompp_fit} object. If
#' a model, then the model is fit according to specifications. If a fit,
#' then the model used to fit the model is recovered and used to continue
#' the MCMC calculations where the previous one left off.
#' @param background A matrix where the rows are the grid cells for the studied
#' region and the columns are the covariates. \code{NA}s must be removed. If
#' the function is being used on a \code{pompp_fit} object, the background
#' must be exactly the same as the one used in the original fit.
#' @param mcmc_setup A list containing \code{iter} to inform the model how
#' many iterations are to be run. The list may optionally contain the objects.
#' @param verbose Set to \code{FALSE} to suppress all messages to console.
#' @param ... Parameters passed on to specific methods.
#' \code{burnin} and \code{thin} to inform these instructions as well.
#' @return An object of class \code{"pompp_fit"}.
#' @details The background is kept outside of the
#' @seealso \code{\link{pompp_model}} and \code{\link{pompp_fit-class}}.
#' @examples
#' \donttest{
#' beta <- c(-1, 2) # Intercept = -1. Only one covariate
#' delta <- c(3, 4) # Intercept = 3. Only one covariate
#' lambdaStar <- 1000
#' gamma <- 2
#' mu <- 5
#'
#' total_points <- rpois(1, lambdaStar)
#' random_points <- cbind(runif(total_points), runif(total_points))
#' grid_size <- 50
#'
#' # Find covariate values to explain the species occurrence.
#' # We give them a Gaussian spatial structure.
#' reg_grid <- as.matrix(expand.grid(seq(0, 1, len = grid_size), seq(0, 1, len = grid_size)))
#' Z <- MASS::mvrnorm(1, rep(0, total_points + grid_size * grid_size),
#'   3 * exp(-as.matrix(dist(rbind(random_points, reg_grid))) / 0.2))
#' Z1 <- Z[1:total_points]; Z2 <- Z[-(1:total_points)]
#'
#' # Thin the points by comparing the retaining probabilities with uniforms
#' # in the log scale to find the occurrences
#' occurrences <- log(runif(total_points)) <= -log1p(exp(-beta[1] - beta[2] * Z1))
#' n_occurrences <- sum(occurrences)
#' occurrences_points <- random_points[occurrences,]
#' occurrences_Z <- Z1[occurrences]
#'
#' # Find covariate values to explain the observation bias.
#' # Additionally create a regular grid to plot the covariate later.
#' W <- MASS::mvrnorm(1, rep(0, n_occurrences + grid_size * grid_size),
#'   2 * exp(-as.matrix(dist(rbind(occurrences_points, reg_grid))) / 0.3))
#' W1 <- W[1:n_occurrences]; W2 <- W[-(1:n_occurrences)]
#' S <- MASS::mvrnorm(1, rep(0, n_occurrences),
#'   2 * exp(-as.matrix(dist(occurrences_points)) / 0.3))
#'
#' # Find the presence-only observations.
#' marks <- exp(mu + S + rnorm(n_occurrences, 0, 0.3))
#' po_sightings <- log(runif(n_occurrences)) <= -log1p(exp(-delta[1] - delta[2] * W1 - gamma * S))
#' n_po <- sum(po_sightings)
#' po_points <- occurrences_points[po_sightings, ]
#' po_Z <- occurrences_Z[po_sightings]
#' po_W <- W1[po_sightings]
#' po_marks <- marks[po_sightings]
#'
#' jointPrior <- prior(
#'   NormalPrior(rep(0, 2), 10 * diag(2)), # Beta
#' NormalPrior(rep(0, 3), 10 * diag(3)), # Delta
#' GammaPrior(0.00001, 0.00001), # LambdaStar
#' NormalPrior(0, 100), GammaPrior(0.001, 0.001) # Marks
#' )
#'
#' model <- pompp_model(po = cbind(po_Z, po_W, po_points, po_marks),
#' intensitySelection = 1, observabilitySelection = 2, marksSelection = 5,
#'                     coordinates = 3:4,
#'                     intensityLink = "logit", observabilityLink = "logit",
#'                     initial_values = 2, joint_prior = jointPrior)
#'
#' bkg <- cbind(Z2, W2, reg_grid) # Create background
#'
#' # Be prepared to wait a long time for this
#' fit <- fit_pompp(model, bkg, neighborhoodSize = 20, area = 1,
#'   mcmc_setup = list(burnin = 1000, iter = 2000), cores = 1)
#'
#' summary(fit)
#'
#' # Rhat upper CI values are above 1.1. More iterations are needed...
#' }
#' @importFrom coda mcmc mcmc.list
#' @export
methods::setGeneric("fit_pompp", function(object, background,
                                            mcmc_setup = list(iter = 5000),
                                            verbose = TRUE, ...)
  standardGeneric("fit_pompp"))

#' @rdname fit_pompp
#' @param area A positive number with the studied region's area.
#' @param neighborhoodSize Number of neighbors to use in the NNGP method.
#' @param cores Number of cores to pass to OpenMP.
#' @importFrom parallel detectCores
#' @importFrom geoR likfit
methods::setMethod("fit_pompp",
                   methods::signature(object = "pompp_model",
                                      background = "matrix"),
                   function(object, background, neighborhoodSize = 20, mcmc_setup, verbose = TRUE, area = 1, cores = parallel::detectCores(), ...){
                     ## Verifying background names if columns are selected by column name. Crewating background selection variables
                     backConfig <- checkFormatBackground(object, background)

                     # Helper function
                     s <- function(n) methods::slot(object, n)

                     ## Verifying area
                     stopifnot(length(area) == 1, area > 0)
                     ## Verifying mcmc_setup - begin
                     stopifnot(is.list(mcmc_setup))
                     if (is.null(mcmc_setup$burnin)) mcmc_setup$burnin <- 0
                     if (is.null(mcmc_setup$thin)) mcmc_setup$thin <- 1
                     mcmc_setup$burnin <- as.numeric(mcmc_setup$burnin)
                     mcmc_setup$thin <- as.numeric(mcmc_setup$thin)
                     mcmc_setup$iter <- as.numeric(mcmc_setup$iter)
                     stopifnot("iter" %in% names(mcmc_setup), !is.na(mcmc_setup$burnin),
                               length(mcmc_setup$burnin) == 1, mcmc_setup$burnin == floor(mcmc_setup$burnin),
                               mcmc_setup$burnin >= 0, !is.na(mcmc_setup$thin), length(mcmc_setup$thin) == 1,
                               mcmc_setup$thin == floor(mcmc_setup$thin), mcmc_setup$thin > 0,
                               !is.na(mcmc_setup$iter), length(mcmc_setup$iter) == 1,
                               mcmc_setup$iter == floor(mcmc_setup$iter), mcmc_setup$iter > 0,
                               mcmc_setup$iter >= mcmc_setup$thin,
                               cores > 0, cores == floor(cores), length(cores) == 1,
                               cores <= parallel::detectCores())
                     ## Verifying mcmc_setup - end

                     # Helper parameters
                     nb <- length(s("intensitySelection")) + 1
                     nd <- length(s("observabilitySelection")) + 1
                     npar <- nb + nd + 8 - 1 # -1 because log posterior is not working yet.
                     chains <- length(s("init"))

                     # Set up parameter names to construct the coda::mcmc objects
                     if (length(s("iSelectedColumns") > 0))
                       betanames <- c("(Intensity intercept)",
                                      s("iSelectedColumns"))
                     else betanames = paste0("beta_", 1:nb - 1)
                     if (length(s("oSelectedColumns") > 0))
                       deltanames <- c("(Observability intercept)",
                                       s("oSelectedColumns"))
                     else deltanames = paste0("delta_", 1:nd - 1)
                     deltanames <- c(deltanames, "gamma")
                     parnames <- c(betanames, deltanames,
                                   "lambdaStar", "marksMean", "marksNugget",
                                   "n_U", "n_Xprime", "sumLatentMarks",
                                   "varianceLatentMarks",
                                   "allMarksSum", "allMarksVariance"
#                                   ,"log_Posterior" # Not working yet.
                     )

                     ## Determining GP hyperparameters
                     if (verbose) cat("Preparing model parameters.\n")
                     fitted <- geoR::likfit(
                       geoR::as.geodata(s("po")[,
                         c(s("coordinates"), s("marksSelection"))
                       ]), lambda = 0, ini.cov.pars = c(1, 1),
                       cov.model = "exponential", fix.nugget = FALSE,
                       messages = FALSE
                     )

                     time <- Sys.time()
                     mcmcRun <- list()
                     heatMap <- rep(0, nrow(background))
                     for (c in 1:chains){
                       if (chains > 1 && verbose) cat("Starting chain ", c, ".\n",sep="")
                       temp <- cppPOMPP(
                         methods::slot(s("init")[[c]], "beta"),
                         methods::slot(s("init")[[c]], "delta"), # Init delta
                         methods::slot(s("init")[[c]], "lambdaStar"), # Init lambdaStar
                         "",
                         "",
                         "",
                         retrievePars(methods::slot( # beta prior parameters
                           s("prior"),"beta")),
                         retrievePars(methods::slot( # delta prior parameters
                           s("prior"),"delta")),
                         methods::slot(s("prior"), "lambdaStar")$shape,
                         methods::slot(s("prior"), "lambdaStar")$rate,
                         "",
                         background, # background data
                         area, # region area
                         "",
                         methods::slot(s("init")[[c]], "marksMean"),
                         methods::slot(s("init")[[c]], "marksPrecision"),
                         methods::slot(s("prior"), "marksMean")$mu,
                         methods::slot(s("prior"), "marksMean")$Sigma,
                         methods::slot(s("prior"), "marksPrecision")$shape,
                         methods::slot(s("prior"), "marksPrecision")$rate,
                         s("po"), s("po")[, s("marksSelection")], # po data
                         s("po")[, s("coordinates")],
                         backConfig$bIS - 1, # background intensity columns
                         backConfig$bOS - 1, # background observability colmns
                         s("intensitySelection") - 1, # po intensity columns
                         s("observabilitySelection") - 1, # po obserability columns
                         1, fitted$sigmasq, fitted$phi,
                         neighborhoodSize,
                         s("coordinates")[1], s("coordinates")[2],
                         mcmc_setup$burnin, # MCMC burn-in
                         mcmc_setup$thin, # MCMC thin
                         mcmc_setup$iter, # MCMC iterations
                         cores, verbose)
                       #heatMap <- heatMap + temp$xPrimePred # Not available yet
                       mcmcRun[[c]] <- do.call(cbind, temp[-length(temp)])
                       colnames(mcmcRun[[c]]) <- parnames
                       mcmcRun[[c]] <- coda::mcmc(mcmcRun[[c]], thin = mcmc_setup$thin)
                       if (chains > 1 && verbose) cat("Finished chain ",c,".\n\n",sep="")
                     }
                     if (chains > 1 && verbose) cat("Total computation time: ", format(unclass(Sys.time()-time),
                                                                                       digits = 2), " ",
                                                    attr(Sys.time() - time, "units"), ".\n", sep="")

                     return(methods::new("pompp_fit",
                                         fit = do.call(coda::mcmc.list, mcmcRun),
                                         # heatMap = heatMap,
                                         original = object,
                                         neighborhoodSize = neighborhoodSize,
                                         backgroundSummary = summary(background),
                                         area = area,
                                         parnames = parnames,
                                         mcmc_setup = mcmc_setup))
                   })


##' @rdname fit_pompp
##' @param cores Number of cores to pass to OpenMP.
# methods::setMethod("fit_pompp", signature(object = "pompp_fit",
#                                             background = "matrix"),
#                    function(object, background, mcmc_setup = list(iter = object$mcmc_setup$iter), verbose = TRUE,
#                             cores = 1, ...){
#                      # Helper function
#                      s <- function(n) methods::slot(object, n)
#                      so <- function(n) methods::slot(s("original"), n)
#                      backConfig <- checkFormatBackground(s("original"), background)
#
#                      # Check background differences
#                      if (verbose) cat("Checking if everything's ok...")
#                      stopifnot(all.equal(s("backgroundSummary"), summary(background)),
#                                "iter" %in% names(mcmc_setup), !is.na(mcmc_setup$iter),
#                                length(mcmc_setup$iter) == 1,
#                                mcmc_setup$iter == floor(mcmc_setup$iter), mcmc_setup$iter > 0,
#                                mcmc_setup$iter >= mcmc_setup$thin,
#                                cores > 0, cores == floor(cores), length(cores) == 1,
#                                cores <= parallel::detectCores())
#                      if ("thin" %in% names(mcmc_setup))
#                        stopifnot(mcmc_setup$thin == s("mcmc_setup")$thin)
#                      else
#                        mcmc_setup$thin = s("mcmc_setup")$thin
#                      if ("burnin" %in% names(mcmc_setup) && mcmc_setup$burnin > 0)
#                        warning("\nBurnin is disabled when continuing MCMC procedure.")
#                      else
#                        cat(" Everything seems ok.\n")
#
#                      # Helper parameters
#                      betaPos <- 1:length(methods::slot(
#                        methods::slot(so("prior"), "beta"), "mu"))
#                      deltaPos <- (max(betaPos) + 1):(max(betaPos) +
#                                                        length(methods::slot(
#                                                          methods::slot(so("prior"),
#                                                                        "delta"), "mu")))
#                      lambdaPos <- max(deltaPos) + 1
#                      chains <- length(s("fit"))
#                      lastPoint <- nrow(s("fit")[[1]])
#                      lastPoint <- s("fit")[[1]][lastPoint, ]
#
#                      time <- Sys.time()
#                      mcmcRun <- list()
#                      heatMap <- s("heatMap")
#                      for (c in 1:chains){
#                        if (chains > 1 && verbose) cat("Starting chain ",c,".\n",sep="")
#                        temp <- cppPOMPP(
#                          lastPoint[betaPos], # Starting at last stored point
#                          lastPoint[deltaPos], # Starting at last stored point
#                          lastPoint[lambdaPos], # Starting at last stored point
#                          paste0(so("intensityLink"), "_", # intenisty Link + prior
#                                 methods::slot(
#                                   methods::slot(so("prior"), "beta"), "family")),
#                          paste0(so("observabilityLink"),"_", # observability Link + prior
#                                 methods::slot(
#                                   methods::slot(so("prior"), "delta"), "family")),
#                          methods::slot( # lambdaStar prior
#                            methods::slot(
#                              so("prior"),"lambdaStar"), "family"),
#                          retrievePars(methods::slot( # beta prior parameters
#                            so("prior"),"beta")),
#                          retrievePars(methods::slot( # delta prior parameters
#                            so("prior"),"delta")),
#                          retrievePars( # lambdaStar prior parameters
#                            methods::slot(
#                              so("prior"),"lambdaStar")),
#                          ifelse(is.integer(background), # background class
#                                 "int_mat", "num_mat"),
#                          background, # background data
#                          s("area"), # region area
#                          ifelse(is.integer(so("po")),
#                                 "int_mat", "num_mat"), # po data class
#                          so("po"), # po data
#                          backConfig$bIS - 1, # background intensity columns
#                          backConfig$bOS - 1, # background observability colmns
#                          so("intensitySelection") - 1, # po intensity columns
#                          so("observabilitySelection") - 1, # po obserability columns
#                          0, # MCMC burn-in
#                          mcmc_setup$thin, # MCMC thin
#                          mcmc_setup$iter, # MCMC iterations
#                          cores, verbose)
#                        heatMap <- heatMap + temp$xPrimePred
#                        mcmcRun[[c]] <- do.call(cbind,
#
#                        )
#                        colnames(mcmcRun[[c]]) <- s("parnames")
#                        mcmcRun[[c]] <- coda::mcmc(mcmcRun[[c]], thin = mcmc_setup$thin)
#                        if (chains > 1) cat("Finished chain ",c,".\n\n",sep="")
#                      }
#                      if (chains > 1 && verbose) cat("Total computation time: ", format(unclass(Sys.time()-time),
#                                                                                        digits = 2), " ",
#                                                     attr(Sys.time() - time, "units"), ".\n", sep="")
#
#                      return(methods::new("pompp_fit",
#                                          fit = coda::mcmc.list(
#                                            lapply(1:chains, function(i)
#                                              coda::mcmc(rbind(s("fit")[[i]],
#                                                               mcmcRun[[i]])))
#                                          ),
#                                          heatMap = heatMap,
#                                          original = s("original"),
#                                          backgroundSummary = s("backgroundSummary"),
#                                          area = s("area"),
#                                          parnames = s("parnames"),
#                                          mcmc_setup = list(
#                                            burnin = s("mcmc_setup")$burnin,
#                                            thin = mcmc_setup$thin,
#                                            iter = s("mcmc_setup")$iter + mcmc_setup$iter
#                                          )))
#                    } )

# Auxiliary function
checkFormatBackground <- function(object, background){
  # Helper function
  s <- function(n) methods::slot(object, n)

  # Intensity
  if (length(s("iSelectedColumns")) > 0){
    backIntensitySelection <- c()
    for (col in s("iSelectedColumns")){
      if (col %in% colnames(background))
        backIntensitySelection <- c(backIntensitySelection, which(col == colnames(background)))
      else stop(paste0("Column ",col," not found in the background covariates"))
    }
  } else
    if (max(s("intensitySelection")) > ncol(background))
      stop(paste0("Requested background column ", max(s(object,"intensitySelection")),
                  " not available in background which only has ",ncol(background)," columns."))
  else backIntensitySelection = s("intensitySelection")

  # Observability
  if (length(s("oSelectedColumns")) > 0){
    backObservabilitySelection <- c()
    for (col in s("oSelectedColumns")){
      if (col %in% colnames(background))
        backObservabilitySelection <- c(backObservabilitySelection, which(col == colnames(background)))
      else stop(paste0("Column ",col," not found in the background covariates"))
    }
  } else
    if (max(s("observabilitySelection")) > ncol(background))
      stop(paste0("Requested background column ", max(s("observabilitySelection")),
                  " not available in background which only has ",ncol(background)," columns."))
  else
    backObservabilitySelection <- s("observabilitySelection")

  return(list(bIS = backIntensitySelection, bOS = backObservabilitySelection))
}

#' Build a model to be used in the \code{pompp} fitting function
#'
#' Constructor for \code{pompp_model-class} objects, built to facilitate
#' the use of the fitting function. The output of this function has the
#' necessary signature for the fit_pompp function to start the model fit.
#' @param po A matrix whose rows represent the presence-only data and the
#' columns the covariates observed at each position.
#' @param intensitySelection Either a numeric or character vector and
#' represents the selection of covariates used for the intensity set. If
#' numeric it is the positions of the columns and if character, the names of
#' the columns.
#' @param observabilitySelection Either a numeric or character vector and
#' represents the selection of covariates used for the observability set. If
#' numeric it is the positions of the columns and if character, the names of
#' the columns.
#' @param marksSelection Either a numeric or character vector and
#' represents the selection of column used for the marks. If
#' numeric it is the position of the column and if character, the name of
#' the column.
#' @param coordinates Either a numeric or character vector and
#' represents the columns to be used for the coordinates. If
#' numeric it is the positions of the columns and if character, the names of
#' the columns. They must be in longitude, latitude order.
#' @param intensityLink A string to inform what link function the model has
#' with respect to the intensity covariates. Current version accepts 'logit'.
#' @param observabilityLink A string to inform what link function the model has
#' with respect to the observabilitycovariates. Current version accepts 'logit'.
#' @param initial_values Either a single integer, a single
#' \code{pompp_initial-class} or a list containing
#' \code{pompp_initial-class} objects. The length of the list will inform the
#' model how many independent chains will be run. If an integer, that many
#' initial values will be randomly generated.
#' @param joint_prior A \code{pompp_prior} object.
#' @param verbose Set to \code{FALSE} to suppress all messages to console.
#' @return A \code{pompp_model} object with the requested slots. It is ready
#' to be used in the \code{fit_pompp} function.
#' @seealso \code{\link{initial}}, \code{\link{prior}} and
#' \code{\link{fit_pompp}}.
#' @examples
#' # Let us simulate some data to showcase the creation of the model.
#' beta <- c(-1, 2)
#' delta <- c(3, 4)
#' lambdaStar <- 1000
#' gamma <- 2
#' mu <- 5
#'
#' total_points <- rpois(1, lambdaStar)
#' random_points <- cbind(runif(total_points), runif(total_points))
#'
#' # Find covariate values to explain the species occurrence.
#' # We give them a Gaussian spatial structure.
#' Z <- MASS::mvrnorm(1, rep(0, total_points), 3 * exp(-as.matrix(dist(random_points)) / 0.2))
#'
#' # Thin the points by comparing the retaining probabilities with uniforms
#' # in the log scale to find the occurrences
#' occurrences <- log(runif(total_points)) <= -log1p(exp(-beta[1] - beta[2] * Z))
#' n_occurrences <- sum(occurrences)
#' occurrences_points <- random_points[occurrences,]
#' occurrences_Z <- Z[occurrences]
#'
#' # Find covariate values to explain the observation bias.
#' # Additionally create a regular grid to plot the covariate later.
#' W <- MASS::mvrnorm(1, rep(0, n_occurrences), 2 * exp(-as.matrix(dist(occurrences_points)) / 0.3))
#' S <- MASS::mvrnorm(1, rep(0, n_occurrences), 2 * exp(-as.matrix(dist(occurrences_points)) / 0.3))
#'
#' # Find the presence-only observations.
#' marks <- exp(mu + S + rnorm(n_occurrences, 0, 0.3))
#' po_sightings <- log(runif(n_occurrences)) <= -log1p(exp(-delta[1] - delta[2] * W - gamma * S))
#' n_po <- sum(po_sightings)
#' po_points <- occurrences_points[po_sightings, ]
#' po_Z <- occurrences_Z[po_sightings]
#' po_W <- W[po_sightings]
#' po_marks <- marks[po_sightings]
#'
#' # Now we create the model
#' model <- pompp_model(po = cbind(po_Z, po_W, po_points, po_marks),
#'   intensitySelection = 1, observabilitySelection = 2,
#'   marksSelection = 5, coordinates = 3:4,
#'   intensityLink = "logit", observabilityLink = "logit",
#'   initial_values = 2, joint_prior = prior(
#'     NormalPrior(rep(0, 2), 10 * diag(2)),
#'     NormalPrior(rep(0, 3), 10 * diag(3)),
#'     GammaPrior(1e-4, 1e-4),
#'     NormalPrior(0, 100), GammaPrior(0.001, 0.001)))
#' # Check how it is.
#' model
#' @importFrom methods is show
#' @export
pompp_model = function(po, intensitySelection,
                         observabilitySelection,
                         marksSelection, coordinates,
                         intensityLink = "logit", observabilityLink = "logit",
                         initial_values = 1,
                         joint_prior = prior(
                           beta = NormalPrior(
                             rep(0, length(intensitySelection) + 1),
                             10 * diag(length(intensitySelection) + 1)),
                           delta = NormalPrior(
                             rep(0, length(observabilitySelection) + 1),
                             10 * diag(length(observabilitySelection) + 1)),
                           lambdaStar = GammaPrior(
                             1e-10, 1e-10
                           ),
                           marksMean = NormalPrior(1, 100),
                           marksPrecision = GammaPrior(0.001, 0.001)),
                       verbose = TRUE){
  #if (missing(po)) {return(model_builder())} # Easy and interactive model building

  stopifnot(is.matrix(po),
            is.numeric(intensitySelection) || is.character(intensitySelection),
            is.character(intensityLink), length(intensityLink) == 1,
            is.character(observabilityLink), length(observabilityLink) == 1,
            (is.list(initial_values) && all(do.call(c, lapply(initial_values, is, "pompp_initial"))))
            || (is.numeric(initial_values) && length(initial_values) == 1))

  if (is.numeric(intensitySelection)) icharSel = character() else
    if (is.character(intensitySelection)){ # Create columns positions for intensity
      icharSel = intensitySelection
      intensitySelection = c()
      for (selected in icharSel){
        if (!(selected %in% colnames(po))){
          message(paste0("Column ",selected," not found in the observed points covariates. Ignoring it."))
          next
        }
        intensitySelection = c(intensitySelection,which(selected == colnames(po)))
      }
    }
  if (is.numeric(observabilitySelection)) ocharSel = character() else
    if (is.character(observabilitySelection)){ # Create columns positions for intensity
      ocharSel = observabilitySelection
      observabilitySelection = c()
      for (selected in ocharSel){ # Putting indexes together
        if (!(selected %in% colnames(po))){
          message(paste0("Column ",selected," not found in the observed points covariates. Ignoring it."))
          next
        }
        observabilitySelection = c(observabilitySelection,which(selected==colnames(po)))
      }
    }
  if (is.numeric(marksSelection)) marksSel = character() else
    if (is.character(marksSelection)) {
      marksSel = marksSelection
      if (!(marksSel %in% colnames(po))) stop(paste0("Column",
                                                     marksSel,
                                                     "not found in the observed points covariates."))
      marksSelection = which(marksSel == colnames(po))
    }
  if (is.numeric(coordinates)) coordCharSel = character() else
    if (is.character(coordinates)){ # Create columns positions for intensity
      coordCharSel = coordinates
      coordinates = c()
      for (selected in coordCharSel){ # Putting indexes together
        if (!(selected %in% colnames(po))){
          stop(paste0("Column ",selected," not found in the observed points covariates."))
        }
        coordinates = c(coordinates, which(selected==colnames(po)))
      }
    }
  if (methods::is(initial_values,"pompp_initial"))
    initial_values = list(initial_values)
  if (is.numeric(initial_values))
    initial_values = initial(length(intensitySelection) + 1,
                             length(observabilitySelection) + 2,
                             nrow(po), 0, 1, random=TRUE) * initial_values

  if (verbose) cat("Loading data with", nrow(po), "observed points.\n")
  output <- methods::new("pompp_model", po = po, intensityLink = intensityLink,
                         intensitySelection = intensitySelection,
                         observabilityLink = observabilityLink,
                         observabilitySelection = observabilitySelection,
                         marksSelection = marksSelection,
                         coordinates = coordinates,
                         init = initial_values, prior = joint_prior,
                         iSelectedColumns = icharSel, oSelectedColumns = ocharSel,
                         mSelectedColumns = marksSel, cSelectedColumns = coordCharSel)
  if (verbose) {
    cat("Data loaded successfully with ", length(intensitySelection), " intensity variables and ",
        length(observabilitySelection), " observability variables selected.\n", length(initial_values), " chains ",
        ifelse(length(initial_values) > 1, "were", "was"), " initialized.\n", sep = "")
    if (!length(icharSel)) cat("Intensity covariates selected with column indexes. Make sure the background covariates are in the same position.\n")
    if (!length(ocharSel)) cat("Observability covariates selected with column indexes. Make sure the background covariates are in the same position.\n")
    if (!length(coordCharSel)) cat("Coordinates selected with column indexes. Make sure the background covariates are in the same position.\n")
  }

  return(output)
}
