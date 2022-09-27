methods::setOldClass("mcmc.list")
methods::setOldClass("table")
methods::setOldClass("list")

#' @include model-class.R
NULL

#' Class for the result of the MCMC procedure.
#'
#' Objects of this class are the main objects of this package. They contain
#' much information about the fitted model.
#' @field fit The actual fit from the model. It is an object of class
#' \code{\link[coda]{mcmc.list}}, as generated from the \code{coda} package.
#' @field original The model used to generate the chains, an object with class
#' \code{pompp_model}.
#' @field backgroundSummary A small summary of the original background
#' covariates. This is to ensure that continuing the chains will use the
#' identical background matrix. Only the summary is kept for storage efficiency.
#' @field area A positive number indicating the area measure of the region being
#' studied.
#' @field parnames The names of the parameters. If the model used selects the
#' covariates with column names, they are replicated here. If they are the
#' column indexes, names are generated for identification.
#' @field mcmc_setup The original mcmc setup used.
#' @seealso \code{\link{fit_pompp}}
#' @export
#' @exportClass pompp_fit
methods::setClass("pompp_fit",
                  representation(fit = "mcmc.list",
                                 original = "pompp_model",
                                 neighborhoodSize = "numeric",
                                 backgroundSummary = "table",
                                 area = "numeric",
                                 parnames = "character",
                                 mcmc_setup = "list"))

#### Basic methods ####
#' @rdname pompp_fit-class
#'
#' @param object A pompp_fit object.
#' @return \strong{\code{show}} and \strong{\code{print}}: The invisible object.
#' @export
#' @exportMethod show
methods::setMethod("show","pompp_fit",function(object){
  cat("Fit of a pompp model.\n")

  ## data
  cat("Presence-only dataset size:",nrow(methods::slot(methods::slot(object,"original"),"po")),"\n\n")
  sc <- methods::slot(methods::slot(object,"original"), "iSelectedColumns")
  if (length(sc))
    cat(length(sc), " intensity covariates selected:\n", sc, "\n")
  else{
    sc <- methods::slot(methods::slot(object,"original"),"intensitySelection")
    cat(length(sc), " intensity covariates selected. Columns: ", sc, "\n")
  }
  sc <- methods::slot(methods::slot(object,"original"),"oSelectedColumns")
  if (length(sc))
    cat(length(sc), " observability covariates selected:\n", sc, "\n")
  else{
    sc <- methods::slot(methods::slot(object,"original"),"observabilitySelection")
    cat(length(sc), " observability covariates selected. Columns: ", sc, "\n")
  }
  sc <- methods::slot(methods::slot(object, "original"), "marksSelection")
  if (length(sc)) cat("Marks column: ", sc, "\n")
  cat("\n")

  ## Link function
  links <- c(methods::slot(methods::slot(object,"original"), "intensityLink"),
             methods::slot(methods::slot(object,"original"), "observabilityLink"))
  names(links) <- c("intensity", "observability")
  cat("Link functions chosen:\n")
  print(links)
  cat("\n")

  ## Prior
  cat("Prior selection:\n")
  methods::show(methods::slot(methods::slot(object, "original"), "prior"))
  cat("\n")

  ## MCMC configuration
  chains = length(methods::slot(methods::slot(object, "original"), "init"))
  setup = methods::slot(object, "mcmc_setup")
  cat(chains,ifelse(chains > 1," chains"," chain"), " of MCMC ",
      ifelse(chains > 1,"were", "was"), " configured with ",
      format(setup$burnin, scientific = FALSE), " warmup ",
      ifelse(setup$burnin > 1, "iterations", "iteration"),
      " and ", format(setup$iter, scientific = FALSE), " valid ",
      ifelse(setup$iter>1, "iterations","iteration"), ", storing one in every ",
      ifelse(setup$thin>1, paste(setup$thin, "steps"), "step"), ".\n\n", sep="")

  ## Results
  print(round(summary(object), digits = 3))
  cat("\n")

  ## Comments
  cat("The effective sample size represents the sample size of an independent",
      "sample which would yield equivalent estimates as the autocorrelated",
      "MCMC result.\n")
  if (chains > 1)
    cat("Rhat has been calculated from multiple chains. Lower, closer to 1",
        "values indicate better convergence of the chain. For posterior",
        "estimates to be trusted, the upper CI limit should be below 1.1. If",
        "they are not, run more iterations. See help('fit_pompp') to see",
        "how to do utilize the iterations already run.\n\n")
  else
    cat("Rhat cannot be estimated with only 1 chain. Run more chains for this",
        "statistic to be displayed.")

  invisible(object)
})

#' @rdname pompp_fit-class
#'
#' @param x A pompp_fit object.
#' @param ... Ignored.
#' @export
#' @exportMethod print
methods::setMethod("print", "pompp_fit", function(x, ...) methods::show(x))

#' @method print pompp_fit
#' @export
#' @rdname pompp_fit-class
print.pompp_fit <- function(x,...) methods::show(x)

#' @rdname pompp_fit-class
#'
#' @param object A pompp_fit object.
#' @param ... Ignored.
#' @export
#' @exportMethod summary
methods::setMethod("summary", "pompp_fit", function(object,...) summary.pompp_fit(object, ...))

#' @method summary pompp_fit
#' @return \strong{\code{summary}}: A matrix with the summary statistics of the
#' fit. It is also printed in the \code{print} method. The summary can be
#' treated as a matrix, such as retrieving rows/columns and creating tables
#' with the \code{xtable} package.
#' @rdname pompp_fit-class
#' @export
summary.pompp_fit <- function(object, ...){
  chains <- length(methods::slot(methods::slot(object, "original"), "init"))
  nb <- length(methods::slot(methods::slot(object, "original"), "intensitySelection")) + 1
  nd <- length(methods::slot(methods::slot(object, "original"), "observabilitySelection")) + 2
  npar <- nb + nd + 7 # +1 = lambdaStar + marksMu + marksTau2 + nU + nX' + marksSum + marksVariance
  result <- matrix(0, npar, 6) # Mean, median, sd, lower CI bound, upper CI bound, effective sample size
  colnames(result) <- c("mean", "median", "sd", "2.5%", "97.5%", "eff. sample size")
  rownames(result) <- methods::slot(object, "parnames")[1:npar]
  fitToMatrix <- as.matrix(methods::slot(object, "fit"))[, 1:npar]
  result[,1] <- colMeans(fitToMatrix)
  result[,2] <- apply(fitToMatrix, 2, stats::median)
  result[,3] <- apply(fitToMatrix, 2, stats::sd)
  result[,4] <- apply(fitToMatrix, 2, stats::quantile, 0.025)
  result[,5] <- apply(fitToMatrix, 2, stats::quantile, 0.975)
  result[,6] <- coda::effectiveSize(methods::slot(object,"fit"))[1:npar]
  if (chains > 1){
    cols <- colnames(result)
    result <- cbind(result, coda::gelman.diag(methods::slot(object, "fit"))$psrf[1:npar, ])
    colnames(result) <- c(cols, "Estimated Rhat","Upper CI Rhat")
  }
  result
}

#' @rdname pompp_fit-class
#'
#' @param x A pompp_fit object.
#' @return \strong{\code{names}}: A character vector with the available options
#' for the \code{`$`} and \code{`[[`} methods.
#' @export
#' @exportMethod names
methods::setMethod("names", "pompp_fit", function(x){
  nn <- c("parameters", "covariates_importance", "mcmc_chains", "model",
          "eff_sample_size", "area", "initial_values", "mcmc_setup")
  if (length(methods::slot(methods::slot(x, "original"), "init")) > 1)
    nn <- c(nn, "Rhat", "Rhat_upper_CI")
  nn
})

#' @method names pompp_fit
#' @rdname pompp_fit-class
#' @export
names.pompp_fit <- function(x) names(x)

#' @rdname pompp_fit-class
#'
#' @param x A pompp_fit object.
#' @param i The requested slot.
#' @return \strong{\code{`$`}} and \strong{\code{`[[`}}: The requested slot.
#' Available options are not necessarily the class slots, and can be checked
#' with the \code{names} method.
#' @export
#' @exportMethod [[
methods::setMethod("[[", "pompp_fit", function(x, i){
  # Helper function
  s <- function(n) methods::slot(x, n)

  nb <- length(methods::slot(s("original"), "intensitySelection")) + 1
  nd <- length(methods::slot(s("original"), "observabilitySelection")) + 2
  npar <- nb + nd + 7 # +1 from lambdaStar + marksMu + marksTau2 + nU + nX' + marksSum + marksVariance

  if (i == "parameters"){
    summ <- summary(x)
    output <- summ[, 1]
    names(output) <- rownames(summ)
  } else
    if (i == "covariates_importance"){
      data <- as.data.frame(x)
      names(data) <- namesAid(names(data))
      obsInterceptName <- names(data)[nb + 1]

      intensity <- as.matrix(t(rbind(rep(1, nrow(data)), apply(
        data[2:(which(names(data) == obsInterceptName) - 1)], 1,
        function(chain) {c2 <- chain * chain; c2 / sum(c2)}
      )))[, -1])
      observability <- as.matrix(t(rbind(rep(1, nrow(data)), apply(
        data[(which(names(data) == obsInterceptName) + 1):(which(names(data) == "lambdaStar") - 1)], 1,
        function(chain) {c2 <- chain * chain; c2 / sum(c2)}
      )))[, -1])
      colnames(intensity) <- names(data)[2:(which(names(data) == obsInterceptName) - 1)]
      colnames(observability) <- names(data)[(which(names(data) == obsInterceptName) + 1):(which(names(data) == "lambdaStar") - 1)]
      output <- list(intensity = intensity, observability = observability)
      class(output) <- "covariates_importance"
    } else
      if (i == "eff_sample_size"){
        summ <- summary(x)
        output <- summ[, 6]
        names(output) <- rownames(summ)
      } else
        if (i == "Rhat"){
          summ <- summary(x)
          output <- summ[, 7]
          names(output) <- rownames(summ)
        } else
          if (i == "Rhat_upper_CI"){
            summ <- summary(x)
            output <- summ[, 8]
            names(output) <- rownames(summ)
          } else
            if (i == "mcmc_chains") output <- s("fit") else
              if (i == "model") output <- s("original") else
                if (i == "initial_values") output <- methods::slot(s("original"), "init") else
                  if (i == "mcmc_setup") output <- s("mcmc_setup") else
#                    if (i == "log_posterior") output <- as.data.frame(x)$log_Posterior else
                      if (i == "area") output <- s("area")

  output
})

#' @rdname pompp_fit-class
#'
#' @param x A pompp_fit object.
#' @param name The requested slot.
#' @export
#' @exportMethod $
methods::setMethod("$", "pompp_fit", function(x, name) x[[name]])

#' @rdname pompp_fit-class
#' @param x A pompp_fit object.
#' @param ... Ignored.
#' @export
#' @exportMethod as.array
methods::setMethod("as.array", "pompp_fit", function(x, ...) as.array.pompp_fit(x, ...))

namesAid <- function(string){
  new_string <- string

  intInt <- "(Intensity intercept)"
  obsInt <- "(Observability intercept)"
  obsStart <- max(which(string == obsInt), which(string == "delta_0"))
  obsEnd <- which(string == "lambdaStar") - 1

  # Find same covariates in both sets
  for (i in 1:(obsStart - 1))
    searching <- string[i] == string[obsStart:obsEnd]
  if (any(searching)){
    new_string[i] <- paste0(string[i], ".int")
    new_string[obsStart - 1 + which(searching)] <- paste0(string[i], ".obs")
  }
  new_string[1] <- ifelse(string[1] == intInt, "Intensity_Intercept", "beta_0")
  new_string[obsStart] <- ifelse(string[obsStart] == obsInt, "Observability_Intercept", "delta_0")

  new_string
}

#' @rdname pompp_fit-class
#' @param x A pompp_fit object.
#' @param ... Ignored in this version.
#' @return \strong{\code{as.array}}: An \code{array} with dimensions I x C x P,
#' where I stands for number of iterations, C for number of chains and P for
#' total number of parameters. P is actually larger than the number of
#' parameters in the model, as the the generated sizes of the latent processes
#' and the log-posterior are also included. This is organized so that is ready
#' for the \code{bayesplot} package functions.
#' @method as.array pompp_fit
#' @export
as.array.pompp_fit <- function(x, ...){
  nchains <- length(methods::slot(x, "fit"))
  chains <- do.call(rbind, methods::slot(x, "fit"))
  chains <- chains[, -((ncol(chains) - 1):ncol(chains))]
  iterations <- nrow(methods::slot(x, "fit")[[1]])
  npar <- ncol(chains)

  ## Format to be used with bayesplot:: functions
  parnames <- methods::slot(x, "parnames")
  return(
    array(chains,
          dim = c(iterations, nchains, npar),
          dimnames = list(iterations = NULL,
                          chains = paste0("chain:", 1:nchains),
                          parameters = namesAid(parnames[1:(length(parnames) - 2)])))
  )
}

#' @export
#' @rdname pompp_fit-class
#' @exportMethod as.matrix
methods::setMethod("as.matrix", "pompp_fit", function(x, ...) as.matrix.pompp_fit(x, ...))

#' @rdname pompp_fit-class
#' @param x A pompp_fit object.
#' @param ... Ignored in this version.
#' @return \strong{\code{as.matrix}}: The dimension of the output is
#' I * C x (P + 2), where I stands for number of iterations, C for number of
#' chains and P for total number of parameters. P is actually larger than the
#' number of parameters in the model, as the generated sizes of the latent
#' processes and the log-posterior are also included.
#'
#' Two extra columns are included to indicate to which chain and to which
#' iteration that draw belongs.
#' @method as.matrix pompp_fit
#' @export
as.matrix.pompp_fit <- function(x, ...){
  nchains <- length(methods::slot(x, "fit"))
  chains <- do.call(rbind, methods::slot(x, "fit"))
  chains <- chains[, -((ncol(chains) - 1):ncol(chains))]
  iterations <- nrow(methods::slot(x, "fit")[[1]])
  parnames <- namesAid(colnames(chains))
  chains <- cbind(chains, rep(factor(1:nchains), each = iterations))
  chains <- cbind(chains, rep(1:iterations, nchains))
  colnames(chains) <- c(parnames, "chain", "iteration")

  return(chains)
}

#' @export
#' @rdname pompp_fit-class
#' @exportMethod as.data.frame
methods::setMethod("as.data.frame","pompp_fit",function(x, row.names = NULL, optional = FALSE, ...) as.data.frame.pompp_fit(x, row.names = NULL, optional = FALSE, ...))

#' @rdname pompp_fit-class
#' @param x A pompp_fit object.
#' @param row.names NULL or a character vector giving the row names for the
#' data frame. Missing values are not allowed.
#' @param optional logical. If TRUE, setting row names and converting column
#' names to syntactic names is optional. See help('as.data.frame') for more.
#' Leaving as \code{FALSE} is recommended.
#' @param ... Ignored in this version.
#' @return \strong{\code{as.data.frame}}: The dimension of the output is
#' I*C x P + 2, where I stands for number of iterations, C for number of chains
#' and P for total number of parameters. P is actually larger than the number
#' of parameters in the model, as the generated sizes of the latent processes
#' and the log-posterior are also included.
#'
#' Two extra columns are included to indicate to which chain and to which
#' iteration that draw belongs. This is to facilitate the use of plotting
#' results via the \code{ggplot2} package if desired.
#'
#' If \code{row.names} is left at \code{NULL} then row names are created as
#' CcIi where c is the chain and i is the iteration of that row.
#' @method as.data.frame pompp_fit
#' @export
as.data.frame.pompp_fit = function(x, row.names = NULL, optional = FALSE, ...){
  nchains <- length(methods::slot(x, "fit"))
  chains <- do.call(rbind, methods::slot(x, "fit"))
  chains <- chains[, -((ncol(chains) - 1):ncol(chains))]
  parnames <- namesAid(colnames(chains))
  iterations <- nrow(methods::slot(x, "fit")[[1]])

  colsList <- list()
  for (pp in 1:length(parnames)) colsList[[parnames[pp]]] <- chains[, pp]

  if (!is.null(NULL))
    row_names <- row.names else
      row_names <- paste0(rep(paste0("C", 1:nchains),
                              each=iterations),
                          rep(paste0("I", 1:iterations), nchains))
  output <- do.call(data.frame, c(colsList, list(check.names = TRUE,
                                                 fix.empty.names = TRUE),
                                  list(row.names = row_names)))
  output$chain <- factor(rep(1:nchains, each=iterations))
  output$iteration <- rep(1:iterations, nchains)

  return(output)
}

#### Interaction methods ####
# Adding chains into a single object
#' @rdname pompp_fit-class
#' @param e1 A pompp_fit object.
#' @param e2 A pompp_fit object with the same slots (except for initial
#' values) as \code{e1}.
#' @return \strong{\code{+}}: A new \code{pompp_fit} object where the chains
#' are combined into a new multi-chain object. This can be used if chains are
#' run in separate occasions or computers to combine them into a single object
#' for analysis.
#' @importFrom methods new
methods::setMethod("+", methods::signature(e1 = "pompp_fit", e2 = "pompp_fit"),
                   function(e1, e2){
                     s1 <- function(n) methods::slot(e1, n)
                     so1 <- function(n) methods::slot(s1("original"), n)
                     s2 <- function(n) methods::slot(e2, n)
                     so2 <- function(n) methods::slot(s2("original"), n)
                     stopifnot(all.equal(so1("po"), so2("po")),
                               all.equal(so1("intensityLink"), so2("intensityLink")),
                               all.equal(so1("intensitySelection"), so2("intensitySelection")),
                               all.equal(so1("observabilityLink"), so2("observabilityLink")),
                               all.equal(so1("observabilitySelection"), so2("observabilitySelection")),
                               all.equal(so1("marksSelection"), so2("marksSelection")),
                               all.equal(so1("prior"), so2("prior")),
                               all.equal(so1("iSelectedColumns"), so2("iSelectedColumns")),
                               all.equal(so1("oSelectedColumns"), so2("oSelectedColumns")),
                               all.equal(so1("mSelectedColumns"), so2("mSelectedColumns")),
                               all.equal(s1("backgroundSummary"), s2("backgroundSummary")),
                               all.equal(s1("area"), s2("area")),
                               all.equal(s1("parnames"), s2("parnames")),
                               all.equal(s1("mcmc_setup"), s2("mcmc_setup")))

                     fff <- coda::mcmc.list(c(s1("fit"), s2("fit")))

                     or <- s1("original")
                     methods::slot(or, "init") <- c(methods::slot(s1("original"), "init"),
                                                    methods::slot(s2("original"), "init"))

                     return(methods::new("pompp_fit",
                                         fit = fff,
                                         original = or,
                                         backgroundSummary = s1("backgroundSummary"),
                                         area = s1("area"),
                                         parnames = s1("parnames"),
                                         mcmc_setup = s1("mcmc_setup")))
                   })

# Combining multiple chains
#' @rdname pompp_fit-class
#' @return \strong{\code{c}}: A new \code{pompp_fit} object where the chains
#' are combined into a new multi-chain object. The \strong{\code{+}} method is
#' used for that, with the same arguments restrictions and results.
#' @export
#' @exportMethod c
methods::setMethod("c", "pompp_fit", function(x, ...) {
  ll <- list(...)
  res <- x
  for (i in 1:length(ll))
    res <- res + ll[[i]]

  res
})

