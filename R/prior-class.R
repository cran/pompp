#### Individual priors ####
methods::setGeneric("retrievePars", function(object){standardGeneric("retrievePars")})

#### Beta-Delta priors ####
#' Generic class for the beta and delta parameters.
#'
#' @field family The family of distributions of the prior.
#' @export
#' @exportClass BetaDeltaPrior
methods::setClass("BetaDeltaPrior", representation(family="character"))

#' @rdname BetaDeltaPrior-class
#' @param object The BetaDeltaPrior object.
#' @return \strong{\code{show}} and \strong{\code{print}}: The invisible object.
#' @export
#' @exportMethod show
methods::setMethod("show", "BetaDeltaPrior", function(object){
  cat("Prior for Beta and Delta from the", methods::slot(object, "family"), "family.\n")
  invisible(object)
})

#' @rdname BetaDeltaPrior-class
#' @param x The BetaDeltaPrior object.
#' @param ... Ignored.
#' @export
#' @exportMethod print
methods::setMethod("print", "BetaDeltaPrior", function(x, ...) methods::show(x))

#' @rdname BetaDeltaPrior-class
#' @method print BetaDeltaPrior
#' @export
print.BetaDeltaPrior <- function(x, ...) methods::show(x)

methods::setMethod("retrievePars", "BetaDeltaPrior", function(object){
  stop("Please specify a specific prior distribution for Beta and Delta.")
})

#' Normal prior class for Beta and Delta parameters.
#'
#' This is used to represent the prior for Beta and Delta individually. They
#' still need to be joined to be used in a model.
#' @field mu The mean vector for the prior.
#' @field Sigma The covariance matrix for the prior.
#' @export
#' @exportClass NormalPrior
methods::setClass("NormalPrior", contains="BetaDeltaPrior",
                  representation = methods::representation(mu = "numeric", Sigma = "matrix"),
                  validity = function(object){
                    if (!isSymmetric(methods::slot(object, "Sigma"))) stop("Covariance matrix is not symmetric.")
                    if (any(eigen(methods::slot(object, "Sigma"), symmetric=TRUE)$values <= 0)) stop("Covariance matrix is not positive-definite.")
                    if (length(methods::slot(object, "mu")) != nrow(methods::slot(object, "Sigma"))) stop(paste0("Prior mean vector and covariance matrix have incompatible sizes. Mean length is ", length(methods::slot(object, "mu")), " and covariance matrix is ", nrow(methods::slot(object, "Sigma")), " x ", ncol(methods::slot(object, "Sigma")), "."))
                    TRUE
                  })

#' Create a Normal prior object for model specification.
#'
#' Constructor for \code{NormalPrior-class} objects
#' @param mu The mean vector for the Normal distribution.
#' @param Sigma The covariance matrix for the Normal distribution.
#' @details Matrix Sigma must be square and positive definite. Its dimensions
#' must match mu's length.
#' @return A \code{NormalPrior} object with adequate slots.
#' @seealso \code{\link{prior}}
#' @examples
#' NormalPrior(rep(0, 10), diag(10) * 10)
#' @export
NormalPrior <- function(mu, Sigma){
  mu <- as.numeric(mu)
  Sigma <- as.matrix(Sigma)
  methods::new("NormalPrior", mu = mu, Sigma = Sigma, family = "normal")
}

## Methods
#' @rdname NormalPrior-class
#' @param x The NormalPrior object.
#' @return \strong{\code{names}}: A character vector with the prior parameters.
#' @export
#' @exportMethod names
methods::setMethod("names","NormalPrior",function(x) c("mu","Sigma"))

#' @rdname NormalPrior-class
#' @param x The NormalPrior object.
#' @param name The requested slot.
#' @return \strong{\code{`$`}}: The requested slot's value.
#' @export
#' @exportMethod $
methods::setMethod("$","NormalPrior",function(x,name) methods::slot(x, name))

#' @rdname NormalPrior-class
#' @param x The NormalPrior object.
#' @param name The requested slot.
#' @param value New value.
#' @return \strong{\code{`$<-`}}: The new object with the updated slot.
#' @export
#' @exportMethod $<-
methods::setMethod("$<-","NormalPrior",function(x,name,value){
  methods::slot(x, name) <- value
  methods::validObject(x)
  x
})

#' @rdname NormalPrior-class
#' @param object The NormalPrior object.
#' @return \strong{\code{show}} and \strong{\code{print}}: The invisible object.
#' @export
#' @exportMethod show
methods::setMethod("show","NormalPrior",function(object){
  cat("Normal prior\n\nMu mean vector:\n")
  print(methods::slot(object, "mu"))
  s = methods::slot(object, "Sigma")
  cat("Sigma covariance matrix:\n")
  if (length(s) == 1)
    print(s)
  else if (max(s - diag(diag(s)))) # If Sigma is not a diagonal matrix
    print(methods::slot(object, "Sigma"))
  else if (min(diag(s)) != max(diag(s))) # If Sigma is a diagonal matrix of different values
    cat("diag(c(", paste(diag(s), collapse = ","), "))\n", sep="")
  else
    cat(s[1], "* I, where I is an identity matrix.\n")
  invisible(object)
})

#' @rdname NormalPrior-class
#' @param x The NormalPrior object.
#' @param ... Ignored.
#' @export
#' @exportMethod print
methods::setMethod("print","NormalPrior",function(x,...) methods::show(x))

#' @rdname NormalPrior-class
#' @param x The NormalPrior object.
#' @param ... Ignored.
#' @method print NormalPrior
#' @export
print.NormalPrior <- function(x, ...) methods::show(x)

methods::setMethod("retrievePars", "NormalPrior", function(object){
  list(mean = methods::slot(object, "mu"),
       covariance = methods::slot(object, "Sigma"))
})

#### LambdaStar priors ####
#' Generic class for the LambdaStar parameters.
#'
#' @field family The family of distributions of the prior.
#' @export
#' @exportClass LambdaStarPrior
methods::setClass("LambdaStarPrior", methods::representation(family = "character"))

#' @rdname LambdaStarPrior-class
#' @param object The LambdaStarPrior object.
#' @return \strong{\code{show}} and \strong{\code{print}}: The invisible object.
#' @export
#' @exportMethod show
methods::setMethod("show", "LambdaStarPrior", function(object){
  cat("Prior for LambdaStar from the", methods::slot(object, "family"), "family.\n")
  invisible(object)
})

#' Gamma prior class for the LambdaStar parameter.
#'
#' This is used to represent the prior for lambdaStar individually. It
#' still needs to be joined with the prior for Beta and Delta to be used
#' in a model.
#' @field shape The shape parameter of the Gamma distribution.
#' @field rate The rate parameter of the Gamma distribution.
#' @seealso \code{\link{prior}}
#' @examples
#' GammaPrior(0.0001, 0.0001)
#' @export
#' @exportClass GammaPrior
methods::setClass("GammaPrior", contains="LambdaStarPrior",
                  representation = methods::representation(shape = "numeric", rate = "numeric"),
                  validity = function(object){
                    for (par in c("shape", "rate")){
                      if (length(methods::slot(object, par)) > 1) stop(paste0("Prior parameter ", par, " for lambdaStar must have length 1."))
                      if (methods::slot(object, par) <= 0) stop("Prior parameters must be positive")
                    }
                    TRUE
                  })

methods::setMethod("retrievePars", "GammaPrior", function(object){
  stop("Please specify a specific prior distribution for LambdaStar.")
})

## Constructor
#' Create a Gamma prior object for model specification.
#'
#' Constructor for \code{GammaPrior-class} objects
#' @param shape A positive number.
#' @param rate A positive number.
#' @return A \code{GammaPrior} object with adequate slots.
#' @export
GammaPrior <- function(shape, rate) {
  stopifnot(shape > 0, rate > 0)
  new("GammaPrior",shape = shape, rate = rate, family = "gamma")
}

## Methods
#' @rdname GammaPrior-class
#' @param x The GammaPrior object.
#' @return \strong{\code{names}}: A character vector with the prior parameters.
#' @export
#' @exportMethod names
methods::setMethod("names", "GammaPrior", function(x) c("shape", "rate"))

#' @rdname GammaPrior-class
#' @param x The Gamma<- object.
#' @param name The requested slot.
#' @return \strong{\code{`$`}} The requested slot's value.
#' @export
#' @exportMethod $
methods::setMethod("$", "GammaPrior", function(x, name) methods::slot(x, name))

#' @rdname GammaPrior-class
#' @param x The GammaPrior object.
#' @param name The requested slot.
#' @param value New value.
#' @return \strong{\code{`$<-`}}: The new object with the updated slot.
#' @export
#' @exportMethod $<-
methods::setMethod("$<-", "GammaPrior", function(x, name, value){
  methods::slot(x, name) <- value
  methods::validObject(x)
  x
})

#' @rdname GammaPrior-class
#' @param object The GammaPrior object.
#' @return \strong{\code{show}} and \strong{\code{print}}: The invisible object.
#' @export
#' @exportMethod show
methods::setMethod("show", "GammaPrior", function(object){
  cat("Gamma prior\n\n")
  vec <- c(methods::slot(object, "shape"), methods::slot(object, "rate"))
  names(vec) <- names(object)
  print(vec)
  invisible(object)
})

#' @rdname GammaPrior-class
#' @param x The GammaPrior object.
#' @param ... Ignored.
#' @export
#' @exportMethod print
methods::setMethod("print", "GammaPrior", function(x, ...) methods::show(x))

#' @rdname GammaPrior-class
#' @param x The GammaPrior object.
#' @param ... Ignored.
#' @method print GammaPrior
#' @export
print.GammaPrior <- function(x, ...) methods::show(x)

methods::setMethod("retrievePars", "GammaPrior", function(object){
  list(a = methods::slot(object, "shape"),
       b = methods::slot(object, "rate"))
})

#### Joint prior ####
#' Joint prior class for the pompp package parameters
#'
#' Objects of this class are the joining of independent priors for Beta, Delta
#' and LambdaStar. They can be used in the \code{fit_pompp} function.
#' @field beta An object of a class which inherits the \code{BetaDeltaPrior} S4
#' class with the appropriate Beta prior.
#' @field delta An object of a class which inherits the \code{BetaDeltaPrior} S4
#' class with the appropriate Delta prior.
#' @field lambdaStar An object of a class which inherits the
#' \code{LambdaStarPrior} S4 class with the appropriate LambdaStar prior.
#' @field marksMean An object of S4 class \code{NormalPrior} with the chosen
#' prior for the marks mean
#' @field marksPrecision An object of S4 class \code{GammaPrior} with the chosen
#' prior for the marks precision
#' @export
#' @exportClass pompp_prior
methods::setClass("pompp_prior",
                  methods::representation(
                    beta = "BetaDeltaPrior",
                    delta = "BetaDeltaPrior",
                    lambdaStar = "LambdaStarPrior",
                    marksMean = "NormalPrior",
                    marksPrecision = "GammaPrior"))

#' @rdname pompp_prior-class
#' @param x The pompp_prior object.
#' @export
#' @exportMethod names
methods::setMethod("names", "pompp_prior", function(x) c("beta", "delta", "lambdaStar", "marksMean", "marksPrecision"))

#' @rdname pompp_prior-class
#' @param x The pompp_prior object.
#' @param name The requested slot.
#' @export
#' @exportMethod $
methods::setMethod("$", "pompp_prior", function(x,name) methods::slot(x, name))

#' Build a joint prior for pompp model parameters
#'
#' Constructor for \code{pompp_prior} objects, which is used in the
#' \code{pompp_fit} function. The generated prior is so that Beta, Delta
#' and LambdaStar are indepdendent a priori.
#' @param beta An S4 object whose class inherits from \code{BetaDeltaPrior}.
#' @param delta An S4 object whose class inherits from \code{BetaDeltaPrior}.
#' @param lambdaStar An S4 object whose class inherits from \code{LambdaStarPrior}.
#' @param marksMean An S4 object of class \code{NormalPrior}.
#' @param marksPrecision An S4 object of class \code{GammaPrior}.
#' @return A \code{pompp_prior} object with the adequate slots. It is ready to
#' be included in a model via the \code{pompp_model} function.
#' @seealso \code{\link{fit_pompp}}, \code{\link{NormalPrior}},
#' \code{\link{GammaPrior}} and \code{\link{pompp_model}}.
#' @examples
#' # Let us say there are 3 intensity covariates and 4 observability covariates.
#' # One more element is included in both sets due to the intercepts.
#' new_prior <- prior(
#'   NormalPrior(rep(0, 4), 10 * diag(4)),
#'   NormalPrior(rep(0, 5), 10 * diag(5)),
#'   GammaPrior(0.0001, 0.0001),
#'   NormalPrior(0, 100), GammaPrior(0.001, 0.001)
#' )
#' @export
prior <- function(beta, delta, lambdaStar, marksMean, marksPrecision){
  if (length(marksMean$mu) > 1) stop("marksMeans must have a one-dimensional Normal prior.")
  methods::new("pompp_prior", beta = beta, delta = delta, lambdaStar = lambdaStar,
               marksMean = marksMean, marksPrecision = marksPrecision)
}

#' @rdname pompp_prior-class
#' @param object The pompp_prior object.
#' @return \strong{\code{names}}: A character vector with the model parameters
#' names.
#' @export
#' @exportMethod show
setMethod("show", "pompp_prior", function(object){
  cat("Joint prior for a pompp model. Individual components:\n\nBeta:\n")
  show(methods::slot(object, "beta"))
  cat("\nDelta:\n")
  show(methods::slot(object, "delta"))
  cat("\nlambdaStar:\n")
  print(methods::slot(object, "lambdaStar"))
  cat("\nMarks mean:\n")
  print(methods::slot(object, "marksMean"))
  cat("\nMarks precision:\n")
  print(methods::slot(object, "marksPrecision"))
})

#' @rdname pompp_prior-class
#' @param x The pompp_prior object.
#' @param ... Ignored.
#' @export
#' @exportMethod print
setMethod("print", "pompp_prior", function(x, ...) methods::show(x))

#' @rdname pompp_prior-class
#' @param x The pompp_prior object.
#' @param ... Ignored.
#' @method print pompp_prior
#' @export
print.pompp_prior <- function(x, ...) methods::show(x)

#' @rdname pompp_prior-class
#' @param x The pompp_prior object.
#' @param name The requested slot.
#' @return \strong{\code{`$`}}: The requested slot's value.
#' @export
#' @exportMethod $
setMethod("$", "pompp_prior", function(x, name) methods::slot(x, name))

#' @rdname pompp_prior-class
#' @param x The pompp_prior object.
#' @param name The requested slot.
#' @param value New value.
#' @return \strong{\code{`$<-`}}: The new object with the updated slot.
#' @export
#' @exportMethod $<-
setMethod("$<-", "pompp_prior", function(x, name, value){
  methods::slot(x, name) <- value
  methods::validObject(x)
  x
})
