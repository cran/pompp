#' Class for the initial values for the MCMC for the pompp package
#'
#' @field beta Initial values for beta.
#' @field delta Initial values for delta.
#' @field lambdaStar Initial values for lambdaStar.
#' @field tag Indicates the source of the initial values.
#' @export
#' @exportClass pompp_initial
methods::setClass("pompp_initial",
                  methods::representation(beta = "numeric",
                                          delta = "numeric",
                                          lambdaStar = "numeric",
                                          marksMean = "numeric",
                                          marksPrecision = "numeric",
                                          tag = "character"),
                  validity = function(object){
                    if (length(methods::slot(object,"beta")) < 1) stop("Beta must have at least 1 initial value")
                    if (length(methods::slot(object,"delta")) < 1) stop("Delta must have at least 1 initial value")
                    if (length(methods::slot(object,"lambdaStar")) != 1) stop("lambdaStar must have exactly 1 initial value")
                    if (length(methods::slot(object,"marksMean")) != 1) stop("marksMean must have exactly 1 initial value")
                    if (length(methods::slot(object,"marksPrecision")) != 1) stop("marksPrecision must have exactly 1 initial value")
                    TRUE
                  })

#' @rdname pompp_initial-class
#' @param x The pompp_initial object.
#' @return \strong{\code{names}}: A character vector with the initialized
#' parameter names.
#' @export
#' @exportMethod names
methods::setMethod("names","pompp_initial", function(x) c("beta", "delta", "lambdaStar", "marksMean", "marksPrecision"))

#' @rdname pompp_initial-class
#' @param x The pompp_initial object.
#' @param name The requested slot.
#' @return \strong{\code{`$`}}: The requested initial value (in case of
#' LambdaStar) or values (in case of Beta or Delta).
#' @export
#' @exportMethod $
methods::setMethod("$","pompp_initial",function(x,name) methods::slot(x, name))

#' @rdname pompp_initial-class
#' @param e1 A pompp_initial object.
#' @param e2 Another pompp_initial object or a list with pompp_initial
#' objects for \strong{+} and a positive integer for \strong{*}. e1 and e2
#' can be switched (+ and * are commutative).
#' @return \strong{\code{+}}: A list with the objects. Useful to start the
#' \code{fit_pompp} function, as it requires a list of initial values.
#' @export
#' @exportMethod +
methods::setMethod("+", "pompp_initial", function(e1, e2) list(e1, e2))

#' @rdname pompp_initial-class
#' @export
#' @exportMethod +
methods::setMethod("+", methods::signature(e1 = "list", e2 = "pompp_initial"),
                   function(e1, e2){
                     for (i in 1:length(e1))
                       if (!methods::is(e1[[i]], "pompp_initial"))
                         stop("Initial values can only be added to a list of other initial values.")
                     e1[[i + 1]] = e2
                     e1
                   })

#' @export
#' @exportMethod +
#' @rdname pompp_initial-class
methods::setMethod("+", methods::signature(e1 = "pompp_initial", e2 = "list"),
                   function(e1, e2){
                     l = list(e1)
                     for (i in 1:length(e2)){
                       if (!methods::is(e2[[i]], "pompp_initial")) stop("Initial values can only be added to a list of other initial values.")
                       l[[i + 1]] = e2[[i]]
                     }
                     l
                   })

#' @rdname pompp_initial-class
#' @param e1 A pompp_initial object.
#' @return \strong{\code{*}}: A list with \code{e2} random initial values.
#' @export
#' @exportMethod *
methods::setMethod("*", methods::signature(e1 = "pompp_initial", e2 = "numeric"),
                   function(e1, e2) {
                     if (e2 <= 0) stop("Can only multiply by a positive number.")
                     if (as.integer(e2) != e2) stop("Con only multiply by an integer.")
                     if (methods::slot(e1, "tag") != "random") {
                       message("Identical initial values are not recommended for independent Markov Chains.")
                       l = list()
                       for (i in 1:e2) l[[i]] = e1
                     } else {
                       l = list()
                       nb = length(methods::slot(e1, "beta"))
                       nd = length(methods::slot(e1, "delta"))
                       for (i in 1:e2) l[[i]] =
                         initial(nb, nd, methods::slot(e1, "lambdaStar"),
                                 methods::slot(e1, "marksMean"),
                                 methods::slot(e1, "marksPrecision"), TRUE)
                     }
                     l
                   })

#' @export
#' @exportMethod *
#' @rdname pompp_initial-class
methods::setMethod("*", methods::signature(e1 = "numeric", e2 = "pompp_initial"), function(e1, e2) e2 * e1)

#' @rdname pompp_initial-class
#' @param object A pompp_initial object.
#' @return \strong{\code{show}} and \strong{\code{print}}: The invisible object.
#' @export
#' @exportMethod show
methods::setMethod("show", "pompp_initial", function(object){
  cat("Initial values for a pompp model.")
  if (methods::slot(object, "tag") == "supplied") cat(" Values were supplied by user.\n\n")
  else if (methods::slot(object, "tag") == "random") cat(" Values were randomly generated.\n\n")

  cat("beta:\n")
  print(methods::slot(object, "beta"))
  cat("\ndelta:\n")
  print(methods::slot(object, "delta"))
  cat("\nlambdaStar:\n")
  print(methods::slot(object, "lambdaStar"))
  cat("\nmarksMu:\n")
  print(methods::slot(object, "marksMu"))
  cat("\nmarksPrecision:\n")
  print(methods::slot(object, "marksPrecision"))

  invisible(object)
})

#' @export
#' @exportMethod print
#' @param ... Currently unused.
#' @rdname pompp_initial-class
methods::setMethod("print", "pompp_initial", function(x, ...) methods::show(x))

#' @method print pompp_initial
#' @export
#' @rdname pompp_initial-class
print.pompp_initial <- function(x, ...) methods::show(x)

#' Initial values constructor for pompp modeling
#'
#' Helper function to create a valid set of initial values to be used with the
#' fit_pompp function.
#' @param beta Either a vector or a single integer. The vector is used if the
#' initial values are provided and the integer is used as the vector size to
#' be randomly generated.
#' @param delta Either a vector or a single integer. The vector is used if the
#' initial values are provided and the integer is used as the vector size to
#' be randomly generated.
#' @param lambdaStar A positive number.
#' @param marksMean Any real number. If random, defines the mean of the random
#' value.
#' @param marksPrecision A positive number. If random, defines the mean of the
#' random value.
#' @param random A logical value. If \code{TRUE}, then the initial values are
#' generated from standard normal distribution for \code{beta} and \code{delta}
#' and from a \code{Beta(lambdaStar, 1)} for \code{lambdaStar}. The latter is
#' generated as a low value due to potential explosive values resulting from
#' background area scaling.
#' @return A \code{pompp_initial} object. It can be used in the
#' \code{fit_pompp} function by itself, but must be in a list if multiple
#' initial values are supplied. Initial values can be combined by adding them
#' (with the use of `+`).
#' @seealso \code{\link{pompp_initial-class}}.
#' @examples
#' # Let us create initial values for a model with, say, 3 intensity covariates
#' # and 4 observability covariates. We add an initial values for both these
#' # cases due to the intercepts.
#'
#' # This first one is
#' in1 <- initial(rep(0, 4), c(0, 2, -1, -2, 3), 100, 0, 1)
#'
#' # Then we initalize some randomly.
#' in2 <- initial(4, 5, 100, 0, 1, random = TRUE)
#'
#' # We can even multiply the random one to generate more. Let us join them all
#' # to include in a model.
#' initial_values <- in1 + in2 * 3
#' # 4 chains are initialized.
#' @export
initial <- function(beta = numeric(), delta = numeric(), lambdaStar = numeric(),
                    marksMean = numeric(), marksPrecision = numeric(),
                    random = FALSE){
  if (random) methods::new("pompp_initial", beta = stats::rnorm(beta),
                           delta = stats::rnorm(delta),
                           lambdaStar = stats::rbeta(1, lambdaStar, 1), # Simulated small for safety
                           marksMean = stats::rnorm(1, marksMean),
                           marksPrecision = stats::rexp(1, 1 / marksPrecision),
                           tag = "random")
  else methods::new('pompp_initial',
                    beta = as.numeric(beta),
                    delta = as.numeric(delta),
                    lambdaStar = as.numeric(lambdaStar),
                    marksMean = as.numeric(marksMean),
                    marksPrecision = as.numeric(marksPrecision),
                    tag = "supplied")
}
