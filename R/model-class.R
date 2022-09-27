#' @include initial-class.R prior-class.R
NULL

#' Class that defines a model for the pompp package.
#'
#' The model includes the presence-only data, all selected variables, the link
#' functions for \eqn{q} and \eqn{p}, the initial values and the prior
#' distribution.
#' @field po The matrix containing the covariates values for the data.
#' @field intensityLink A string informing about the chosen link for the
#' intensity covariates. Current acceptable choice is only \code{"logit"}.
#' @field intensitySelection A vector containing the indexes of the selected
#' intensity columns in the \code{po} matrix.
#' @field observabilityLink A string informing about the chosen link for the
#' observability covariates. Current acceptable choice is only \code{"logit"}.
#' @field observabilitySelection A vector containing the indexes of the selected
#' observability columns in the \code{po} matrix.
#' @field marksSelection A single value containing the index of the selected
#' marks column in the \code{po} matrix.
#' @field coordinates A vector of two values containing the column positions
#' of the longitude and latitude in the \code{po} matrix.
#' @field init A list with objects of class \code{pompp_initial} indicating
#' the initial values for each chain. The length of this list tells the program
#' how many chains are requested to be run.
#' @field prior An object of class \code{pompp_prior} which indicates the
#' joint prior distribution for the model parameters.
#' @field iSelectedColumns If the intensity covariates selection was made with
#' the name of the columns, they are stored in this slot.
#' @field oSelectedColumns If the observability covariates selection was made
#' with the name of the columns, they are stored in this slot.
#' @field mSelectedColumns If the marks selection was made
#' with the name of the column, it is stored in this slot.
#' @seealso \code{\link{pompp_initial-class}} and
#' \code{\link{pompp_prior-class}} and \code{\link{pompp_model}}
#' @export
#' @exportClass pompp_model
methods::setClass("pompp_model",
                  methods::representation(po = "matrix",
                                          intensityLink = "character",
                                          intensitySelection = "numeric",
                                          observabilityLink = "character",
                                          observabilitySelection = "numeric",
                                          marksSelection = "numeric",
                                          coordinates = "numeric",
                                          init = "list",
                                          prior = "pompp_prior",
                                          iSelectedColumns = "character",
                                          oSelectedColumns = "character",
                                          mSelectedColumns = "character",
                                          cSelectedColumns = "character"
                  ),
  validity = function(object){
    s <- function(n) methods::slot(object, n)
    if (length(s("intensityLink"))>1)
      stop("Argument intensityLink must have length 1.")
    if (length(s("observabilityLink"))>1)
      stop("Argument observabilityLink must have length 1.")
    validLinks = c("logit")
    if (!(s("intensityLink") %in% validLinks))
      stop(paste0(s("intensityLink")," is not a valid link. Accepted options for current version are\n",validLinks))
    if (!(s("observabilityLink") %in% validLinks)) stop(paste0(s("observabilityLink")," is not a valid link. Accepted options for current version are\n", validLinks))
    if (any(s("intensitySelection") != as.integer(s("intensitySelection")))) stop("intensitySelection must be integers.")
    if (any(s("observabilitySelection") != as.integer(s("observabilitySelection")))) stop("observabilitySelection must be integers.")
    if (any(s("coordinates") != as.integer(s("coordinates")))) stop("coordinates must be integers.")
    if (length(s("coordinates")) != 2) stop("coordinates must have two elements.")
    nbSel = length(s("intensitySelection")) + 1; ndSel = length(s("observabilitySelection")) + 2
    if ((length(s("marksSelection")) != 1) || (s("marksSelection") != as.integer(s("marksSelection")))) stop("marksSelection must be a single integer.")
    if (length(s("coordinates")) != 2) stop("coordinates must contain 2 values: the columns of the coordinates. Longitude and latitude, in this order.")
    for (i in s("init")){
      if (!is(i,"pompp_initial")) stop("Initial values must be constructed with the initial function.")
      if (length(methods::slot(i,"beta")) != nbSel) stop(paste("\nInitial values for beta has the wrong size. Expected size:", nbSel))
      if (length(methods::slot(i,"delta")) != ndSel) stop(paste("\nInitial values for delta has the wrong size. Expected size:", ndSel))
    }
    if (length(methods::slot(methods::slot(s("prior"),"beta"),"mu")) != nbSel) stop(paste("Prior for beta has wrong sized parameters. Expected size:", nbSel))
    if (length(methods::slot(methods::slot(s("prior"),"delta"),"mu")) != ndSel) stop(paste("Prior for delta has wrong sized parameters. Expected size:", ndSel))
    TRUE
                  })

#' @rdname pompp_model-class
#'
#' @param x The pompp_model object.
#' @return \strong{\code{names}}: A character vector with possible options
#' for the \code{`$`} and \code{`$<-`} methods.
#' @export
#' @exportMethod names
methods::setMethod("names", "pompp_model",
                   function(x) c("po", "intensityLink", "intensitySelection",
                                 "observabilityLink", "observabilitySelection",
                                 "marksSelection", "coordinates", "initial_values", "prior"))

#' @rdname pompp_model-class
#'
#' @param x The pompp_model object.
#' @param name The requested slot. Available options are not necessarily the
#' class slots. They can be checked with the \code{names} method.
#' @return \strong{\code{`$`}}: The requested slot's value.
#' @export
#' @exportMethod $
methods::setMethod("$","pompp_model",function(x,name){
  if (name %in% c("initial_values", "initial")) methods::slot(x, "init") else
    if (name == "coordinates") methods::slot(x, "po")[, methods::slot(x, "coordinates")] else
      methods::slot(x, name)
})

#' @rdname pompp_model-class
#'
#' @param x The pompp_model object.
#' @param name The requested slot.
#' @param value New value.
#' @return \strong{\code{`$<-`}}: The new object with the updated slot.
#' @export
#' @exportMethod $<-
methods::setMethod("$<-","pompp_model",function(x, name, value){
  if (name %in% c("initial values", "initial")) methods::slot(x, "init") = value else methods::slot(x, name) = value
  methods::validObject(x)
  x
})

#' @rdname pompp_model-class
#' @param object The pompp_model object.
#' @return \strong{\code{show}} and \strong{\code{print}}: The invisible object.
#' @export
#' @exportMethod show
methods::setMethod("show", "pompp_model", function(object){
  s <- function(name) methods::slot(object, name)
  cat("Model for the pompp package.\n\n")
  cat(nrow(s("po")), "presence-only locations are included.\n")
  cat(length(s("intensitySelection")), "variables were selected for the intensity set. ")
  if (length(s("iSelectedColumns"))){
    cat("They were:\n")
    cat(s("iSelectedColumns"), sep = ", ")
  }
  cat("\nThe intensity link is", s("intensityLink"), "and the effects prior is",
      methods::slot(methods::slot(s("prior"), "beta"), "family"), "\b.\n")
  cat(length(s("observabilitySelection")), "variables were selected for the observability set. ")
  if (length(s("oSelectedColumns"))){
    cat("They were:\n")
    cat(s("oSelectedColumns"), sep = ", ")
  }
  cat("\nThe observability link is", s("intensityLink"), "and the effects prior is",
      methods::slot(methods::slot(s("prior"), "delta"), "family"), "\b.\n\n")
  if (length(s("mSelectedColumns")))
    cat("The column selected for the marks was", s("mSelectedColumns"))
  supplied <- 0; random <- 0
  for (i in 1:length(s("init"))) {
    if (methods::slot(s("init")[[i]], "tag") == "supplied") supplied <- supplied + 1
    if (methods::slot(s("init")[[i]], "tag") == "random") random <- random + 1
  }
  cat(length(s("init")), " chains were initialized where ",
      ifelse(supplied, paste(supplied, "are supplied values"), ""),
      ifelse(supplied * random, " and ", ""),
      ifelse(random, paste(random, "are random values"), ""), ".\n", sep = "")
})

#' @export
#' @exportMethod print
#' @param ... Currently unused.
#' @rdname pompp_model-class
methods::setMethod("print", "pompp_model", function(x, ...) methods::show(x))

#' @method print pompp_model
#' @export
#' @rdname pompp_model-class
print.pompp_model <- function(x, ...) methods::show(x)


