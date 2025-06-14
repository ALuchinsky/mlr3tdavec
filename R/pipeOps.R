library(TDAvec)
library(R6)
library(mlr3)
library(mlr3pipelines)
library(mlr3misc)
library(mlr3learners)
library(paradox)
library(mlr3viz)
library(mlr3tuning)
library(TDA)
library(data.table)

make_clear_PD <- function(X, maxdimension = 1, maxscale = 2) {
  PD = ripsDiag(X, maxdimension = 1, maxscale = 2)$diagram
  PD = PD[-1,]
  return(PD)
}


#' Percistence Diagram Creation PipeOp
#' 
#' @export
PipeOpTDA_PD <- R6Class(
  inherit = PipeOpTaskPreprocSimple,
  
  public = list(
    featureName = "X",
    pdName = "PD",
    homDim = 1,
    maxScale = 1,
    #  
    initialize = function(id = "custom", featureName = "X", pdName = "PD") {
      self$featureName = featureName
      self$pdName = pdName
      ps = ParamSet$new(params = list(
        homDim = p_int(lower = 0, upper = 1),
        maxScale = p_dbl(lower = 0)
      ))
      super$initialize(
        id=id,
        param_set = ps,
        packages = "mlr3pipelines",
        feature_types = c("numeric")
      )
      self$featureName = featureName
    }
  ),
  
  private = list(
    .transform = function(task) {
      X = as.list(task$data(cols = self$featureName))[[1]]
      pds = lapply(X, make_clear_PD, maxdimension=self$homDim, maxscale = self$maxScale)
      TT = data.table(col = pds)
      colnames(TT) = self$pdName
      task$cbind(TT)
      return(task)
    }
  )
)
