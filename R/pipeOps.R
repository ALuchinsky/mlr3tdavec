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

#' @export
PipeOpTDAVec = R6Class("PipeOpTDAVec",
                       inherit = PipeOpTaskPreprocSimple,
                       public = list(
                         pdsName = "PD",
                         outName = NULL,
                         scaleSeq = NULL,
                         featureName = "PD",
                         nSeq = 10,
                         initialize = function(id = "custom", outName = NULL) {
                           ps = ParamSet$new(params = list(
                             homDim = p_int(lower = 0, upper = 1),
                             nGrid = p_int(lower = 3, upper = 200),
                             vectName = p_fct(c("PL", "PS", "BC"))
                           ))
                           super$initialize(
                             id = id,
                             param_set = ps,
                             packages = "mlr3pipelines",
                             feature_types = c("numeric")
                           )
                           self$outName = ifelse(is.null(outName), id, outName)
                         }
                       ),
                       private = list(
                         .transform = function(task) {
                           homDim = self$param_set$values$homDim
                           nSeq = self$param_set$values$nGrid
                           pds = task$data(cols = self$pdsName)[[1]]
                           maxD = computeLimits(pds, homDim = homDim)[3]
                           self$scaleSeq = seq(0, maxD, length.out = nSeq)
                           func = switch (self$param_set$values$vectName,
                                          "PL" = computePersistenceLandscape,
                                          "PS" = computePersistenceSilhouette,
                                          "BC" = computeBettiCurve
                           )
                           # cat("PipeOpTDAVec: transfer: homDim=", homDim, 
                           #     ", nSeq=", nSeq, "[scaleSeq]=", length(self$scaleSeq), "\n")
                           TT = t(sapply(pds, func, homDim = homDim, scaleSeq = self$scaleSeq))
                           TT = data.table(col = TT)
                           colnames(TT) = gsub("col", self$outName, colnames(TT))
                           task$cbind(TT)
                           return(task)
                         }
                       )
)


