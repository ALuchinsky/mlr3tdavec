# library(TDAvec)
# library(R6)
# library(mlr3)
# library(mlr3pipelines)
# library(mlr3misc)
# library(mlr3learners)
# library(paradox)
# library(mlr3viz)
# library(mlr3tuning)
# library(TDA)
# library(data.table)

.onLoad <- function(libname, pkgname) {
  packageStartupMessage("Loading mlr3tdavec package")
  if (!requireNamespace("mlr3", quietly = TRUE)) return(invisible())
  refl <- get("mlr_reflections", envir = asNamespace("mlr3"))
  # IMPORTANT: keep task_feature_types a named *character* vector (atomic)
  refl$task_feature_types["list"] <- "list"
  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Attaching mlr3tdavec package")
}


#=========== Elipse Generation
gen_ellipse = function(r=1, nSize=100, sdDev=0.1) {
  # generate point cloud
  theta <- stats::runif(nSize,min = 0,max = 2*pi) # generate angles from Unif(0,2*pi)
  X <- matrix(nrow = nSize,ncol = 2)

  X[,1] <- cos(theta)+stats::rnorm(nSize,0,sdDev)
  X[,2] <- r*sin(theta)+stats::rnorm(nSize,0,sdDev)
  return(X)
}

make_clear_PD <- function(X, maxdimension = 1, maxscale = 2) {
  PD = TDA::ripsDiag(X, maxdimension = 1, maxscale = 2)$diagram
  PD = PD[-1,]
  return(PD)
}

tsk_ellipse <- function(N=100, sdDev = 0.1, nSize = 100) {
  rList = stats::runif(N)
  data <- data.table::data.table(r=rList)
  data$X = lapply(rList, gen_ellipse, sdDev = sdDev, nSize = nSize)
  task = mlr3::as_task_regr(data, target = "r")
  task$characteristics = list(sdDev = sdDev)
  return(task)
}


gen_CR <- function(type, nPoints = 100, r0 = 1, dr = 0.1, sdDev = 0.1) {
  theta <- stats::runif(nPoints, min = 0,max = 2*pi) # generate angles from Unif(0,2*pi)
  X <- matrix(nrow = nPoints, ncol = 2)

  if (type==1){
    r <-  stats::rnorm(nPoints,mean=r0,sd=sdDev)
  } else{
    means <- sample(c(r0-dr, r0+dr), nPoints,replace = TRUE)
    r <-  stats::rnorm(nPoints, mean=means,sd=sdDev)
  }
  X[,1] <- r*cos(theta)
  X[,2] <- r*sin(theta)
  return(X)
}

tsk_CvsR <- function(N=100, nPoints = 100, r0 = 1, dr = 0.1, sdDev = 0.1) {
  typeList = sample( c(1, 2), N, replace = TRUE)
  data = data.table::data.table(type=typeList)
  data$X = lapply(typeList, gen_CR, nPoints = nPoints, r0 = r0, dr = dr, sdDev = sdDev)
  task = mlr3::as_task_classif(data, target = "type")
  task$characteristics = list(sdDev = sdDev)
  return(task)
}



generate_circle_ring  <- function(nSamples = 100, sdDev = 0.15, r0 = 1, dr = 0.1, nPoints = 100) {
  PD_circle_ring <- function(nSize = nPoints,grp,sdDev){
    # generate point cloud
    #.  circle if grp=1, ring (2 concentric circles) if grp=2
    theta <- stats::runif(nSize,min = 0,max = 2*pi) # generate angles from Unif(0,2*pi)
    X <- matrix(nrow = nSize,ncol = 2)

    if (grp==1){
      r <-  stats::rnorm(nSize,mean=r0,sd=sdDev)
    } else{
      means <- sample(c(r0-dr, r0+dr),nSize,replace = TRUE)
      r <-  stats::rnorm(nSize,mean=means,sd=sdDev)
    }

    X[,1] <- r*cos(theta)
    X[,2] <- r*sin(theta)
    # compute PD
    TDAstats::calculate_homology(X,dim = 1,threshold = 2.2+3*sdDev)
  }

  # creates PDs of n1 circles followed by n2 rigns
  generatePD <- function(nSize=nPoints,n1 = nSamples/2,n2=nSamples/2,sdDev){
    # generate PDs
    n <- n1+n2
    D <- list()
    for (k in 1:n)
      if (k<=n1) D[[k]] <- PD_circle_ring(nSize,grp = 1,sdDev) else
        D[[k]] <- PD_circle_ring(nSize,grp = 2,sdDev)
    return(D)
  }

  Y <- factor( c(rep(1, nSamples/2), rep(2, nSamples/2)))
  PDs <- generatePD(sdDev = sdDev)
  # generate data
  return(list(
    Y = factor(Y),
    PDs = PDs,
    IDs = rep(1, length(Y)),
    params = list(data="circle_ring", sdDev=sdDev, r0=r0, dr=dr)
  ))
}


#============= Pipelines

PipeOpTDA_PD <- R6::R6Class(
  inherit = mlr3pipelines::PipeOpTaskPreprocSimple,

  public = list(
    featureName = "X",
    pdName = "PD",
    homDim = 1,
    maxScale = 1,
    #
    initialize = function(id = "custom", featureName = "X", pdName = "PD") {
      self$featureName = featureName
      self$pdName = pdName
      ps = paradox::ParamSet$new(params = list(
        homDim = paradox::p_int(lower = 0, upper = 1),
        maxScale = paradox::p_dbl(lower = 0)
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
      TT = data.table::data.table(col = pds)
      colnames(TT) = self$pdName
      task$cbind(TT)
      return(task)
    }
  )
)



PipeOpTDAVec = R6::R6Class("PipeOpTDAVec",
                       inherit = mlr3pipelines::PipeOpTaskPreprocSimple,
                       public = list(
                         pdsName = "PD",
                         outName = NULL,
                         scaleSeq = NULL,
                         featureName = "PD",
                         nSeq = 10,
                         initialize = function(id = "custom", vectName = NULL, outName = NULL, homDim = NULL, nGrid=NULL, K=NULL, w=0) {
                           ps = paradox::ParamSet$new(params = list(
                              homDim = paradox::p_int(lower = 0, upper = 1),
                              nGrid = paradox::p_int(lower = 3, upper = 200),
                              vectName = paradox::p_fct(c("PL", "PS", "BC", "FDA")),
                              K=paradox::p_int(lower = 1, upper = 100),
                              w=paradox::p_int(lower=0, upper = 5)
                           ))
                           if(!is.null(vectName)) {
                             ps$values$vectName = vectName
                             outName = ifelse(is.null(outName), vectName, outName)
                           }
                           if(!is.null(K)) {
                             ps$values$K = K
                           }
                           if(!is.null(homDim)) { ps$values$homDim = homDim}
                           if(!is.null(nGrid)) { ps$values$nGrid = nGrid}
                           if(!is.null(w)) { ps$values$w = w}
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
                         #==================================
                         .findLimits = function(D,homDim){
                           # D: a list of PDs
                           n = length(D)
                           findMax = function(pd,homDim){
                             pd = pd[pd[,1]==homDim,2:3,drop=F] # extracting PD of dimension d
                             if (nrow(pd)>0) c(min(pd[,1]),
                                               max(pd[is.finite(pd[,2]),2]))
                             else rep(NA,2)
                           }
                           # body
                           minB = numeric(length = n)
                           maxD = numeric(length = n)
                           for (k in 1:n){
                             ans = findMax(D[[k]],homDim)
                             minB[k] = ans[1]
                             maxD[k] = ans[2]
                           }
                           ans = c(min(minB,na.rm = T),
                                    max(maxD,na.rm = T))
                           names(ans) = c('minB','maxD')
                           return(ans)
                         },
                         #=========================================================
                         # Computes FDA coefficients for the given PD
                         .generateForierFdata = function(D, homDim, K,
                                                         w = 0,
                                                         type='betti'){

                           wFunc = function(b,d) (d-b)**w
                           n = length(D)
                           X = array(dim = c(n, 2*K + 1))
                           bound = private$.findLimits(D, homDim)
                           for (k in 1:n) {
                             pd = D[[k]]
                             pd = pd[pd[,1]==homDim,2:3,drop=F]
                             b = pd[,1]/bound['maxD']; d = pd[,2]/bound['maxD']
                             alpha = d-b # persistence

                             X[k,1] = sum(wFunc(b, d)*alpha)

                             for (m in 1:K){
                               c = 2*m*pi
                               if(type == "betti") {
                                 alphasin = sin(c*d)-sin(c*b)
                                 alphacos = cos(c*d)-cos(c*b)
                                 X[k,2*m] = -sqrt(2)/c*sum(wFunc(b, d)*alphacos)
                                 X[k,2*m+1] = sqrt(2)/c*sum(wFunc(b, d)*alphasin)
                               }
                               else if(type == "silhouette") {
                                 alphacos = sin(c*b) + sin(c*d) - 2*sin(c*(b+d)/2)
                                 alphasin = cos(c*b) + cos(c*d) - 2*cos(c*(b+d)/2)
                                 cc = -sqrt(2)/(c*c)
                                 X[k,2*m] = cc*sum(wFunc(b, d)*alphacos)
                                 X[k,2*m+1] = cc*sum(wFunc(b, d)*alphasin)
                               }
                               else {
                                 cat("ERROR: unknown type ", type,"\n")
                               }
                             }
                           }
                           return(X)
                         }
                         ,
                         .transform = function(task) {
                           homDim = self$param_set$values$homDim
                           pds = task$data(cols = self$pdsName)[[1]]
                           if(self$param_set$values$vectName =="FDA") {
                             TT = private$.generateForierFdata(pds, homDim = homDim, K=self$param_set$values$K, w=self$param_set$values$w)
                           } else {
                             func = switch (self$param_set$values$vectName,
                                            "PL" = TDAvec::computePersistenceLandscape,
                                            "PS" = TDAvec::computePersistenceSilhouette,
                                            "BC" = TDAvec::computeBettiCurve
                             )
                             nSeq = self$param_set$values$nGrid
                             maxD = TDAvec::computeLimits(pds, homDim)[3]
                             self$scaleSeq = seq(0, maxD, length.out = nSeq)
                             TT = t(sapply(pds, func, homDim = homDim, scaleSeq = self$scaleSeq))
                           }
                           TT = data.table::data.table(col = TT)
                           colnames(TT) = gsub("col", self$outName, colnames(TT))
                           task$cbind(TT)
                           return(task)
                         }
                       )
)





