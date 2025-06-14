library(data.table)
library(mlr3)
library(mlr3misc)


named_union = function(x, y) set_names(union(x, y), union(names(x), names(y)))
mlr_reflections$task_feature_types = named_union(mlr_reflections$task_feature_types, c(list="list"))


gen_ellipse = function(r=1, nSize=100, sdDev=0.1) {
  # generate point cloud
  theta = runif(nSize,min = 0,max = 2*pi) # generate angles from Unif(0,2*pi)
  X = matrix(nrow = nSize,ncol = 2)
  X[,1] = cos(theta)+rnorm(nSize,0,sdDev)
  X[,2] = r*sin(theta)+rnorm(nSize,0,sdDev)
  return(X)
}


#' Creates ellepse regression task
#' 
#' This function creates MLR3 task with squzzed ellipse regression problem
#' 
#' @param N number of points
#' @param sdDev specifies how far are points from the ellipse
#' @return MLR3 task object with r (list of squzzed ratios) and X (list of point clouds) data fields
#' @export
tsk_ellipse <- function(N=100, sdDev = 0.1) {
  rList = runif(N)
  data = data.table(r=rList)
  data$X = lapply(rList, gen_ellipse, sdDev = sdDev)
  task = as_task_regr(data, target = "r")
  task$characteristics = list(sdDev = sdDev)
  return(task)
}


gen_CR <- function(type, nPoints = 100, r0 = 1, dr = 0.1, sdDev = 0.1) {
  theta <- runif(nPoints, min = 0,max = 2*pi) # generate angles from Unif(0,2*pi)
  X <- matrix(nrow = nPoints, ncol = 2)
  
  if (type==1){
    r <-  rnorm(nPoints,mean=r0,sd=sdDev)
  } else{
    means <- sample(c(r0-dr, r0+dr), nPoints,replace = TRUE)
    r <-  rnorm(nPoints, mean=means,sd=sdDev)
  }
  X[,1] <- r*cos(theta)
  X[,2] <- r*sin(theta)
  return(X)
}

#' Creates Circle vs Ring Problem
#' 
#' This ffunction MLR3 task with circle vs ring classification problem
#' 
#' Each data entry is a cloud of point scattered around either a single circle with radius r0
#'  or a ring (i.e. two circles with radia r0+dr and r0-dr)
#' The problem would be to guess type of the point cloud
#' 
#' @param N number of point clouds
#' @param nPoints number of points in each point cloud
#' @param r0 radius of the enter circle
#' @param dr half of the distance between two circles in the ring
#' @param sdDev measures the closeness of points to original shapes
tsk_CvsR <- function(N=100, nPoints = 100, r0 = 1, dr = 0.1, sdDev = 0.1) {
  pb = progress::progress_bar$new(total = N)
  typeList = sample( c(1, 2), N, replace = TRUE)
  data = data.table(type=typeList)
  data$X = lapply(typeList, gen_CR, nPoints = nPoints, r0 = r0, dr = dr, sdDev = sdDev)
  task = as_task_classif(data, target = "type")
  task$characteristics = list(sdDev = sdDev)
  return(task)
}

