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



