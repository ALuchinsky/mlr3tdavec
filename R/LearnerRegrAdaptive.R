LearnerRegrAdaptive = R6::R6Class(
  "LearnerRegrMy",
  inherit = LearnerRegr,

  public = list(
    initialize = function(lambda = 0) {
      ps = ps(
        lambda = paradox::p_dbl(lower = 0, upper = 10, default = 0)
        # add hyperparameters here (paradox 1.0.1 style)
        # alpha = p_dbl(lower = 0, upper = 1, default = 0.5)
      )
      if(!is.null(lambda)) {
        ps$values$lambda = lambda
      }
      super$initialize(
        id = "regr.my",
        feature_types = c("numeric", "integer", "logical", "factor", "ordered", "character"),
        predict_types = "response",           # add "se" etc if you support it
        param_set = ps,
        properties = character(0)             # e.g. c("weights") if supported
      )
    }
  ),

  private = list(
    .createF = function(out) {
      target_name = out$target_names[1]
      data = out$data()
      feature_cols = which(colnames(data) != target_name)
      fea_data = data[,feature_cols, with = FALSE]
      cols = colnames(fea_data)
      classes = unique(sapply(cols, function(col) strsplit(col[1], "\\.")[[1]][1]))
      data_list = list()
      for(cl in classes) {
        cols_ = names(fea_data)[startsWith(names(fea_data), cl)]
        XX = fea_data[,cols_, with = FALSE]
        keep <- apply(XX, 2, function(col) length(unique(col)) > 1)
        data_list[cl] <- list(as.matrix(XX[, keep, drop = FALSE, with = FALSE]))
      }
      dims = c(nrow(data_list[[1]]), max(sapply(data_list, ncol)), length(data_list))
      F = array(0, dims)
      for(j in seq_along(data_list)) {
        D = data_list[[j]]
        F[,1:ncol(D), j] = D
      }
      return(F)
    },

    .computeLoss = function(par,Loss=0,nIter=1,maxIter=100,relTol=0.001,lambda, Ftrain,rtrain,R){
      par = par/sum(abs(par))
      dims <- dim(Ftrain)
      X <- matrix(1,nrow = dims[1],ncol=dims[2]+1)
      for (m in 1:dims[2]) X[,m+1] <- Ftrain[,m,]%*%par
      betaHat <- solve(t(X)%*%X+lambda*R)%*%t(X)%*%rtrain
      yHat <- X %*% betaHat

      C <- crossprod(t(betaHat[-1]))
      LossUpd <- sum((rtrain - yHat)^2)+lambda*sum(C*R[-1,-1])


      Z <- matrix(nrow = dims[1],ncol = dims[3])
      for (j in 1:dims[3]) Z[,j] <- Ftrain[,,j]%*%betaHat[-1]
      gammaHat <- solve(t(Z)%*%Z)%*%t(Z)%*%(rtrain-betaHat[1])

      h = list(gammaHat = gammaHat, betaHat = betaHat, loss = LossUpd)
      self$state$history_ <- append(self$state$history_, list(h))

      if ((abs(1-LossUpd/Loss)<relTol)|(nIter>=maxIter))
        return(gammaHat)
      else{
        gammaNorm=sqrt(sum(gammaHat^2))
        betaNorm=sqrt(sum(betaHat[-1]^2))
        private$.computeLoss(par=gammaHat,Loss=LossUpd,nIter=nIter+1,maxIter,relTol,lambda,Ftrain,rtrain,R)
      }
    },


    .train = function(out) {
      F <- private$.createF(out)
      target_name = out$target_names[1]
      target = as.numeric(out$data()[, target_name[[1]], with = FALSE])
      self$state$history_ =  c()
      private$.computeLoss( runif(dim(F)[3]), Loss = 0, nIter = 1, maxIter = 100, relTol = 0.001,
                            Ftrain = F, rtrain = target, lambda = self$param_set$values$lambda, R = diag(dim(F)[2]+1))
      fit = self$state$history_[[ length(self$state$history_)]]
      invisible(fit)
    },

    .predict = function(task) {
      Ftrain <- private$.createF(task)
      dims = dim(Ftrain)
      X <- matrix(1,nrow = dims[1],ncol=dims[2]+1)
      for (m in 1:dims[2]) X[,m+1] <- Ftrain[,m,] %*% self$model$gammaHat
      yHat <- X %*% self$model$betaHat
      PredictionRegr$new(task = task, response = as.numeric(yHat))
    }
  )
)
