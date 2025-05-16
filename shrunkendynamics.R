library(Rcpp)
sourceCpp("llik.cpp")
library(tictoc)
library(foreach)
library(doParallel)
library(doRNG)
library(optimParallel)
nr_free_cores = max(1, detectCores() - 1)
registerDoParallel(cores = nr_free_cores)
library(latex2exp)
library(ggplot2)
library(ggpubr)

## set dimensions and nr of snapshots
nr_replications = 10
dim_T = 1500
q_proportion = 0.5
meanidx = 4 # 1: mean + sd, 4: median + mad, 7: quantiles
dim_N = round(q_proportion * dim_T)
dim_snap = min(dim_T, 50)

## shrink routine
shrink4 = function(lambda, q, h = NA) {
  lambdai = 1/lambda
  if (q <= 1) {
    # following Ledoit and Wolf (2022, Bernoulli 28, 1519-1547)
    if (is.na(h)) h = min(q*q, 1/q/q)^0.35 / length(lambda)^0.35
    h2 = h * h

    aid = outer(lambdai, lambdai, "-")
    aidnum = aid * aid + h2 * lambdai * lambdai
    ## eq (4.13)
    thetahat = colMeans(lambdai * (aid / aidnum))
    ## eq (4.14)
    calA = thetahat^2 + colMeans((h * lambdai * lambdai) / aidnum)^2
    mui = (1 - q)^2 * lambdai +
      ## eq (4.12)+(4.13)
      2 * q * (1 - q) * lambdai * thetahat+
      ## eq (4.12)+(4.14)
      q^2 * lambdai * calA
    ## add refinement (4.15)
    mu = (1 / mui); mu = mu / sum(mu) * sum(lambda)
    return(mu)
  } else {
    stop("nog niet geimplementeerd")
  }
}

## generate true spectra
real_spectra_high = function(x) 5 * (2 + sin(2*pi * 2 * x))
real_spectra_low  = function(x) 2 * (1.1 - sin(2*pi * 2 * x))
real_spectra_fn = function(x, dim_N) outer(real_spectra_high(x) - real_spectra_low(x), 1:dim_N / dim_N) + real_spectra_low(x)




llik = function(Wyt, omega, a, b, f1) {
  dim_T = nrow(Wyt)
  dim_N = ncol(Wyt)
  fnow = f1
  # as = a[1] + (a[2] - a[1]) * (0:(dim_N-1))/(dim_N - 1)
  # bs = b[1] + (b[2] - b[1]) * (0:(dim_N-1))/(dim_N - 1)
  # omega = omega[1] + (omega[2] - omega[1]) * (0:(dim_N-1))/(dim_N - 1)
  ft = matrix(0, nrow = dim_T, ncol = dim_N)
  likt = rep(0, dim_T)
  for (i1 in 1:dim_T) {
    ft[i1, ] = fnow
    aid1 = Wyt[i1, ]^2
    likt[i1] = - 0.5 * sum(log(fnow)) - 0.5 * sum(aid1 / fnow)
    fnow = (1-b) * omega + b * fnow + a * (aid1 - fnow)
    if (min(fnow) <= 0) return(list(llik = -1e10))
  }
  likt = likt - 0.5 * dim_N * log(2*pi)
  llik = mean(likt)/dim_N
  if (is.infinite(llik) | is.nan(llik)) return(list(llik = -1e10))
  return(list(
    llik = llik,
    likt = likt,
    ft = ft
  ))
}



llik1 = function(Wyt, omega, a, b, f1) {
  if (max(abs(c(a,b))) > 10) return(list(llik = -1e10))
  a = 1 / (1+exp(-a))
  b = 1 / (1+exp(-b))
  dim_T = length(Wyt)
  fnow = f1
  likt = ft = rep(0, dim_T)
  for (i1 in 1:dim_T) {
    ft[i1] = fnow
    aid1 = Wyt[i1]^2
    likt[i1] = -0.5 * log(2*pi*fnow) - 0.5 * aid1/fnow
    fnow = (1-b) * omega + b * fnow + a * (aid1 - fnow)
    if (fnow <= 0) return(list(llik = -1e10))
  }
  llik = mean(likt)
  if (is.infinite(llik) | is.nan(llik)) llik = -1e10
  return(list(
    llik = llik,
    likt = likt,
    ft = ft
  ))
}




summarize_replications = function(in_matrix, upto = NULL) {
  if (is.null(upto)) upto = ncol(in_matrix) else upto = median(c(1, upto, ncol(in_matrix)))
  dim_T = nrow(in_matrix)
  foreach (iter0=1:dim_T, .combine = rbind) %dopar% {
    aid1 = in_matrix[iter0, 1:upto]
    c(
      mean(aid1) + c(0,-1.5,1.5) * sd(aid1), 
      median(aid1) + c(0,-1.5,1.5) * mad(aid1), 
      quantile(aid1, c(0.5, 0.2, 0.8))
    )
  }
}



make_frame = function(dim_snap, estimated_spectra, 
                      Gas_ft_1_summ, Gas_ft_N_summ, 
                      GasS_ft_1_summ, GasS_ft_N_summ, 
                      meanidx = 1, withplot = FALSE, withband = FALSE,
                      withsqrtn = TRUE, whichone = "max", ymax = 50) {
  idx = 0:dim_snap / dim_snap
  my_df = rbind(
    data.frame(x = idx, y = real_spectra_high(idx), ylw = NA, yup = NA, name = "population", col = "black", lty = "dashed", mm = "max"),
    data.frame(x = idx, y = real_spectra_low(idx), ylw = NA, yup = NA, name = "population", col = "black", lty = "solid", mm = "min"),
    data.frame(x = idx, y = estimated_spectra[ , 1], ylw = NA, yup = NA, name = "population biased", col = "red", lty = "dashed", mm = "max"),
    data.frame(x = idx, y = estimated_spectra[ , dim_N], ylw = NA, yup = NA, name = "population biased", col = "red", lty = "solid", mm = "min"),
    # data.frame(x = idx, y = shrunken_spectra[ , 1], name = "shr_max", col = "blue", lty = "dashed", mm = "max"),
    # data.frame(x = idx, y = shrunken_spectra[ , dim_N], name = "shr_min", col = "blue", lty = "solid", mm = "min"),
    data.frame(x = 1:dim_T / dim_T, y   = Gas_ft_1_summ[ , meanidx + 0], ylw = Gas_ft_1_summ[ , meanidx + 1], yup = Gas_ft_1_summ[ , meanidx + 2], name = "shrink intercept", col = "forestgreen", lty = "dashed", mm = "max"),
    data.frame(x = 1:dim_T / dim_T, y   = Gas_ft_N_summ[ , meanidx + 0], ylw = Gas_ft_N_summ[ , meanidx + 1], yup = Gas_ft_N_summ[ , meanidx + 2], name = "shrink intercept", col = "forestgreen", lty = "dashed", mm = "min"),
    data.frame(x = 1:dim_T / dim_T, y   = GasS_ft_1_summ[ , meanidx + 0], ylw = GasS_ft_1_summ[ , meanidx + 1], yup = GasS_ft_1_summ[ , meanidx + 2], name = "shrink path", col = "purple", lty = "dashed", mm = "max"),
    data.frame(x = 1:dim_T / dim_T, y   = GasS_ft_N_summ[ , meanidx + 0], ylw = GasS_ft_N_summ[ , meanidx + 1], yup = GasS_ft_N_summ[ , meanidx + 2], name = "shrink path", col = "purple", lty = "dashed", mm = "min")
  )
  if (withplot) {
    gg1 = ggplot(data = subset(my_df, mm == whichone), aes(x = x, y = y, color = name)) +
      geom_line(aes(color = name), size = 1) + xlab("time") + coord_cartesian(ylim = c(0,ymax))
    if (withband) gg1 = gg1 + geom_ribbon(aes(ymin=ylw, ymax=yup, fill = name), alpha = 0.1)
    gg1 = gg1 +
      theme(axis.text=element_text(size=16), axis.title=element_text(size=18,face="bold"),
            legend.text = element_text(size=16), legend.title = element_text(size=25)) +
      ylab(TeX(paste0("$\\lambda^{\\", whichone, "}$"))) + 
      labs(col=unname(TeX(paste0("$\\lambda^{\\", whichone, "}$")))) +
      guides(fill = "none")
    plot(gg1)
  }
  return(gg1)
}


#######################################
#######################################
##
## START OF SIMULATION AND PRESENTATION
##
#######################################
#######################################

Gas_ft_1 = Gas_ft_N = GasS_ft_1 = GasS_ft_N = GasSS_ft_1 = GasSS_ft_N = 
  matrix(0, nrow = dim_T, ncol = nr_replications)
shrunken_spectra = estimated_spectra = 0


for (iter0 in 1:nr_replications) {
  ## generate empirical and shrunken spectra
  tic()
  yt = foreach (i1=1:dim_T, .combine = rbind) %dorng% {
    c(i1, sqrt(real_spectra_fn(i1 / dim_T, dim_N)) * rnorm(dim_N))
  }
  yt = yt[order(yt[ , 1]), 1 + 1:dim_N]
  real_spectra = real_spectra_fn(0:dim_snap / dim_snap, dim_N)
  estimated_spectra_tmp = foreach (i1=0:dim_snap, .combine = rbind) %dorng% {
    aid1 = sqrt(real_spectra[i1+1, ]) * matrix(rnorm(dim_T * dim_N), nrow = dim_N, ncol = dim_T)
    aid1 = aid1 %*% t(aid1) / dim_T
    lambda = eigen(aid1, symmetric = TRUE)$values
    c(i1, lambda, shrink4(lambda, q_proportion))
  }
  estimated_spectra_tmp = estimated_spectra_tmp[ order(estimated_spectra_tmp[ , 1]), ]
  shrunken_spectra = ((iter0 - 1) * shrunken_spectra + estimated_spectra_tmp[ , 1:dim_N + 1 + dim_N]) / iter0
  estimated_spectra = ((iter0 - 1) * estimated_spectra + estimated_spectra_tmp[ , 1:dim_N + 1]) / iter0
  toc()
  
  cat(iter0, " ... Eigenvector projections\n")
  tic()
  aid1 = eigen(cov(yt))
  mW = aid1$vectors
  yW = yt %*% mW
  lambda = aid1$values
  mu = shrink4(lambda, dim_N/dim_T)
  toc()
  
  
  cat(iter0, " ... Univariate fits (debiased targeting)\n")
  tic()
  Gas_ft = foreach(i1=1:dim_N , .combine = cbind) %dopar% {
    x0 = c(-3,3)
    estim = optim(
      x0,
      function(x) {
        return(-llik1(yW[ , i1], mu[i1], x[1], x[2], mu[i1])$llik)
      },
      method = "BFGS",
      control = list(maxit = 1000)
    )
    llik1(yW[ , i1], mu[i1], estim$par[1], estim$par[2], mu[i1])$ft
  }
  Gas_ft_1[ , iter0] = Gas_ft[ , 1]
  Gas_ft_N[ , iter0] = Gas_ft[ , dim_N]
  toc()
  
  
  cat(iter0, " ... Univariate fits (biased targeting; subsequent debiasing)\n")
  tic()
  GasSS_ft = foreach(i1=1:dim_N , .combine = cbind) %dopar% {
    x0 = c(-3,3)
    estim = optim(
      x0,
      function(x) {
        return(-llik1(yW[ , i1], lambda[i1], x[1], x[2], lambda[i1])$llik)
      },
      method = "BFGS",
      control = list(maxit = 1000)
    )
    llik1(yW[ , i1], lambda[i1], estim$par[1], estim$par[2], lambda[i1])$ft
  }
  GasSS_ft_1[ , iter0] = GasSS_ft[ , 1]
  GasSS_ft_N[ , iter0] = GasSS_ft[ , dim_N]
  
  GasS_ft = foreach(i1=1:dim_T, .combine = rbind) %dopar% {
    shrink4(GasSS_ft[i1, ], dim_N/dim_T)
  }
  GasS_ft_1[ , iter0] = GasS_ft[ , 1]
  GasS_ft_N[ , iter0] = GasS_ft[ , dim_N]
  
  ## intermediate sign of life and plot
  Gas_ft_1_summ = summarize_replications(Gas_ft_1, upto = iter0)
  Gas_ft_N_summ = summarize_replications(Gas_ft_N, upto = iter0)
  GasS_ft_1_summ = summarize_replications(GasS_ft_1, upto = iter0)
  GasS_ft_N_summ = summarize_replications(GasS_ft_N, upto = iter0)
  GasSS_ft_1_summ = summarize_replications(GasSS_ft_1, upto = iter0)
  GasSS_ft_N_summ = summarize_replications(GasSS_ft_N, upto = iter0)
  gg1 = make_frame(dim_snap, estimated_spectra, 
                     Gas_ft_1_summ, Gas_ft_N_summ, 
                     GasS_ft_1_summ, GasS_ft_N_summ, 
                     meanidx = meanidx, withplot = TRUE)
}


save(Gas_ft_1, Gas_ft_1_summ, Gas_ft_N, Gas_ft_N_summ, 
     GasS_ft_1, GasS_ft_1_summ, GasS_ft_N, GasS_ft_N_summ, 
     GasSS_ft_1, GasSS_ft_1_summ, GasSS_ft_N, GasSS_ft_N_summ, 
     file = paste0("T", dim_T, "_N", dim_N, "_rep", nr_replications, "_",
                   stri_replace_all(date(), ".", regex = ":")))

gg1 = make_frame(dim_snap, estimated_spectra, 
                   Gas_ft_1_summ, Gas_ft_N_summ, 
                   GasS_ft_1_summ, GasS_ft_N_summ, 
                   meanidx = 7, withplot = TRUE, withband = TRUE,
                   whichone = "max", ymax = 50)
gg2 = make_frame(dim_snap, estimated_spectra, 
                 Gas_ft_1_summ, Gas_ft_N_summ, 
                 GasS_ft_1_summ, GasS_ft_N_summ, 
                 meanidx = 7, withplot = TRUE, withband = TRUE,
                 whichone = "min", ymax = 5)
ggarrange(gg1, gg2, ncol = 2, common.legend = TRUE)

