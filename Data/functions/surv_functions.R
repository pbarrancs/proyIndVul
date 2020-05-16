# author: Roman Velez
# date: 7 May 2020

#==========================================================================
# Functions and Routines for survival analysis
#==========================================================================

#### ========================= FUNCTIONS ================================
conf_int_mean <- function(df_km, tao, alpha = 0.05){
  #quotient
  di_ni.di <- df_km[,'dj'] / (df_km[,'nj'] * (df_km[,'nj'] - df_km[,'dj']))
  di_ni.di[di_ni.di == Inf] <- 0
  
  #step integral
  ind_int <- indiv_intergal_St(tj = df_km[,'tj'], Sj = df_km[,'St'], tao = tao)
  
  # MEAN
  # mean up to tao
  mean_tao <- trunc_mean2(tj = df_km[,'tj'], Sj = df_km[,'St'], tao = tao)
  
  # VARAINCE
  # pop t0 = 0
  ind_int <- ind_int[-1]
  
  # cumm sum backwards
  int_ti.Tao <- ind_int %>% rev() %>% cumsum() %>% rev()
  
  # variance
  var_mu <- (di_ni.di * int_ti.Tao) %>% sum()
  
  #std dev
  std_dev <- sqrt(var_mu)
  
  # CONFIDENCE INTERVALS
  ci_lower  <- mean_tao - qnorm(1-alpha/2) * std_dev
  ci_upper <- mean_tao + qnorm(1-alpha/2) * std_dev
  
  # values
  ans <- c(mean_tao, var_mu, std_dev, ci_lower, ci_upper)
  names(ans) <- c('mean','var','std. dev', 'ci_lower', 'ci_upper')
  
  return(ans)
}

indiv_intergal_St <- function(tj, Sj, tao){
  #last index such that tj <= tao
  k <- length(which(tj <= tao))
  
  #individual integrals
  # times
  tj1 <- c(tj[1:k], tao)
  tj2 <- c(0,tj1[1:k])
  
  # surv
  sj <- c(1,Sj[1:k])
  
  # individual integrals
  ind_int <- (tj1 - tj2) * sj
  
  return(ind_int)
}

trunc_mean2 <- function(tj, Sj, tao){
  #individual areas for each pair of time
  ind_int <- indiv_intergal_St(tj, Sj, tao)
  
  #overal integrla
  area <- sum(ind_int)
  
  return(area)
}

trunc_mean <- function(tj, Sj, tao){
  #size of the array of the times
  n <- length(tj)
  
  #last index such that tj <= tao
  k <- length(which(tj <= tao))
  
  #case when Tao > tmax or Tao <= tmax
  if(k < n){
    # truncated data up to tao
    tj1 <- tj[1:k]
    Sj <- Sj[1:k]
    
    #intervals to make crude integral
    tj2 <- c(0, tj1[-k])
    sj2<-c(1, Sj[-k])
  } else {
    # extended data up to tao
    tj1 <- c(tj,tao)
    
    #intervals for crude integral
    tj2 <- c(0, tj)
    sj2<-c(1,Sj)
  }
  
  #numeric integration
  # overestimated the integral  
  mu <- sum((tj1-tj2)*sj2)
  
  return(mu)
}

deriv_Ht <- function(tj,Hj){
  #len
  n <- length(tj)
  ind1 <- 1:(n-1)
  ind2 <- 2:n
  
  #h1
  h1 <- Hj[1]
  
  #hj
  hj <- c(h1,Hj[ind2] - Hj[ind1])
  
  return(hj)
}

parse_censor  <- function(x){
  #parse
  x <- as.character(x)
  
  #length of string
  n <- nchar(x)
  
  #last char
  c <- substr(x,n,n)
  
  #censoring
  if(c == '+'){
    num <- as.integer(substr(x,1,n-1))
    c <- 0
  }
  else{
    num <- as.integer(x)
    c <- 1
  }
  
  return(c(num,c))
}

ind_risk <- function(N, dj, cj){
  # first risk
  nj <- c(N)
  
  #len
  n <- length(dj)
  # risks
  for(i in 1:(n-1)){
    nj <- c(nj,nj[i] - (dj[i] + cj[i]))
  }
  return(nj)
}


parameter_t0 <- function(df, t0, param = 'St'){
  #indicator
  indi <- df[,'tj'] <= t0
  
  #last tj
  index_tj <- length(which(indi == TRUE))
  
  #S(t0)
  paramt0 <- as.matrix(df[index_tj, param])
  
  return (paramt0)
}


quantile_est <- function(St, ee, tj, p = 0.5){
  # size
  n <- length(St)
  
  #find infimum
  index_max <- length(which(St > 1 - p))
  
  # estimator
  tp <- tj[index_max]
  return(tp)
}

quantile_ci <- function(St, ee, tj, p = 0.5, alpha = 0.05){
  #indicator 
  indi <- abs(St - (1 - p)) / ee <= qnorm(1-alpha/2)
  
  #interval of the condition
  tj_ci <- tj[indi]
  
  #sort and return intervals
  tj_ci <- sort(tj_ci)
  n <- length(tj_ci)
  
  #confidence
  ci_lower <- tj_ci[1]
  ci_upper <- tj_ci[n]
  
  return(c(ci_lower,ci_upper))
}

#### ========================= ROUTINES ===================================
censor_time <- function(df, column = 1){
  # times
  t <- df[,column]
  
  #df t and c
  df <- apply(X = t, MARGIN = 1, FUN = parse_censor)
  
  # transpose
  df <- t(df)
  
  #tidy
  colnames(df) <- c('time', 'censored')
  
  #to dataframe
  df <- data.frame(df)
  
  return(df)
}

table_KM <- function(df){
  #' Product Limit Estimator and its estimated variance and std deviations BY HAND
  #' according to Klein book
  
  #total events
  N <- dim(df)[1]
  
  # times 
  tj <- unique(df[, 'time'])
  
  # censored
  cj <- df %>% group_by(time) %>% summarise(cj = sum(censored == 0))
  cj <- cj[,'cj'] %>% as.matrix()
  
  # events
  dj <-df %>% group_by(time) %>% summarise(dj = sum(censored == 1))
  dj <- dj[,'dj'] %>% as.matrix()
  
  # nj 
  nj <- ind_risk(N = N, dj = dj, cj = cj)
  
  # dj/nj
  dj_nj <- dj/nj
  
  #Kaplan Meier
  St <- cumprod(1-dj_nj)
  
  #sum for variance
  sum_di.ni_ni.di <- cumsum(dj/(nj*(nj-dj)))
  
  #varaince
  var.St <- St*St*sum_di.ni_ni.di
  
  #std dev
  std_dev.St <- sqrt(var.St)
  
  #table
  df <- data.frame(tj, nj, dj, cj, dj_nj, St, sum_di.ni_ni.di, var.St, std_dev.St)
  
  # delete observations with ONLY censored data
  indi <- ifelse(dj == 0 & cj > 0,FALSE,TRUE)
  df <- df[indi,]
  
  return(df)
}

table_KM_surv <- function(df){
  #kaplan meier fit
  sfit <- survfit(Surv(time,censored) ~ 1, data = df)
  
  # basic table
  tj <- sfit$time
  nj = sfit$n.risk 
  dj = sfit$n.event 
  cj = sfit$n.censor
  
  # dj/nj
  dj_nj <- dj/nj
  
  #Kaplan Meier
  St <- cumprod(1-dj_nj)
  
  #sum for variance
  sum_di.ni_ni.di <- cumsum(dj/(nj*(nj-dj)))
  
  #varaince
  var.St <- St*St*sum_di.ni_ni.di
  
  #std dev
  std_dev.St <- sqrt(var.St)
  
  #table
  df <- data.frame(tj, nj, dj, cj, dj_nj, St, sum_di.ni_ni.di, var.St, std_dev.St)
  
  # delete observations with ONLY censored data
  indi <- ifelse(dj == 0 & cj > 0,FALSE,TRUE)
  df <- df[indi,]
  
  return(df)
}


table_NA <- function(df){
  #kaplan meier fit
  sfit <- survfit(Surv(time,censored) ~ 1, data = df)
  
  # basic table
  tj <- sfit$time
  nj = sfit$n.risk 
  dj = sfit$n.event 
  
  #Nelson Aalen
  Ht_na <- cumsum(dj/nj)
  
  #Variance N.A
  var.Ht <- cumsum(dj/(nj*nj))
  
  #std dev
  std.Ht <- sqrt(var.Ht)
  
  #dataframe
  df <- data.frame(tj, nj, dj, Ht_na, var.Ht, std.Ht)
  
  return(df)
}

table_confidence <- function(df, alpha = 0.05){
  #kaplan meier fit
  sfit <- survfit(Surv(time,censored) ~ 1, data = df)
  
  # basic table
  tj <- sfit$time
  nj = sfit$n.risk 
  dj = sfit$n.event 
  cj = sfit$n.censor 
  St = sfit$surv
  
  # extras
  sum_di.ni_ni.di <- cumsum(dj/(nj*(nj-dj)))
  
  #varaince
  var.St <- St*St*sum_di.ni_ni.di
  
  #std dev
  std_dev.St <- sqrt(var.St)
  
  #linear confidence interval
  ci_lower_linear <- St - qnorm(1-alpha/2) * std_dev.St
  ci_lower_linear <- ifelse(ci_lower_linear < 0, 0, ci_lower_linear) #min(0,ci)
  
  ci_upper_linear <- St + qnorm(1-alpha/2) * std_dev.St
  ci_upper_linear <- ifelse(ci_upper_linear > 1, 1, ci_upper_linear) #max(1,ci)
  
  #log-interval
  W <- exp(qnorm(1-alpha/2) * sqrt(var.St/(St*St)) / log(St))
  ci_lower_log <- St^(1/W)
  ci_lower_log <- ifelse(ci_lower_log < 0, 0, ci_lower_log) #min(0,ci)
  ci_upper_log <- St^(W)
  ci_upper_log <- ifelse(ci_upper_log > 1, 1, ci_upper_log) #max(1,ci)
  
  
  #dataframe
  df <- data.frame(tj, St, ci_lower_linear, ci_upper_linear, ci_lower_log, ci_upper_log)
  
  return(df)
}

plot_ci <- function(df_ci, linear = TRUE, log = TRUE){
  
  #subtitle
  if(linear & log){
    subttl <- 'Red: linear conf. int. Green: log conf. int.'
  }
  else{ 
    if(linear){subttl <- 'Red: linear conf. int.'}
    else 
      if(log){subttl <- 'Green: log conf. int.'}
    else{subttl <- ''}
  }
  
  #plot
  plot(x = df_ci$tj, y = df_ci$St, type = 'l', xlab = 'time', ylab = 'St', main = 'Survival Plot', sub = subttl, ylim = c(0,1))
  
  #ci 
  if(linear){
    lines(df_ci$tj, df_ci$ci_lower_linear, col = 'red')
    lines(df_ci$tj, df_ci$ci_upper_linear, col = 'red')
  }
  
  if(log){
    lines(df_ci$tj, df_ci$ci_lower_log, col = 'green')
    lines(df_ci$tj, df_ci$ci_upper_log, col = 'green')
  }
}

table_left_trunc <- function(t, u, c){
  #total events
  N <- length(t)
  
  #sort vecs by times
  u <- u[order(t)]
  c <- c[order(t)]
  t <- sort(t)
  
  
  #individual times tj
  tj <- t %>% unique() %>% sort()
  
  # censored
  cj <- data.frame(t, c) %>% 
    group_by(t) %>% 
    summarise(cj = sum(c == 0)) %>% 
    select(cj) %>% 
    as.matrix()
  
  # obs
  dj <- data.frame(t, c) %>% 
    group_by(t) %>% 
    summarise(dj = sum(c == 1)) %>% 
    select(dj) %>% 
    as.matrix()
  
  # risk n
  nj <- ind_risk(N = N, dj = dj, cj = cj)
  
  # individual times
  n <- length(tj)
  
  #nj'
  nj_prime <- c()
  for(i in 1:n)
    nj_prime <- c(nj_prime, left_trunc(t0 = tj[i], ui = u, ti = t, closed = FALSE))
  
  #nj''
  nj_dprime <- c()
  for(i in 1:n)
    nj_dprime <- c(nj_dprime, left_trunc(t0 = tj[i],ui = u, ti = t, closed = TRUE))
  
  #kaplan Meier estimations 
  St <- (1 - dj/nj) %>% cumprod()
  St_prime <- (1 - dj/nj_prime) %>% cumprod()
  St_dprime <- (1 - dj/nj_dprime) %>% cumprod()
  
  #dataframe
  df <- data.frame(tj, dj, cj, nj, nj_prime, nj_dprime, St, St_prime, St_dprime)
  
  return(df)
}

left_trunc <- function(t0, ui, ti, closed = TRUE){
  # t0 and int
  # u & c are vectors
  if(closed)
    r <-  length( which (ui <= t0 & t0 <= ti ))
  else
    r <-  length( which (ui < t0 & t0 <= ti ))    
  return(r)
}