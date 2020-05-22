# author: Roman, Santiago, Mariana
# date: 18 May 2020

#==========================================================================
# Functions and Routines, Non Parametric 
#==========================================================================


# INDICE DE VULNERABILIDAD - PCA
scores.PCA <- function(X,ind){
  #parse matrix 
  X <- as.matrix(X)
  
  Xs <- X[ind,] 
  R <- cor(Xs)
  z <- ( eigen(R)$vectors[,1] %>% as.vector())
  index <- X %*% z
  
  return(abs(100/index))
}

# REMUESTREO BOOTSTRAP - boot()
boot.PCA <- function(X, B, seed){
  # set seed
  set.seed(seed)
  
  #bootsrap
  boot.scores <- boot(X, statistic = scores.PCA, R = B)
  
  return(boot.scores)
}

# REMUESTREO BOOTSTRAP - manual
bootAux.PCA <- function(X, B, seed){
  #seed
  set.seed(seed)
  
  # as matrix
  X <- as.matrix(X)
  
  # vars
  m <- nrow(X)
  
  # loadings
  scores <- matrix(nrow = B, ncol = m)
  
  # bootsrap
  for ( i in 1:B ){
    ind <- sample(1:m, size = m, replace = TRUE)
    Xs <- X[ind,] 
    R <- cor(Xs)
    z <- ( eigen(R)$vectors[,1] %>% as.vector())
    scores[i, ] <- X %*% z
  }
  
  return(abs(100/scores)) # Se toma el valor absoluto (Explicacion)
}

# IC BOOTSTRAP - basic,norm,perc,bca
boot.IC <- function(X){
  # métodos
  types <- c("basic","norm","perc","bca")
  
  # Crea vector de nombres (method_L,method_U)
  ics <- c()
  for (n in types) {
    ics <- c(ics,paste0(n,"_L"),paste0(n,"_U"))
  }
  
  # Matriz para guardar ICs
  mBootCI <- matrix(NA,16,2*length(types), dimnames = list(delegaciones,ics))
  for (i in 1:16) {
    # boot.ci con "mBoot"
    BootCI <- boot.ci(mBoot, index = i, type=types, conf = .90)
    # Get CI
    mBootCI[i,] <- c(BootCI$basic[4:5], 
                       BootCI$normal[2:3], 
                       BootCI$percent[4:5],
                       BootCI$bca[4:5])
  }
  return(mBootCI)
}

# IC BOOTSTRAP - studentized a mano
boot.IC.t <- function(X, B, seed){
  #seed
  set.seed(seed)
  # as matrix
  X <- as.matrix(X)
  # vars
  m <- nrow(X)
  
  # store theta_hat_b
  scores <- matrix(nrow = B, ncol = m)
  
  # store se(theta_hat_b)
  se_b <- matrix(nrow = B,ncol = m)
  
  # bootsrap
  i <- 1
  for ( i in 1:B ){
    # muestra bootstrap (elige 16 delegaciones)
    ind <- sample(1:m, size = m, replace = TRUE) 
    Xs <- X[ind,]
    
    # compute and store index (theta_hat_b para cada del)
    R <- cor(Xs)
    z <- ( eigen(R)$vectors[,1] %>% as.vector())
    scores[i, ] <- X %*% z
    
    # compute standard error for each theta_hat_b: second bootstrap
    boot.scores <- boot(Xs, statistic = scores.PCA, R = B)
     # matriz de error estandar para cada theta_hat_b
    for (j in 1:16) {
      se_b[i,j] <- sd(boot.scores$t[,j])
    }
  }
  theta_hat <- 100/abs(scores)
  theta_hat <- apply(scores, 2, mean)
}




