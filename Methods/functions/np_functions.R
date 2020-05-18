# author: Roman, Santiago, Mariana
# date: 18 May 2020

#==========================================================================
# Functions and Routines, Non Parametric 
#==========================================================================

# funcion para calcular el índice de vulnerabilidad usando PCA
scores.PCA <- function(X,ind){
  #parse matrix 
  X <- as.matrix(X)
  
  Xs <- X[ind,] 
  R <- cor(Xs)
  z <- ( eigen(R)$vectors[,1] %>% as.vector())
  index <- X %*% z
  
  return(abs(100/index))
}

# función para obtener el remuestreo Bootstrap usando boot()
boot.PCA <- function(X, B, seed){
  # set seed
  set.seed(seed)
  
  #bootsrap
  boot.scores <- boot(X, statistic = scores.PCA, R = B)
  
  return(boot.scores)
}

# función para obtener el remuestreo Bootstrap a mano
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
  
  return(abs(scores)) # Se toma el valor absoluto (Explicación)
}








