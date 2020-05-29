# author: Roman, Santiago, Mariana
# date: 18 May 2020

#==========================================================================
# Functions and Routines, Non Parametric Statistics
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
boot.IC <- function(data,c = 0.95){
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
    BootCI <- boot.ci(data, index = i, type=types, conf = c)
    # Get CI
    mBootCI[i,] <- c(BootCI$basic[4:5], 
                       BootCI$normal[2:3], 
                       BootCI$percent[4:5],
                       BootCI$bca[4:5])
  }
  return(mBootCI)
}

# IC BOOTSTRAP - stud a mano
boot.IC.t <- function(data, B, seed, c = .95){
  #seed
  set.seed(seed)
  # as matrix
  X <- as.matrix(data)
  # vars
  m <- nrow(X)
  # alpha
  alpha <- 1-c
  
  # store theta_hat_b
  scores <- matrix(nrow = B, ncol = m)
  
  # store se(theta_hat_b)
  se_b <- matrix(nrow = B,ncol = m)
  
  # Bootsrap
  for ( i in 1:B ){
    # muestra bootstrap (elige 16 delegaciones)
    ind <- sample(1:m, size = m, replace = TRUE) 
    Xs <- X[ind,]
    
    # compute and store index (theta_hat_b para cada del)
    R <- cor(Xs)
    z <- (eigen(R)$vectors[,1] %>% as.vector())
    scores[i, ] <- X %*% z
    
    # compute standard error for each theta_hat_b: second bootstrap
    boot.scores <- boot(Xs, statistic = scores.PCA, R = B)
    # matriz de error estandar para cada theta_hat_b
    for (j in 1:16) {
      se_b[i,j] <- sd(boot.scores$t[,j])
    }
  }
  # Obtenemos estimadores bootstrap (200 para cada delegacion)
  scores <- 100/abs(scores) # ajuste del indice
  
  # Cuantiles
  t <- matrix(nrow = B, ncol = m)
  q <- matrix(nrow = m, ncol = 2, dimnames = list(delegaciones,c("t_1-(a/2)","t_(a/2)")))
  for (d in 1:m) {
    # Estadistica t_b's para cada delegacion
    t[,d] <- (scores[,d] - indice[d])/se_b[,d] # (theta_hat - theta)/se(theta_hat)
    # obtener cuantiles al 90% de confianza
    q[d,] <- quantile(t[,d],probs = c(1-alpha/2,alpha/2))
  }
  
  # Intervalos
  IC_t <- matrix(nrow = m, ncol = 2)
  for (d in 1:m) {
    IC_t[d,] <- (indice[d] - q[d,]*sd(scores[,d]))
  }
  dimnames(IC_t) <- list(delegaciones,c("stud_L","stud_U"))
  return(IC_t)
}

# FORMATO ICS - solo sirve para este trabajo ok.
formato.ICS <- function(data){
  matriz <- data %>% as.matrix()
  basic <- c()
  norm <- c()
  perc <- c()
  bca <- c()
  stud <- c()
  for (d in 1:16) {
    basic <- c(basic,matriz[d,1:2])
    norm <- c(norm,matriz[d,3:4])
    perc <- c(perc,matriz[d,5:6])
    bca <- c(bca,matriz[d,7:8])
    stud <- c(stud,matriz[d,9:10])
  }
  m <- as.data.frame(rbind(basic,norm,perc,bca,stud))
  colnames(m) <- c("A.Obregon_L","A.Obregon_U",
                   "Azcapotzalco_L","Azcapotzalco_U",
                   "B.Juarez_L","B.Juarez_U",
                   "Coyoacan_L","Coyoacan_U",
                   "Cuajimalpa_L" ,"Cuajimalpa_U" , 
                   "Cuauhtemoc_L","Cuauhtemoc_U",
                   "G.Madero_L","G.Madero_U",
                   "Iztacalco_L","Iztacalco_U",
                   "Iztapalapa_L","Iztapalapa_U",
                   "M.Contreras_L","M.Contreras_U",
                   "M.Hidalgo_L","M.Hidalgo_U",
                   "MilpaAlta_L","MilpaAlta_U",
                   "Tlahuac_L","Tlahuac_U",
                   "Tlalpan_L","Tlalpan_U",
                   "V.Carranza_L","V.Carranza_U",
                   "Xochimilco_L","Xochimilco_U")
  return(m)
}

# JACKKIFE AFTER BOOTSTRAP
data <- df
B <- 200
seed
jackknife.after.boot <- function(data, B, seed){
  noms <- c("AO","AZ","BJ","CO","CM","CU","GAM","IZC","IZP","MC","MH","MA","TLH","TLP","VC","XO")
  #seed
  set.seed(seed)
  # as matrix
  X <- as.matrix(data)
  # vars
  m <- nrow(X)
  rownames(X) <- noms
  
  # loadings
  scores <- matrix(nrow = B, ncol = m)
  colnames(scores) <- noms
  indices <- matrix(nrow = B,ncol = m)
  colnames(indices) <- noms

  # Bootstrap
  for (b in 1:B ){
    # muestra bootstrap (elige 16 delegaciones)
    ind <- sample(1:m, size = m, replace = TRUE)
    Xs <- X[ind,]
    
    # J(i) indices de muestras sin cada una de las delegaciones
    for (d in 1:m) {
      if(!(noms[d] %in% rownames(Xs))) indices[b,d] <- b
    }
    
    # compute and store index (theta_hat_-j)
    R <- cor(Xs)
    z <- (eigen(R)$vectors[,1] %>% as.vector())
    scores[b, ] <- X %*% z
  }
  # Obtenemos estimadores bootstrap (200 para cada delegacion)
  scores <- 100/abs(scores) # ajuste del indice
  
  # Jackknife After Bootstrap: error estandar para cada muestra boot sin la delegacion "XX"
    se.B <- matrix(ncol = m, nrow = m) # error estandar para cada indice, para cada muestra sin del "XX"
    # se_jack_B(i)
    for (d in 1:m) { # Delegacion indice
      for (j in 1:m) { # Delegacion que no aparece
      # calculamos se para J(i)
        # Cada renglón es la del que NO aparece en el score para la delegacion de la columna
        s <- scores[!is.na(indices[,j]),d]
        se.B[j,d] <- sqrt(mean((s - mean(s))^2))
      }
    }
  
    # se_jack(se_jack_B(i))
    se <- function(x) sqrt((length(x)-1)*mean((x-mean(x))^2))
    
    se.jack <- apply(se.B, 2, se)
    
    return(se.jack)
}

