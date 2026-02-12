# Calcula a probabilidade de cada genótipo ser o máximo em cada ponto do gradiente
# model: objeto retornado por FW_bglr
# P_pred: vetor nomeado com os valores de gradiente (l) a avaliar
# tie.method: "share" divide a contagem igualmente entre empatados; "first" escolhe o primeiro
prob_best_by_gradient <- function(model, l_new, tie.method = c("share","first")) {
  tie.method <- match.arg(tie.method)
  
  # 1) Reescalar o gradiente para [-1, 1] usando o range do treino
  rng <- model$range
  if (diff(rng) > 0) {
    x_pred <- 2 * (l_new - rng[1]) / diff(rng) - 1
  } else {
    x_pred <- rep(0, length(l_new))
  }
  
  # 2) Base de Legendre ortonormal para os pontos de gradiente
  D <- model$order
  P_env_pred <- legendre_basis(x_pred, order = D, orthonormal = TRUE)
  rownames(P_env_pred) <- names(l_new)
  colnames(P_env_pred) <- paste0("deg", 0:D)
  # garantir deg0 = 1 (constante)
  P_env_pred[, "deg0"] <- 1
  
  # 3) Preparar as cadeias Beta por genótipo e iteração
  INT <- model$INT            # [gen x iter]
  COEF <- model$COEF          # lista de D matrizes [gen x iter]
  gid <- rownames(INT)
  nGen <- length(gid)
  nIter <- ncol(INT)
  nEnvNew <- nrow(P_env_pred)
  
  # 4) Inicializar matriz de probabilidades (gen x gradiente)
  prob_mat <- matrix(0, nrow = nGen, ncol = nEnvNew,
                     dimnames = list(gid, rownames(P_env_pred)))
  
  # 5) Loop sobre iterações para contar vencedores
  for (i in seq_len(nIter)) {
    # Beta_i: [gen x (1+D)] = [INT, COEF_1, COEF_2, ..., COEF_D]
    Beta_i <- cbind(INT[, i],
                    do.call(cbind, lapply(COEF, function(M) M[, i])))
    rownames(Beta_i) <- gid
    colnames(Beta_i) <- paste0("deg", 0:D)
    
    # Predições: Yhat_i = P_env_pred %*% t(Beta_i) => [envNew x gen]
    Yhat_i <- P_env_pred %*% t(Beta_i)
    
    # Para cada ponto do gradiente, identificar o(s) máximo(s)
    for (j in seq_len(nEnvNew)) {
      yj <- Yhat_i[j, ]
      mx <- max(yj, na.rm = TRUE)
      winners <- which(yj == mx)
      
      if (tie.method == "share") {
        # dividir igualmente a contagem em caso de empate
        prob_mat[winners, j] <- prob_mat[winners, j] + (1 / length(winners))
      } else {
        # atribui 1 ao primeiro vencedor
        prob_mat[winners[1], j] <- prob_mat[winners[1], j] + 1
      }
    }
  }
  
  # 6) Normalizar pela quantidade de iterações para obter probabilidades
  prob_mat <- prob_mat / nIter
  
  # Também podemos reportar o vencedor mais provável por gradiente
  winner <- apply(prob_mat, 2, function(p) names(p)[which.max(p)])
  list(
    prob = prob_mat,      # [gen x gradiente] com probabilidades
    winner = winner       # nome do genótipo mais provável por ponto de gradiente
  )
}


# Calcula produtividade média e desvio-padrão por genótipo em cada ponto do gradiente
# model: objeto retornado por FW_bglr
# l_new: vetor nomeado com valores do gradiente (ex.: saída de predict_env_gradient_pca)
# retorna: list(mean = [gen x grad], sd = [gen x grad])
prod_summary_by_gradient <- function(model, l_new) {
  # 1) Reescalar gradiente para [-1, 1] usando o range do treino
  rng <- model$range
  if (diff(rng) > 0) {
    x_pred <- 2 * (l_new - rng[1]) / diff(rng) - 1
  } else {
    x_pred <- rep(0, length(l_new))
  }
  
  # 2) Base de Legendre ortonormal
  D <- model$order
  P_env_pred <- legendre_basis(x_pred, order = D, orthonormal = TRUE)
  rownames(P_env_pred) <- names(l_new)
  colnames(P_env_pred) <- paste0("deg", 0:D)
  P_env_pred[, "deg0"] <- 1
  
  # 3) Preparar cadeias
  INT  <- model$INT      # [gen x iter]
  COEF <- model$COEF     # lista de D matrizes [gen x iter]
  gid  <- rownames(INT)
  nGen <- length(gid)
  nIter <- ncol(INT)
  nEnvNew <- nrow(P_env_pred)
  
  # 4) Acumuladores (Welford) para média e variância
  mean_mat <- matrix(0, nrow = nGen, ncol = nEnvNew,
                     dimnames = list(gid, rownames(P_env_pred)))
  m2_mat <- matrix(0, nrow = nGen, ncol = nEnvNew,
                   dimnames = list(gid, rownames(P_env_pred)))
  
  # 5) Loop das iterações: calcular Yhat_i e atualizar média/variância
  for (i in seq_len(nIter)) {
    Beta_i <- cbind(INT[, i],
                    do.call(cbind, lapply(COEF, function(M) M[, i])))
    rownames(Beta_i) <- gid
    colnames(Beta_i) <- paste0("deg", 0:D)
    
    # Yhat_i: [env x gen]
    Yhat_i <- P_env_pred %*% t(Beta_i)
    # Transpor para [gen x env] para combinar com nossos acumuladores
    Yi <- t(Yhat_i)
    
    # Atualização de Welford por célula
    # Para cada gen (linha) e grad (col), atualiza mean e m2
    delta <- Yi - mean_mat
    mean_mat <- mean_mat + delta / i
    m2_mat   <- m2_mat + delta * (Yi - mean_mat)
  }
  
  # 6) Desvio-padrão (amostral): sqrt(m2 / (nIter - 1))
  sd_mat <- sqrt(m2_mat / pmax(nIter - 1, 1))
  list(mean = mean_mat, sd = sd_mat)
}

















