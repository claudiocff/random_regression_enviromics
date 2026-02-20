#' Orthogonal/Orthonormal Legendre Polynomial Basis on [-1, 1]
#'
#' Constructs a Legendre polynomial basis evaluated at a numeric vector scaled to the interval
#' [-1, 1], up to a specified order. Optionally returns the orthonormal basis via the standard
#' normalization constant.
#'
#' Details:
#' - The basis is generated using the three-term recurrence for Legendre polynomials:
#'   #$# P_{0}(x) = 1, \quad P_{1}(x) = x #$#
#'   #$# P_{k+1}(x) = \frac{(2k+1)\,x\,P_{k}(x) - k\,P_{k-1}(x)}{k+1} \quad \text{for } k \ge 1 #$#
#' - Orthonormal scaling uses:
#'   #$# \tilde{P}_{n}(x) = \sqrt{\frac{2n+1}{2}}\, P_{n}(x) #$#
#'   so that #$# \int_{-1}^{1} \tilde{P}_{n}(x)\, \tilde{P}_{m}(x)\, dx = \delta_{nm} #$#.
#'
#' @param x Numeric vector (already scaled to [-1, 1]) where the basis is evaluated.
#' @param order Integer. Maximum polynomial order (degree) to compute. Returns columns 0..order.
#' @param orthonormal Logical. If TRUE, applies orthonormal scaling to each polynomial degree.
#'
#' @return
#' A numeric matrix of size length(x) x (order+1), where column j corresponds to #$# P_{j}(x) #$#
#' for j = 0..order. If `orthonormal = TRUE`, columns are the orthonormally scaled polynomials.
#'
#' @examples
#' x <- seq(-1, 1, length.out = 5)
#' B <- legendre_basis(x, order = 3, orthonormal = FALSE)  # orthogonal basis
#' B_ortho <- legendre_basis(x, order = 3, orthonormal = TRUE)  # orthonormal basis
#'
#' @export
legendre_basis <- function(x, order, orthonormal = FALSE) {
  # x: vector ya escalado a [-1,1]
  n <- length(x)
  # Pre-allocate the basis matrix with columns for degrees 0..order
  B <- matrix(0, nrow = n, ncol = order + 1)
  # P0(x) = 1
  B[,1] <- 1
  # P1(x) = x
  if (order >= 1) B[,2] <- x
  # Use three-term recurrence to build up to the desired order
  if (order >= 2) {
    for (k in 1:(order-1)) {
      # Recurrence:
      # P_{k+1}(x) = ((2k+1) * x * P_k(x) - k * P_{k-1}(x)) / (k+1)
      B[,k+2] <- ((2*k+1) * x * B[,k+1] - k * B[,k]) / (k+1)
    }
  }
  if (orthonormal) {
    # Escala columnas (excepto P0) para hacer ortonormal
    # tilde(P_n) = sqrt((2n+1)/2) * P_n
    for (d in 0:order) {
      # Apply the normalization factor to each degree
      B[,d+1] <- sqrt((2*d+1)/2) * B[,d+1]
    }
  }
  # Return basis matrix
  B
}

#' Bayesian Multi-Environment Trial Analysis with BGLR via Random Regression (Finlay–Wilkinson)
#'
#' Fits the Finlay–Wilkinson model in two steps using the BGLR package, with genotype
#' effects through a genomic relationship matrix, environment main effects, and a random-regression
#' slope per genotype based on the environment index estimated in Step 1.
#'
#' Conceptually, the model can be written as:
#' #$# y_{ge} = \mu + l_e + \alpha_g + \beta_g \, l_e + \varepsilon_{ge} #$#
#' where y is the phenotype for genotype g in environment e, l_e is the environment index,
#' alpha_g is the genotype-specific intercept, beta_g is the genotype-specific slope
#' (sensitivity to environment), and epsilon is the residual.
#'
#' @param data A data frame containing the multi-environment trial data.
#' @param gen Character. Column name for genotype IDs.
#' @param env Character. Column name for environment IDs.
#' @param trial Character (optional). Column name for experiment/trial factor (if multiple trials within environments).
#' @param block Character (optional). Column name for replication/block factor.
#' @param check Character (optional). Column name with logical values (TRUE/FALSE) marking check genotypes.
#' @param trait Character. Column name for the response variable (phenotype).
#' @param res.het Logical. Heterogeneous residuals? (Present for compatibility; not used directly). Default = FALSE.
#' @param nIter Integer. Number of MCMC iterations for BGLR. Default = 2000.
#' @param burnIn Integer. Burn-in period. Default = 500.
#' @param verbose Logical. Print progress from BGLR? Default = TRUE.
#' @param saveAt Character. Path/prefix to save binary files from BGLR (used later by `readBinMat`). Default = "" (current working directory).
#' @param kinship.matrix Matrix (required). Genomic relationship matrix with row/column names matching genotype IDs in `data[[gen]]`.
#' @param seed Integer. Random seed for reproducibility. Default = 123.
#'
#' @details
#' Workflow:
#' - Step 1: Fit genotype (G) and environment (E) main effects via BGLR; obtain the environment projections to build an environmental index.
#' - Step 2: Fit genotype-specific intercept and slope (G + G×E) using the index derived from Step 1 through a random-regression formulation.
#'
#' Implementation notes:
#' - The genomic relationship matrix is aligned internally to the genotype incidence matrix (Zg) order, ensuring consistency between model inputs.
#' - The environment design matrix is centered to represent the environmental index used in the slope computation.
#' - Posterior samples for intercept and slope are read from the binary files saved by BGLR (`ETA_int_b.bin`, `ETA_slope_b.bin`), which depend on `saveAt`.
#' - BLAS threading is controlled via RhpcBLASctl for reproducibility/performance.
#'
#' Checks performed by the code:
#' - Verifies that specified columns (`gen`, `env`, `trait`, and optionally `trial`, `block`, `check`) exist in `data`.
#' - Filters `data` to genotypes present in `rownames(kinship.matrix)`.
#' - Aligns `kinship.matrix` to match the genotype incidence matrix (`Zg`) ordering/content.
#'
#' @return
#' An object of class `bglr_met` containing:
#' - `path`: the `saveAt` path used to write/read BGLR binaries.
#' - `GGE`: wide matrix (by genotype:environment rows) with raw predictions (G + GE) across iterations.
#' - `GE`: wide matrix (by genotype:environment rows) with centered predictions (mu + XB) across iterations.
#' - `INT`: matrix of genotype intercepts by iterations.
#' - `SLOPE`: matrix of genotype slopes by iterations.
#' - `fit1`: BGLR fit object for Step 1 (G + E).
#' - `fit2`: BGLR fit object for Step 2 (G + GE).
#' - `data`: curated/aligned data used in the model.
#' - `trait`, `gen`, `env`, `trial`, `block`: column names used in the analysis.
#' - `check`: vector of check genotype IDs (if `check` column was provided).
#'
#' @examples
#' \donttest{
#' # Example using EnvRtype data
#' library(EnvRtype)
#' data("maizeG")     # genomic relationship matrix
#' data("maizeYield") # phenotypic data
#'
#' data <- maizeYield
#' kinship.matrix <- maizeG
#'
#' # Select 5 random genotypes as checks (as in your setup)
#' set.seed(123)
#' check_ids <- as.character(sample(unique(data$gid), 5))
#' data$check <- data$gid %in% check_ids
#'
#' # Identify a candidate trait column (exclude ID-like columns)
#' candidate_traits <- setdiff(names(data), c("gid","env","trial","block","check"))
#' trait_col <- candidate_traits[1] # use the first available trait column
#'
#' # Create a temporary output directory for BGLR binary files
#' out_dir <- tempfile(pattern = "bglr_fw_")
#' dir.create(out_dir)
#'
#' # Run the Finlay–Wilkinson two-step model
#' fit_fw <- FW_bglr(
#'   data = data,
#'   gen = "gid",
#'   env = "env",
#'   trait = trait_col,
#'   check = "check",
#'   kinship.matrix = kinship.matrix,
#'   nIter = 2000,
#'   burnIn = 500,
#'   saveAt = paste0(out_dir, "/"),
#'   verbose = TRUE
#' )
#'
#' # Inspect results
#' names(fit_fw)
#' dim(fit_fw$INT)   # intercepts by genotype x iterations
#' dim(fit_fw$SLOPE) # slopes by genotype x iterations
#' dim(fit_fw$GE)    # genotype:environment rows x iterations
#' dim(fit_fw$GGE)   # same, uncentered
#'
#' # Check genotypes
#' fit_fw$check
#' }
#'
#' @references
#' Finlay, K. W., & Wilkinson, G. N. (1963). The analysis of adaptation in a plant-breeding programme.
#' Australian Journal of Agricultural Research, 14(6), 742–754.
#'
#' @export
#' @importFrom BGLR BGLR Multitrait readBinMat
#' @importFrom tidyr expand_grid pivot_wider
#' @importFrom dplyr group_by mutate ungroup left_join select
#' @importFrom stats model.matrix
#' @importFrom reshape2 melt
#' @importFrom tibble column_to_rownames
FW_bglr <- function(data, # data frame containing the data
                    gen, # column name for genotypes
                    env, # column name for environments
                    trial = NULL, # optional, column name for trials
                    block = NULL, # optional, column name for replications
                    check = NULL, # optional, column containing TRUE or FALSE for checks
                    trait, # column name for the response variable
                    res.het = FALSE, # boolean for heterogeneous residuals
                    nIter = 2000, # number of iterations for Bayesian sampling
                    burnIn = 500, # number of burn-in iterations
                    verbose = TRUE, # boolean to print progress
                    saveAt = "", # path to save results
                    kinship.matrix = NULL, # genomic relationship matrix (required for some models)
                    W = NULL, #matrix containing environmental covariates (env x ec)
                    order = 1, #polynomial order
                    seed = 123) { # number of factors for FA models (required for some models)
  
  # Fix random seed for reproducibility across runs
  set.seed(seed)
  # BLAS threading controls: verify and then set number of threads to 1
  RhpcBLASctl::blas_get_num_procs() ## verification
  RhpcBLASctl::blas_set_num_threads(1) ## limit to 12 cores
  RhpcBLASctl::blas_get_num_procs() ## verification
  # Ensure BGLR is available
  library(BGLR)
  # Drop unused factor levels to avoid misalignment issues later
  data <- droplevels(data)
  # Initial data curation: check required columns exist in data
  stopifnot(
    gen %in% colnames(data),
    env %in% colnames(data),
    trait %in% colnames(data),
    is.null(trial) || trial %in% colnames(data),
    is.null(block) || block %in% colnames(data),
    is.null(check) || check %in% colnames(data)
  )
  
  # Cache check genotype IDs (if provided)
  if(!is.null(check)){
    check_names <- unique(data[data[[check]]==TRUE,gen])
  }
  
  # Keep only relevant columns depending on presence of trial/block
  if(!is.null(trial) & !is.null(block)){
    data <- data[, c(gen, env, trial, block, trait)]
  } else if (is.null(trial) & !is.null(block)){
    data <- data[, c(gen, env,block, trait)]
  } else if (!is.null(trial) & is.null(block)){
    data <- data[, c(gen, env, trial, trait)]
  }else{
    data <- data[, c(gen, env, trait)]
  }
  
  # Order rows by environment then genotype for stable design matrices
  data <- data[order(data[[env]], data[[gen]]), ]
  
  # If a kinship (G) matrix is provided, filter data to genotypes present there
  if(!is.null(kinship.matrix)){
    
    data <- data[data[[gen]]%in%rownames(kinship.matrix),]
    
  }
  
  # Response vector
  y <- data[[trait]]
  n <- nrow(data)
  
  # List of model components (ETA) for BGLR Step 1
  ETA <- list()
  
  # Genotype main effect: either via genomic relationship (GBLUP) or simple BLUP
  if(!is.null(kinship.matrix)){
    cat("kinship matrix provided: Running GBLUP \n")
    # Align kinship to current set/order of genotypes
    kinship.matrix <- kinship.matrix[rownames(kinship.matrix)%in%data[[gen]],colnames(kinship.matrix)%in%data[[gen]]]
    
    # Incidence matrix for genotypes
    Zg <- model.matrix(~ -1 + factor(data[[gen]]))
    colnames(Zg) <- gsub("factor\\(data\\[\\[gen\\]\\]\\)", "", colnames(Zg))
    # Reorder kinship rows/cols to match Zg columns
    kinship.matrix <- kinship.matrix[match(colnames(Zg), rownames(kinship.matrix)), match(colnames(Zg), colnames(kinship.matrix))]
    stopifnot(all(rownames(kinship.matrix) == colnames(Zg)))
    
    # Spectral decomposition to build principal components of kinship
    evdG <- eigen(kinship.matrix)
    PC <- sweep(evdG$vector, STATS = sqrt(evdG$values), MARGIN = 2, FUN = "*")
    # Project Zg onto PC space (orthogonalized random effect design)
    Xg <- Zg %*% PC
    colnames(Xg) <- colnames(Zg)
    
    # Bayesian ridge regression on projected genotype effects
    ETA$g <- list(X = Xg, model = "BRR", saveEffects = TRUE)
    
  } else{
    
    cat("kinship matrix not provided: Running BLUP \n")
    # Basic genotype incidence matrix without genomic information
    Zg <- model.matrix(~ -1 + factor(data[[gen]]))
    colnames(Zg) <- gsub("factor\\(data\\[\\[gen\\]\\]\\)", "", colnames(Zg))
    Xg <- Zg
    ETA$g <- list(X = Xg, model = "BRR", saveEffects = TRUE)
    
  }
  
  
  # Environment main effect:
  # Either modeled via environmental covariates W (PCA-based index) or via centered environment dummies
  if(!is.null(W)){
    cat("W provided: Environmental gradient will be predicted using environmental covariates \n")
    # Environment incidence matrix
    Xl <- model.matrix(~ -1 + factor(data[[env]]))
    colnames(Xl) <- gsub("factor\\(data\\[\\[env\\]\\]\\)", "", colnames(Xl))
    # Align W to environment order in Xl
    W <- W[match(colnames(Xl),rownames(W)),]
    stopifnot(all(rownames(W) == colnames(Xl)))
    
    # Standardize W (mean 0, sd 1) and keep parameters for out-of-sample prediction
    W_std <- scale(W, center = TRUE, scale = TRUE)
    mu_W <- attr(W_std, "scaled:center")
    sigma_W <- attr(W_std, "scaled:scale")
    
    # PCA on standardized W; retain enough PCs to explain >= 95% variance (or minimum 3)
    pr <- prcomp(W_std, center = FALSE, scale. = FALSE)
    expl <- cumsum(pr$sdev^2) / sum(pr$sdev^2)
    k <- which(expl >= 0.95)[1]
    if (is.na(k) || k < 1) k <- min(3, ncol(W_std))
    S     <- pr$x[, 1:k, drop = FALSE]              # [env x k]
    S_std <- scale(S, center = TRUE, scale = TRUE)
    cS <- attr(S_std, "scaled:center")
    sS <- attr(S_std, "scaled:scale")
    # Design for environment effect based on standardized scores
    XW <- Xl %*% S_std
    
    # Bayesian ridge regression on environment scores to infer environment index
    ETA$l <- list(X = XW, model = "BRR", saveEffects = TRUE)
    
  }else{
    
    cat("W not provided: Environmental gradient will be only based on average performance \n")
    # Centered environment dummies (remove global mean from environment effect)
    Xl <- model.matrix(~ -1 + factor(data[[env]]))
    colnames(Xl) <- gsub("factor\\(data\\[\\[env\\]\\]\\)", "", colnames(Xl))
    Xl <- scale(Xl, center = TRUE, scale = FALSE)
    ETA$l <- list(X = Xl, model = "BRR", saveEffects = TRUE)
    
  }
  
  
  # Optional trial and/or block effects (nested within environments)
  if (!is.null(trial) & !is.null(block)) {
    ETA$trial <- list(X = scale(model.matrix(~ -1 + factor(paste(data[[env]], data[[trial]]))), center = TRUE, scale = FALSE),
                      model = "BRR", saveEffects = TRUE)
    
    ETA$block <- list(X = scale(model.matrix(~ -1 + factor(paste(data[[env]], data[[trial]], data[[block]]))), center = TRUE, scale = FALSE),
                      model = "BRR", saveEffects = TRUE)
    
  }else if(is.null(trial) & !is.null(block)){
    ETA$block <- list(X = scale(model.matrix(~ -1 + factor(paste(data[[env]], data[[block]]))), center = TRUE, scale = FALSE),
                      model = "BRR", saveEffects = TRUE)
    
  }else if(!is.null(trial) & is.null(block)){
    ETA$trial <- list(X = scale(model.matrix(~ -1 + factor(paste(data[[env]], data[[trial]]))),center = TRUE, scale = FALSE),
                      model = "BRR", saveEffects = TRUE)
  }
  
  # Step 1: fit G + E model to estimate environment index
  cat("Running First Step (G + E) \n")
  
  fit1 <- BGLR::BGLR(y = y,
                     ETA = ETA,
                     nIter = nIter,
                     burnIn = burnIn,
                     saveAt = saveAt,
                     verbose = verbose)
  
  
  
  # Extract environment index: either via W-based PCA scores or directly from environment effects
  if(!is.null(W)){
    
    beta_w <- fit1$ETA[["l"]]$b
    # lHat1 is the inferred environment gradient (index) per environment
    lHat1 <-as.vector(S_std %*% beta_w)
    names(lHat1) <- rownames(W)
    
  }else{
    beta_w <- NULL
    lHat1 <- fit1$ETA[["l"]]$b
    
  }
  
  
  ## Build environment Legendre basis (orthonormal) at env level
  # Scale lHat1 to [-1, 1] to evaluate Legendre polynomials robustly
  rng_env <- range(lHat1, na.rm = TRUE)
  if (diff(rng_env) > 0) {
    x_env <- 2 * (lHat1 - rng_env[1]) / diff(rng_env) - 1
  } else {
    x_env <- rep(0, length(lHat1))
  }
  P_env <- legendre_basis(x_env, order, orthonormal = TRUE)
  rownames(P_env) <- names(lHat1)
  
  # Fixed-effect part for environment basis (excluding degree 0)
  P_env_fix <- P_env[, 2:(order+1), drop = FALSE]          # [nEnv x D]
  X_env_fix <- Xl %*% P_env_fix 
  
  # Observed per-row degree values (mapped through Xl), to build per-observation slopes
  P_obs_deg <- lapply(seq_len(order), function(d) as.vector(Xl %*% P_env[, d+1]))
  names(P_obs_deg) <- paste0("deg", 1:order)
  
  
  ## Step 2 on raw y: genotype intercept + env fixed + random regression per degree
  # Base ETA2: genotype intercept (random) + environment fixed Legendre basis
  ETA2 <- list(
    int     = list(X = Xg, model = "BRR",   saveEffects = TRUE),
    env_fix = list(X = X_env_fix, model = "FIXED", saveEffects = TRUE)  # or model="BRR", R2=0.1 if you prefer penalization
  )
  
  # Add random-regression slope terms per Legendre degree
  for (d in seq_len(order)) {
    # Interaction design: Zg * degree value at each observation
    Xd <- sweep(Zg, 1L, P_obs_deg[[d]], "*")   # [nObs x nGen]
    # If GBLUP, project to PC space for each slope
    if (!is.null(kinship.matrix)) Xd <- Xd %*% PC
    ETA2[[paste0("slope", d)]] <- list(X = Xd, model = "BRR", saveEffects = TRUE)
  }
  
  # Optional trial/block effects again in Step 2
  if (!is.null(trial) & !is.null(block)) {
    ETA2$trial <- list(X = scale(model.matrix(~ -1 + factor(paste(data[[env]], data[[trial]]))), center = TRUE, scale = FALSE),
                       model = "BRR", saveEffects = TRUE)
    
    ETA2$block <- list(X = scale(model.matrix(~ -1 + factor(paste(data[[env]], data[[trial]], data[[block]]))), center = TRUE, scale = FALSE),
                       model = "BRR", saveEffects = TRUE)
    
  }else if(is.null(trial) & !is.null(block)){
    ETA2$block <- list(X = scale(model.matrix(~ -1 + factor(paste(data[[env]], data[[block]]))), center = TRUE, scale = FALSE),
                       model = "BRR", saveEffects = TRUE)
    
  }else if(!is.null(trial) & is.null(block)){
    ETA2$trial <- list(X = scale(model.matrix(~ -1 + factor(paste(data[[env]], data[[trial]]))),center = TRUE, scale = FALSE),
                       model = "BRR", saveEffects = TRUE)
  }
  
  # Step 2: fit G + GE model with random regression
  cat("Running Second Step (G + GE) \n")
  
  fit2 <- BGLR(y = y,
               ETA = ETA2,
               nIter = nIter,
               burnIn = burnIn,
               saveAt = saveAt,
               verbose = verbose)
  
  # Collect environment fixed-effect coefficients (mu + basis coefficients)
  env_fix_coef <- c(fit2$mu,fit2$ETA$env_fix$b)
  names(env_fix_coef) <- paste0("deg", 0:order)
  
  
  ## Read effects back from binaries
  # Intercept (int) effects across MCMC iterations
  BInt <- readBinMat(paste0(saveAt, "ETA_int_b.bin"))
  # Environment fixed Legendre coefficients per iteration
  fix_coef_chain <- read.table(paste0(saveAt, "ETA_env_fix_b.dat"))
  colnames(fix_coef_chain) <- paste0("deg", 1:order)
  # Posterior mu chain
  mu_chain <- scan(paste0(saveAt, "mu.dat"))
  
  
  
  # Map intercept coefficients back to genotype scale (PC space if GBLUP)
  if (!is.null(kinship.matrix)){
    INT <- tcrossprod(PC, BInt) + mu_chain
  } else{
    INT <- t(BInt) + mu_chain
  }                          
  
  ## Degree-specific coefficients
  # For each slope degree, read coefficients and add env fixed basis contribution
  COEF <- vector("list", order)
  for (d in seq_len(order)) {
    Bd <- readBinMat(paste0(saveAt, "ETA_", paste0("slope", d), "_b.bin")) + fix_coef_chain[,paste0("deg",d)]
    
    if(!is.null(kinship.matrix)){
      COEF[[d]] <- tcrossprod(PC, Bd)
    }else{
      COEF[[d]] <- t(Bd)
    }                        
    
    rownames(COEF[[d]]) <- colnames(Xg)
  }
  
  # Name rows of intercepts by genotype IDs
  rownames(INT) <- colnames(Xg)
  
  ## Build environment gradient basis with Legendre polynomials (degrees 0..order)
  rng_env <- range(lHat1, na.rm = TRUE)
  if (diff(rng_env) > 0) {
    x_env <- 2 * (lHat1 - rng_env[1]) / diff(rng_env) - 1
  } else {
    x_env <- rep(0, length(lHat1))
  }
  P_env <- legendre_basis(x_env, order, orthonormal = TRUE)
  rownames(P_env) <- names(lHat1)
  
  # X_gradient holds degree columns (deg0..degD) at environment level
  X_gradient <- P_env
  L <- lHat1
  rownames(X_gradient) <- names(L)
  colnames(X_gradient) <- paste0("deg", 0:order)
  
  # Ensure deg0 is the constant 1 (redundant but explicit)
  X_gradient[,"deg0"] <- 1
  gid <- rownames(INT)
  
  
  gid <- rownames(INT)
  
  ## Get G + GE effects (mu + XB), centered vs uncentered
  GE <- list()
  GGE <- list()
  nIter2 <- ncol(INT)
  
  for (i in seq_len(nIter2)) {
    ## Beta per genotype for iteration i: [int, coef1, coef2, ..., coefD]
    Beta_i <- cbind(INT[, i],
                    do.call(cbind, lapply(COEF, function(M) M[, i])))
    rownames(Beta_i) <- gid
    ## Centered GE: double-centering by env and by genotype
    GE_i <- scale(X_gradient %*% t(Beta_i), center = TRUE, scale = FALSE)
    GE_i <- scale(t(GE_i), center = TRUE, scale = FALSE)
    
    # Reshape to wide format with "gen:env" rows, one column per iteration
    GE_i <- reshape2::melt(GE_i) |>
      dplyr::mutate(comb = paste(Var1, Var2, sep = ":")) |>
      dplyr::select(-c(Var1, Var2)) |>
      dplyr::mutate(iter = i) |>
      tidyr::pivot_wider(names_from = iter, values_from = value) |>
      tibble::column_to_rownames("comb")
    GE[[i]] <- GE_i
    
    ## Uncentered G+GE (GGE): raw mu + XB
    GGE_i <- X_gradient %*% t(Beta_i)
    GGE_i <- reshape2::melt(GGE_i) |>
      dplyr::mutate(comb = paste(Var1, Var2, sep = ":")) |>
      dplyr::select(-c(Var1, Var2)) |>
      dplyr::mutate(iter = i) |>
      tidyr::pivot_wider(names_from = iter, values_from = value) |>
      tibble::column_to_rownames("comb")
    GGE[[i]] <- GGE_i
  }
  
  # Concatenate iterations into single wide matrices
  GE <- do.call("cbind", GE)
  GGE <- do.call("cbind", GGE)
  
  # Return a structured object with all components
  structure(
    list(
      path = saveAt,
      GGE = GGE,
      GE = GE,
      INT = INT,
      ## keep the degree-1 slope for backward compatibility; also return all degrees
      SLOPE = if (order >= 1) COEF[[1]] else NULL,
      COEF = COEF,
      fit1 = fit1,
      fit2 = fit2,
      data = data,
      trait = trait,
      gen = gen,
      env = env,
      trial = trial,
      block = block,
      check = check_names,
      beta_w = beta_w,
      obs = "obs",
      W = if(!is.null(W)) W else NULL,
      XW = if(!is.null(W)) XW else NULL,
      order = order,
      X_gradient = X_gradient,
      range = rng_env,
      pca = if (!is.null(W)) list(
        rotation    = pr$rotation[, 1:k, drop = FALSE],
        k           = k,
        mu_W        = mu_W,
        sigma_W     = sigma_W,
        score_center = cS,
        score_scale  = sS
      ) else NULL
    ),
    class = "bglr_met"
  )
}

#' Predict Environment Gradient (Index) for New Environments via PCA Projection
#'
#' Uses the PCA mapping and coefficients learned in `FW_bglr` (when environmental covariates `W`
#' were provided) to predict the environment gradient/index #$# l_e #$# for a new set of environments.
#'
#' Procedure:
#' 1) Standardize the new `Wp` covariate matrix with the original training means and sds.
#' 2) Project standardized `Wp` onto the retained PCA loadings (scores).
#' 3) Standardize the new scores using the training score center/scale.
#' 4) Multiply by the learned regression coefficients #$# \beta #$# from Step 1 to obtain #$# l_e #$#.
#'
#' @param model A fitted object from `FW_bglr` with `pca` and `beta_w` components populated
#'   (i.e., `W` was provided during training).
#' @param Wp Numeric matrix of environmental covariates for new environments (rows = environments,
#'   columns = same covariates used in training). Row names should be the environment IDs.
#'
#' @return
#' A named numeric vector of predicted environment gradient values #$# l_e #$#, with names taken from
#' `rownames(Wp)`.
#'
#' @examples
#' # Suppose fit_fw was obtained with W provided:
#' # new_W is a matrix with the same columns as fit_fw$W and rownames as environment IDs
#' # l_new <- predict_env_gradient_pca(fit_fw, new_W)
#'
#' @export
predict_env_gradient_pca <- function(model, Wp) {
  # Ensure the model has PCA info (only available when W was provided during training)
  stopifnot(!is.null(model$pca))
  rot <- model$pca$rotation
  k   <- model$pca$k
  muW <- model$pca$mu_W
  sdW <- model$pca$sigma_W
  cS  <- model$pca$score_center
  sS  <- model$pca$score_scale
  beta <- model$beta_w
  
  
  # Standardize using training means (muW) and sds (sdW) to ensure compatibility
  Wp_std <- scale(Wp, center = muW, scale = sdW)
  
  
  # Project standardized covariates onto retained PCA loadings
  S_new <- as.matrix(Wp_std) %*% rot   # [env_new x k]
  
  
  # Standardize new scores using the same center/scale used during training PCA scores
  S_new_std <- scale(S_new, center = cS, scale = sS)
  
  
  # Environment gradient prediction as linear combination of standardized scores
  l_pred <- as.vector(S_new_std %*% beta)
  names(l_pred) <- rownames(Wp)
  l_pred
}