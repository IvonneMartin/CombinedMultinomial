###############################################################################
##  DataGen() – simulasi data multivariate count + random-effect
###############################################################################
DataGen <- function(N,                         # banyak subjek
                    output_levels,             # scalar (Q)
                    input_levels,              # vector (p input)
                    beta,                      # koefisien tanpa intercept
                    S       = 2000,            # total count per sampel
                    theta   = 0.2,             # over-dispersion
                    s.u     = 0.8,             # sd random intercept
                    method  = c("DMM", "UNBM"),
                    T       = 2,               # time-point per subjek
                    baseline_out = 1L,
                    baseline_in  = rep(1L, length(input_levels))) {

  method <- match.arg(method)

  ## ------------------------------------------------------------------------
  ## 1)  Design matrix & koefisien
  ## ------------------------------------------------------------------------
  Des   <- DesMatIO(output_levels,
                    input_levels,
                    baseline_out,
                    baseline_in)

  Dmat  <- Des$Des.mat                     # intercept + dummy + interaksi
  Q     <- output_levels
  p_beta<- ncol(Dmat) - 1L                 # tanpa intercept

  if (length(beta) != p_beta)
    stop("length(beta) must be ", p_beta,
         " to match columns of design matrix (excluding intercept).")

  b_vec <- c(0, beta)                      # prepend intercept = 0 (log-scale)

  ## Helper untuk fixed-effect per baris
  row_fixef <- function(x_row)
    FixEf(x_row, Q, as.list(input_levels), b_vec, Des)

  ## ------------------------------------------------------------------------
  ## 2)  Generate covariate levels & random effect
  ## ------------------------------------------------------------------------
  n_in   <- length(input_levels)
  n_time <- T
  Ntot   <- N * n_time

  ## (a) matriks level kovariat (Ntot × n_in)
  Xmat <- matrix(nrow = Ntot, ncol = n_in)
  for (j in seq_len(n_in))
    Xmat[, j] <- sample(seq_len(input_levels[j]), Ntot, TRUE)

  ## (b) random intercept individu
  Usubj <- rnorm(N, 0, s.u)
  Urep  <- rep(Usubj, each = n_time)       # panjang N total

  ## ------------------------------------------------------------------------
  ## 3)  Linear predictor (log-scale)  =>  LP
  ## ------------------------------------------------------------------------
  LP <- t(apply(Xmat, 1, row_fixef))       # Ntot × Q
  LP <- LP + matrix(Urep, Ntot, Q)         # tambahkan efek acak

  ## ------------------------------------------------------------------------
  ## 4)  Simulasi count
  ## ------------------------------------------------------------------------
  if (method == "DMM") {                   # Dirichlet-Multinomial
    alpha <- (1/theta) * exp(LP)
    counts<- t(apply(alpha, 1, HMP::Dirichlet.multinomial, Nrs = S))
  } else {                                 # UNBM  (NegBin “tak terbatas”)
    size   <- (1/theta) * exp(log(S/10) + LP)
    probNB <- (1/theta) / (1 + 1/theta)
    counts <- matrix(nrow = Ntot, ncol = Q)
    for (i in seq_len(Ntot))
      counts[i, ] <- rnbinom(Q, size = size[i, ], prob = probNB)
  }

  ## ------------------------------------------------------------------------
  ## 5)  Susun data long
  ## ------------------------------------------------------------------------
  out <- data.frame(
    ID   = rep(seq_len(N), each = n_time),
    time = rep(seq_len(n_time),  N),
    counts,
    Xmat
  )
  names(out)[3:(2+Q)]             <- paste0("C", seq_len(Q))
  names(out)[(3+Q):ncol(out)]     <- paste0("X", seq_len(n_in))

  return(out)
}

