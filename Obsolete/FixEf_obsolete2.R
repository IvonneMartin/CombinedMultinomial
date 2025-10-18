###############################################################################
##  FixEf() – menghitung efek-tetap untuk satu baris kovariat
###############################################################################
## ARGUMEN
##   x.vec   : vektor level kovariat untuk 1 observasi (panjang = #input)
##   Q       : banyak kategori respon (output_levels)
##   lvl.cov : list jumlah-level setiap input  (bisa: as.list(input_levels))
##   b       : vektor koefisien lengkap  (intercept + beta)
##   Des     : list hasil DesMat/DesMatIO  (Des[[1]] = LHS, Des[[2]] = Dmat)
## KELUARAN
##   vektor panjang Q  – linear-predictor (log-scale) untuk kategori 1..Q
###############################################################################
FixEf <- function(x.vec, Q, lvl.cov, b, Des) {

  Y.form <- Des[[1]]           # tabel kombinasi level
  D      <- Des[[2]]           # design matrix

  gtbeta <- as.vector(D %*% b) # nilai (xb) utk SEMUA kombinasi level

  ## --- pilih baris di LHS yg level-input-nya cocok dgn x.vec --------------
  idx    <- seq_len(nrow(Y.form))
  for (j in seq_along(x.vec)){            # kolom 2..(p) = input
    idx <- idx[ Y.form[idx, j + 1L] == x.vec[j] ]
  }
  ## baris yg cocok pasti tepat Q buah  (level output = 1..Q)
  return( gtbeta[idx] )
}
