FixEf <- function(x.vec, Q, lvl.cov, b, Des) {
  Y.form <- Des[[1]]           # tabel kombinasi level (kolom 1 = level output)
  D      <- Des[[2]]           # design matrix
  gtbeta <- as.vector(D %*% b) # nilai (xb) utk semua kombinasi level

  # --- Pilih baris di Y.form yang cocok dengan kombinasi input x.vec ---
  input.cols <- 2:(1 + length(x.vec))  # kolom ke-2 dst adalah kovariat input

  # Cari baris di mana semua entri input cocok dengan x.vec
  match.idx <- which(apply(Y.form[, input.cols, drop = FALSE], 1, function(row)
    all(row == x.vec)
  ))

  # Cek apakah jumlah baris cocok sesuai jumlah level output Q
  if (length(match.idx) != Q) {
    warning(sprintf("Jumlah baris cocok (%d) ≠ Q (%d)", length(match.idx), Q))
  }

  return(gtbeta[match.idx])
}
