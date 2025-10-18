#######

## =============================================================================
##  DesMat():  Membuat matriks dummy + interaksi untuk faktor-faktor kategorik
## -----------------------------------------------------------------------------
##  ARGUMEN
##  ---------------------------------------------------------------------------
##    levels   : integer vector, panjang p, banyaknya level tiap faktor.
##               Elemen pertama biasanya dianggap variabel dependen, sisanya
##               variabel independen – tetapi fungsi TIDAK bergantung pada
##               urutan; itu hanya konvensi pemanggil.
##
##    baseline : integer vector, panjang p, level referensi (default = 1L
##               untuk semua faktor). Kolom dummy hanya dibuat untuk level
##               ≠ baseline, demikian pula interaksinya.
##
##  NILAI KEMBALI
##  ---------------------------------------------------------------------------
##    list dengan dua komponen:
##      $LHS      : data.frame semua kombinasi level setiap faktor.
##      $Des.mat  : matriks (mode = integer) berisi
##                  - kolom pertama  : intercept  ("l_0")
##                  - kolom berikut : dummy main-effect & interaksi, nama kolom
##                    mengikuti pola <hurufFactor>_<level> atau
##                    <hurufHuruf>_<levelLevel> dst.
## =============================================================================
DesMat <- function(levels, baseline = rep(1L, length(levels))) {

  ## -- validasi ---------------------------------------------------------------
  stopifnot(length(levels) == length(baseline),
            all(levels    >= 1L),
            all(baseline  >= 1L),
            all(baseline  <=  levels))

  levels   <- as.integer(levels)
  baseline <- as.integer(baseline)
  p        <- length(levels)

  ## huruf penanda faktor: a, b, c, … (ganti jika p > 26)
  var_names <- letters[seq_len(p)]

  ## -- 1. semua kombinasi level (tabel LHS) ----------------------------------
  lhs <- expand.grid(lapply(levels, seq_len), KEEP.OUT.ATTRS = FALSE)
  names(lhs) <- var_names           # beri nama kolom a, b, c, …

  ## -- 2. mulai Design Matrix dengan intercept -------------------------------
  X          <- matrix(1L, nrow(lhs), 1L)
  colnames(X) <- "l_0"

  ## helper pembuat nama kolom
  make_colname <- function(vars, lvls)
    paste0(paste(vars, collapse = ""), "_", paste(lvls, collapse = ""))

  ## -- 3. loop main-effect & semua interaksi ---------------------------------
  for (k in seq_len(p)) {                        # k = 1 (main) ... p-way
    for (idx in combn(p, k, simplify = FALSE)) {

      ## level non-baseline tiap faktor terpilih
      lvl_lists <- Map(function(i) setdiff(seq_len(levels[i]), baseline[i]), idx)
      if (any(lengths(lvl_lists) == 0L)) next    # skip jika ada faktor tanpa level ≠ baseline

      ## semua kombinasi level non-baseline
      combos <- expand.grid(lvl_lists, KEEP.OUT.ATTRS = FALSE)

      for (r in seq_len(nrow(combos))) {
        lvls <- as.integer(combos[r, ])

        ## baris mana di LHS yang cocok dg kombinasi lvls?
        mask <- rep(TRUE, nrow(lhs))
        for (j in seq_along(idx)) {
          v    <- var_names[idx[j]]
          mask <- mask & (lhs[[v]] == lvls[j])
        }

        X <- cbind(X, as.integer(mask))
        colnames(X)[ncol(X)] <- make_colname(var_names[idx], lvls)
      }
    }
  }

  list(LHS = lhs, Des.mat = X)
}

