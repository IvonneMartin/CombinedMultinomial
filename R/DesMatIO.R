## =============================================================================
##  (Opsional) Pembungkus DesMatIO(): memisahkan output vs input  --------------
## -----------------------------------------------------------------------------
##  ARGUMEN
##    output_levels  : integer tunggal – jumlah level variabel dependen.
##    input_levels   : integer vector – jumlah level tiap prediktor.
##    baseline_out   : baseline level untuk output            (default 1).
##    baseline_in    : baseline level untuk masing-masing input
##                     (default 1 untuk semuanya).
## =============================================================================
DesMatIO <- function(output_levels,
                     input_levels,
                     baseline_out = 1L,
                     baseline_in  = rep(1L, length(input_levels))) {

  levels   <- c(output_levels, input_levels)
  baseline <- c(baseline_out,  baseline_in)

  DesMat(levels = levels, baseline = baseline)
}

## =============================================================================
##  Contoh-contoh pemakaian ----------------------------------------------------
## =============================================================================

## -------------------------------------------------------------------
## Contoh 1  : 1 output (3 level); 2 input (masing-masing 3 level)
## -------------------------------------------------------------------
#ex1 <- DesMatIO(output_levels = 3,
#                input_levels  = c(3,3))

#str(ex1$LHS)      # lihat kombinasi level
#str(ex1$Des.mat)  # lihat matriks dummy
