library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# ==== DATA ====
df <- tribble(
  ~Category,        ~INF,                        ~t,                        ~TRTxT,                     ~BhelmXInfXTRTXt,
  "Firmicutes",     "ref",                       "ref",                     "ref",                      "ref",
  "Actinobacteria", "−0.006(−0.218, 0.207)",     "0.050(−0.155, 0.256)",   "0.046(−0.235, 0.326)",     "0.326(−0.042, 0.694)",
  "Bacteroidetes",  "0.220(−0.056, 0.496)",      "−0.119(−0.395, 0.157)",  "−0.012(−0.381, 0.356)",    "−0.916(−1.573, −0.259)",
  "Proteobacteria", "0.171(−0.054, 0.396)",      "0.056(−0.161, 0.273)",   "0.035(−0.256, 0.326)",     "0.026(−0.376, 0.427)",
  "Unclassified",   "−0.024(−0.304, 0.257)",     "0.129(−0.149, 0.407)",   "−0.099(−0.476, 0.277)",    "−0.159(−0.727, 0.410)",
  "pooled",         "0.166(−0.158, 0.490)",      "0.195(−0.124, 0.515)",   "−0.030(−0.449, 0.388)",    "−0.180(−0.814, 0.454)"
)

# ==== PERSIAPAN DATA ====
df_long <- df %>%
  pivot_longer(
    cols = c(INF, t, TRTxT, BhelmXInfXTRTXt),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  filter(Value != "ref") %>%
  mutate(
    # Ambil nilai tengah dan CI
    Value = str_replace_all(Value, "−", "-"),
    est = as.numeric(str_extract(Value, "[-]?[0-9.]+")),
    lower = as.numeric(str_extract(Value, "(?<=\\()[^,]+")),
    upper = as.numeric(str_extract(Value, "(?<=, )[0-9.-]+")),
    Category = factor(Category,
                      levels = c("Actinobacteria", "Bacteroidetes",
                                 "Proteobacteria", "Unclassified", "pooled"))
  )

# ==== FUNGSI TAMBAH KURUNG DI CI ====
add_parentheses <- function(lower, upper) {
  sprintf("(%s – %s)", round(lower, 3), round(upper, 3))
}

# ==== PILIH VARIABEL: INF, t, atau BhelmXInfXTRTXt ====
var_to_plot <- "BhelmXInfXTRTXt"  # ubah ke "t" atau "BhelmXInfXTRTXt" untuk lainnya
df_plot <- df_long %>% filter(Variable == var_to_plot)

# ==== PLOT ====
ggplot(df_plot, aes(x = est, y = Category)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.25,
                 color = "#4e79a7", linewidth = 1.2) +
  geom_point(size = 3.2, color = "#4e79a7") +
  #geom_text(aes(label = add_parentheses(lower, upper)),
  #          hjust = -0.1, size = 3.3, color = "gray30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_y_discrete(limits = rev(levels(df_plot$Category))) +
  labs(
    x = "Log Odds Ratio (95% CI)",
    y = NULL,
    title = paste0("Effect of ", var_to_plot, " (Ref: Firmicutes)")
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(face = "italic"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  xlim(min(df_plot$lower) - 0.5, max(df_plot$upper) + 0.5)
