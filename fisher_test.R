##################################################################################################
df <- read.table(
  "/home/bismark/4acC_A3A_seq_combine.Csite.txt",
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

# ========== 1. Fisher's exact test between biological replicates ==========

# CtoT rep1 vs rep2
df$pvalue_rep_CtoT <- mapply(function(m1, u1, m2, u2) {
  mat <- matrix(c(m1, u1, m2, u2), nrow = 2)
  if (sum(mat) < 20) return(NA)
  fisher.test(mat)$p.value
},
df$CtoT.rep1_X,
df$CtoT.rep1_N - df$CtoT.rep1_X,
df$CtoT.rep2_X,
df$CtoT.rep2_N - df$CtoT.rep2_X
)

# ac4CtoT rep1 vs rep2
df$pvalue_rep_ac4CtoT <- mapply(function(m1, u1, m2, u2) {
  mat <- matrix(c(m1, u1, m2, u2), nrow = 2)
  if (sum(mat) < 20) return(NA)
  fisher.test(mat)$p.value
},
df$ac4CtoT.rep1_X,
df$ac4CtoT.rep1_N - df$ac4CtoT.rep1_X,
df$ac4CtoT.rep2_X,
df$ac4CtoT.rep2_N - df$ac4CtoT.rep2_X
)

# FDR correction
df$adj_p_rep_CtoT    <- p.adjust(df$pvalue_rep_CtoT, method = "BH")
df$adj_p_rep_ac4CtoT <- p.adjust(df$pvalue_rep_ac4CtoT, method = "BH")

# Count consistent sites between replicates
consistent_CtoT    <- sum(!is.na(df$adj_p_rep_CtoT) & df$adj_p_rep_CtoT > 0.05, na.rm = TRUE)
consistent_ac4CtoT <- sum(!is.na(df$adj_p_rep_ac4CtoT) & df$adj_p_rep_ac4CtoT > 0.05, na.rm = TRUE)
total_valid_CtoT    <- sum(!is.na(df$adj_p_rep_CtoT))
total_valid_ac4CtoT <- sum(!is.na(df$adj_p_rep_ac4CtoT))

cat("CtoT consistent sites between replicates:", consistent_CtoT, "/", total_valid_CtoT, "\n")
cat("ac4CtoT consistent sites between replicates:", consistent_ac4CtoT, "/", total_valid_ac4CtoT, "\n")

# ========== 2. Fisher's exact test between treatment groups (ac4CtoT vs CtoT) ==========

# Merge reads from replicates in each group
m_CtoT    <- df$CtoT.rep1_X + df$CtoT.rep2_X
u_CtoT    <- (df$CtoT.rep1_N - df$CtoT.rep1_X) +
  (df$CtoT.rep2_N - df$CtoT.rep2_X)
m_ac4CtoT <- df$ac4CtoT.rep1_X + df$ac4CtoT.rep2_X
u_ac4CtoT <- (df$ac4CtoT.rep1_N - df$ac4CtoT.rep1_X) +
  (df$ac4CtoT.rep2_N - df$ac4CtoT.rep2_X)

df$pvalue_group <- mapply(function(m1, u1, m2, u2) {
  mat <- matrix(c(m1, u1, m2, u2), nrow = 2)
  if (any(mat == 0) || sum(mat) < 20) return(NA)
  fisher.test(mat)$p.value
}, m_CtoT, u_CtoT, m_ac4CtoT, u_ac4CtoT)

df$adj_p_group <- p.adjust(df$pvalue_group, method = "BH")

# Calculate conversion rates and differences
df$rate_CtoT    <- m_CtoT / (m_CtoT + u_CtoT)
df$rate_ac4CtoT <- m_ac4CtoT / (m_ac4CtoT + u_ac4CtoT)

# Significant sites (FDR < 0.05)
sig <- df[
  !is.na(df$adj_p_group) &
    df$adj_p_group < 0.05,]

cat("Significant sites between treatment groups:", nrow(sig), "\n")

# ========== 3. Export results ==========

write.table(df, "/home/bismark/4acC_A3A_seq_combine.Csite.fisher.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)