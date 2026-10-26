##DSS
library(DSS)
require(bsseq)

dat.list <- list(
  CtoT_rep1    = read.table("CtoT_rep1_bismark_bt2_pe.deduplicated.CX_report.dss.txt", header = TRUE)[, c("chr", "pos", "N", "X")],
  CtoT_rep2    = read.table("CtoT_rep2_bismark_bt2_pe.deduplicated.CX_report.dss.txt", header = TRUE)[, c("chr", "pos", "N", "X")],
  ac4CtoT_rep1 = read.table("ac4CtoT_rep1_bismark_bt2_pe.deduplicated.CX_report.dss.txt", header = TRUE)[, c("chr", "pos", "N", "X")],
  ac4CtoT_rep2 = read.table("ac4CtoT_rep2_bismark_bt2_pe.deduplicated.CX_report.dss.txt", header = TRUE)[, c("chr", "pos", "N", "X")]
)

BSobj <- makeBSseqData(dat.list, c("CtoT_rep1", "CtoT_rep2", "ac4CtoT_rep1", "ac4CtoT_rep2"))

dml_test <- DMLtest(
  BSobj,
  group1 = c("CtoT_rep1", "CtoT_rep2"),
  group2 = c("ac4CtoT_rep1", "ac4CtoT_rep2"),
  smoothing = FALSE
)
dml_test_df = as.data.frame(dml_test)


write.table(dml_test_df, "F:/C12_N12_combine.txt", sep = "\t", quote = FALSE, row.names = FALSE)

dml_test = read.table("F:/C12_N12_combine.txt", sep = "\t", header = T)

library(dplyr)

df_CtoT_rep1    <- df_CtoT_rep1    %>% rename_with(~ paste0(.x, "_CtoT_rep1"),    c(N, X))
df_CtoT_rep2    <- df_CtoT_rep2    %>% rename_with(~ paste0(.x, "_CtoT_rep2"),    c(N, X))
df_ac4CtoT_rep1 <- df_ac4CtoT_rep1 %>% rename_with(~ paste0(.x, "_ac4CtoT_rep1"), c(N, X))
df_ac4CtoT_rep2 <- df_ac4CtoT_rep2 %>% rename_with(~ paste0(.x, "_ac4CtoT_rep2"), c(N, X))

library(dplyr)

df_motif = read.table("C12N12_motif.txt", sep = "\t", header = T)

df_merged <- df_CtoT_rep1 %>%
  full_join(df_CtoT_rep2,    by = c("chr", "pos", "strand", "Type", "X3.base")) %>%
  full_join(df_ac4CtoT_rep1, by = c("chr", "pos", "strand", "Type", "X3.base")) %>%
  full_join(df_ac4CtoT_rep2, by = c("chr", "pos", "strand", "Type", "X3.base")) %>%
  full_join(dml_test,        by = c("chr", "pos")) %>%
  tidyr::drop_na() %>%
  relocate(Type, X3.base, .after = strand) %>%
  filter(N_CtoT_rep1 > 50,
         N_CtoT_rep2 > 50,
         N_ac4CtoT_rep1 > 50,
         N_ac4CtoT_rep2 > 50)

df_merged_4acC <- df_merged %>%
  mutate(rep1 = ifelse(X_CtoT_rep1 / N_CtoT_rep1 > 0.10 & X_ac4CtoT_rep1 / N_ac4CtoT_rep1 < 0.05, "T", "F")) %>%
  mutate(rep2 = ifelse(X_CtoT_rep2 / N_CtoT_rep2 > 0.10 & X_ac4CtoT_rep2 / N_ac4CtoT_rep2 < 0.05, "T", "F")) %>%
  mutate(both = ifelse(X_CtoT_rep1 / N_CtoT_rep1 > 0.10 & X_CtoT_rep2 / N_CtoT_rep2 > 0.10
                       & X_ac4CtoT_rep1 / N_ac4CtoT_rep1 < 0.05 & X_ac4CtoT_rep2 / N_ac4CtoT_rep2 < 0.05, "T", "F"))

df_rep1_4acC <- df_merged_4acC[df_merged_4acC$rep1 == "T", ]
df_rep2_4acC <- df_merged_4acC[df_merged_4acC$rep2 == "T", ]
df_both_4acC <- df_merged_4acC[df_merged_4acC$both == "T", ]

df_all_4acC <- bind_rows(df_rep1_4acC, df_rep2_4acC, df_both_4acC) %>%
  distinct(chr, pos, strand, .keep_all = TRUE)

dmls_rep1 <- callDML(df_rep1_4acC, delta = 0.1, p.threshold = 0.15)
dmls_rep2 <- callDML(df_rep2_4acC, delta = 0.1, p.threshold = 0.15)
dmls <- callDML(df_both_4acC, delta = 0.1, p.threshold = 0.15)

dmls_all <- bind_rows(dmls_rep1, dmls_rep2, dmls) %>%
  distinct(chr, pos, strand, .keep_all = TRUE)

both_sites <- dmls %>%
  distinct(chr, pos, strand) %>%
  mutate(both_flag = TRUE)

df_all_4acC <- df_all_4acC %>%
  left_join(both_sites, by = c("chr", "pos", "strand")) %>%
  mutate(
    type = if_else(
      is.na(both_flag),          
      "filtered",
      "Site"                    
    )
  ) %>%
  dplyr::select(-both_flag)

table(df_all_4acC$type)


write.table(dmls_all, "4acCsite_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(dmls_rep1, "4acCsite_rep1.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(dmls_rep2, "4acCsite_rep2.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(dmls, "4acCsite_filtered.txt", sep = "\t", quote = FALSE, row.names = FALSE)