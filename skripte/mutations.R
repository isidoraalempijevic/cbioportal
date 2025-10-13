# ======================================
# MULTI-CANCER FAM46C MUTATION ANALYSIS — SAMPLE-BASED
# ======================================

# ---- Loading the data ----

melanoma <- read.csv("~/cBioPostal data final/melanoma_fin.csv")

# ---- Library ----
library(tidyr)
library(tidyverse)
library(dplyr)
library(stringr)
library(writexl)
library(ggplot2)

# ---- Addint columns about cancer type ----

melanoma        <- melanoma        |> mutate(Cancer = "Melanoma")

# ---- Making one data frame ----
combined <- melanoma

# ---- Tyding data and adding annotations ----

combined[combined == ""] <- NA

total_samples <- combined %>% distinct(Sample.ID) %>% nrow()
write_xlsx(combined, "combined.xlsx")

# ======================================
# Protein analysis (sample-based)
# ======================================

# 1) Mutation mapped on protein domain
protein_domains_df <- data.frame(
  domain = c("Catalytic", "Central"),
  start = c(14, 229),
  end   = c(228, 343)
)

p <- ggplot() +
  geom_rect(data = protein_domains_df,
            aes(xmin = start, xmax = end, ymin = 0.9, ymax = 1.1, fill = domain), alpha = 0.3) +
  geom_point(data = combined,
             aes(x = Protein.Position, y = 1, color = Cancer, shape = Aggressiveness), size = 3) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(
    x = "Protein Position",
    title = "Mutations Mapped on Protein Domains",
    color = "Cancer Type",
    shape = "Aggressiveness"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold") 
  )

print(p)

ggsave("mutations_maped_protein.png", plot = p, width = 14, height = 12, dpi = 300)


# 2) Density plot

p <- ggplot() +
  geom_rect(data = protein_domains_df,
            aes(xmin = start, xmax = end, ymin = 0, ymax = Inf, fill = domain), alpha = 0.2) +
  geom_density(data = combined,
               aes(x = Protein.Position, color = Aggressiveness), size = 1, na.rm = TRUE) +
  labs(x = "Protein Position", y = "Density",
       title = "Density of Mutations Across Protein with Domain Context") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold")
  )

print(p)

ggsave("density_maped_protein.png", plot = p, width = 14, height = 12, dpi = 300)

# 3) Fisher test 
fisher_table_protein <- combined %>%
  distinct(Sample.ID, Protein_domain, Aggressiveness) %>%
  count(Protein_domain, Aggressiveness) %>%
  pivot_wider(names_from = Aggressiveness, values_from = n, values_fill = 0) %>%
  column_to_rownames("Protein_domain") %>%
  as.matrix()
fisher_protein <- fisher.test(fisher_table_protein)

fisher_df <- data.frame(
  p_value = fisher_protein$p.value,
  method = fisher_protein$method,
  alternative = fisher_protein$alternative
)

write_xlsx(fisher_df, "fisher_protein.xlsx")

# Catalytic vs Central
if (all(c("Catalytic","Central") %in% rownames(fisher_table_protein))) {
  sub_tbl_protein <- fisher_table_protein[c("Catalytic","Central"), , drop = FALSE]
  print(fisher.test(sub_tbl_protein))
}

fisher_result <- fisher.test(sub_tbl_protein)

fisher_df <- data.frame(
  p_value = fisher_result$p.value,
  method = fisher_result$method,
  alternative = fisher_result$alternative
)

write_xlsx(fisher_df, "fisher_sub_protein.xlsx")

# 4) Percentage of Samples with ≥1 FAM46C Mutation by Protein Domain

# Count samples per protein domain and aggressiveness
df_percent_protein <- combined %>%
  distinct(Sample.ID, Protein_domain, Aggressiveness) %>%
  count(Protein_domain, Aggressiveness, name = "n_samples")

# Total samples per aggressiveness group
total_samples_group <- combined %>%
  distinct(Sample.ID, Aggressiveness) %>%
  count(Aggressiveness, name = "total_group")

# Normalize within each group
df_percent_protein <- df_percent_protein %>%
  left_join(total_samples_group, by = "Aggressiveness") %>%
  mutate(Percent = n_samples / total_group * 100)

# Plot
p <- ggplot(df_percent_protein, aes(x = Protein_domain, y = Percent, fill = Aggressiveness)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Percentage of Samples with ≥1 FAM46C Mutation by Protein Domain",
    x = "Protein Domain", y = "% within Aggressiveness Group", fill = "Aggressiveness"
  ) +
  scale_fill_manual(values = c("Aggressive" = "#4D4D4D", "Less_aggressive" = "#C0C0C0")) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold")
  )

# Display the plot
print(p)

ggsave("percentage_maped_protein.png", plot = p, width = 14, height = 16, dpi = 300)

# ======================================
# cDNA analysis (sample-based)
# ======================================

# 1) Mutations mapped on gene regions
gene_domains_df <- data.frame(
  domain = c("5'UTR", "CDS"),
  start = c(1, 134),
  end   = c(133, 1309)
)

p <- ggplot() +
  geom_rect(data = gene_domains_df,
            aes(xmin = start, xmax = end, ymin = 0.9, ymax = 1.1, fill = domain), alpha = 0.3) +
  geom_point(data = combined,
             aes(x = cDNA_pos, y = 1, color = Cancer, shape = Aggressiveness), size = 3) +
  scale_y_continuous(NULL, breaks = NULL) +
  labs(
    x = "Gene Position (cDNA)",
    title = "Mutations Mapped on cDNA Regions",
    color = "Cancer Type",
    shape = "Aggressiveness"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 20, face = "bold")
  )

ggsave("mutations_maped_gene.png", plot = p, width = 14, height = 12, dpi = 300)

# 2) Density plot
p <- ggplot() +
  geom_rect(data = gene_domains_df,
            aes(xmin = start, xmax = end, ymin = 0, ymax = Inf, fill = domain), alpha = 0.2) +
  geom_density(data = combined,
               aes(x = cDNA_pos, color = Aggressiveness), size = 1, na.rm = TRUE) +
  labs(x = "Gene Position (cDNA)", y = "Density",
       title = "Density of Mutations Across Gene with Region Context") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold")
  )

ggsave("density_maped_gene.png", plot = p, width = 14, height = 12, dpi = 300)

# 3) Fisher test
fisher_table_cDNA <- combined %>%
  distinct(Sample.ID, Domain, Aggressiveness) %>%
  count(Domain, Aggressiveness) %>%
  pivot_wider(names_from = Aggressiveness, values_from = n, values_fill = 0) %>%
  column_to_rownames("Domain") %>%
  as.matrix()

fisher_gene <- fisher.test(fisher_table_cDNA)

fisher_df <- data.frame(
  p_value = fisher_gene$p.value,
  method = fisher_gene$method,
  alternative = fisher_gene$alternative
)

write_xlsx(fisher_df, "fisher_gene.xlsx")

# 4) Percentage of Samples with ≥1 FAM46C Mutation by Transcript Region

# Count samples per transcript region and aggressiveness group
df_percent_cDNA <- combined %>%
  distinct(Sample.ID, Domain, Aggressiveness) %>%
  count(Domain, Aggressiveness, name = "n_samples")

# Total samples per aggressiveness group
total_samples_group <- combined %>%
  distinct(Sample.ID, Aggressiveness) %>%
  count(Aggressiveness, name = "total_group")

# Normalize within each group
df_percent_cDNA <- df_percent_cDNA %>%
  left_join(total_samples_group, by = "Aggressiveness") %>%
  mutate(Percent = n_samples / total_group * 100)

# Plot
p <- ggplot(df_percent_cDNA, aes(x = Domain, y = Percent, fill = Aggressiveness)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Percentage of Samples with ≥1 FAM46C Mutation by Transcript Region",
    x = "Transcript Region", y = "% within Aggressiveness Group", fill = "Aggressiveness"
  ) +
  scale_fill_manual(values = c("Aggressive" = "#4D4D4D", "Less_aggressive" = "#C0C0C0")) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold")
  )

# Display the plot
print(p)

ggsave("percentage_maped_gene.png", plot = p, width = 14, height = 16, dpi = 300)

# ======================================
# Analysis of C-terminal (sample-based)
# ======================================

# Flag samples with at least one C-terminal mutation
sample_cterm_flags <- combined %>%
  mutate(C_terminal_flag = Protein.Position > 343) %>%
  group_by(Sample.ID, Aggressiveness) %>%
  summarise(
    has_C_terminal = any(C_terminal_flag, na.rm = TRUE),
    .groups = "drop"
  )

# Count samples by aggressiveness and mutation status
c_terminal_counts <- sample_cterm_flags %>%
  count(Aggressiveness, has_C_terminal, name = "n_samples")

# Total samples per aggressiveness group
total_samples_group <- sample_cterm_flags %>%
  count(Aggressiveness, name = "total_group")

# Normalize within each group
df_cterm_plot <- c_terminal_counts %>%
  left_join(total_samples_group, by = "Aggressiveness") %>%
  mutate(Percent = n_samples / total_group * 100)

# View result
print(df_cterm_plot)

# Filter to show only samples with C-terminal mutations
df_cterm_plot_filtered <- df_cterm_plot %>%
  filter(has_C_terminal)

# Create the plot
p <- ggplot(df_cterm_plot_filtered, aes(x = Aggressiveness, y = Percent, fill = Aggressiveness)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("Aggressive" = "#4D4D4D", "Less_aggressive" = "#C0C0C0")) +
  labs(
    title = "Percentage of Samples with ≥1 C-terminal (aa > 343) Mutation",
    x = "Aggressiveness", y = "% within Aggressiveness Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.position = "none"
  )

# Display the plot
print(p)

ggsave("percentage_cterm.png", plot = p, width = 14, height = 16, dpi = 300)

# Fisher test (Aggressive vs Less_aggressive) × (C-terminal vs No)
cterm_pivot <- sample_cterm_flags %>%
  tidyr::complete(Aggressiveness, has_C_terminal = c(FALSE, TRUE), fill = list()) %>%
  count(Aggressiveness, has_C_terminal, name = "n") %>%
  pivot_wider(names_from = has_C_terminal, values_from = n, values_fill = 0)

a <- cterm_pivot %>% filter(Aggressiveness == "Aggressive") %>% pull(`TRUE`)
b <- cterm_pivot %>% filter(Aggressiveness == "Aggressive") %>% pull(`FALSE`)
c <- cterm_pivot %>% filter(Aggressiveness == "Less_aggressive") %>% pull(`TRUE`)
d <- cterm_pivot %>% filter(Aggressiveness == "Less_aggressive") %>% pull(`FALSE`)

contingency_matrix <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
rownames(contingency_matrix) <- c("Aggressive", "Less_aggressive")
colnames(contingency_matrix) <- c("C_terminal", "Not_C_terminal")
print(contingency_matrix)
print(fisher.test(contingency_matrix))


contingency_df <- as.data.frame.matrix(contingency_matrix)
fisher_result <- fisher.test(contingency_matrix)
fisher_df <- data.frame(
  p_value = fisher_result$p.value,
  method = fisher_result$method,
  alternative = fisher_result$alternative
)

write_xlsx(
  list(
    Contingency_Table = contingency_df,
    Fisher_Test_Result = fisher_df
  ),
  path = "cterm_fisher_results.xlsx"
)

# ======================================
# Motifs of interest
# ======================================

# 1) Binding domain
binding <- combined %>% filter(Protein.Position >= 219, Protein.Position <= 228)
write_xlsx(binding, "roi_binding_219_228.xlsx")

# 2) GS motif
gs_motif <- combined %>% filter(Protein.Position >= 73, Protein.Position <= 74)
write_xlsx(gs_motif, "roi_gs_motif_73_74.xlsx")

# 3) DEhDEh motif
DEhDEh_motif <- combined %>% filter(Protein.Position >= 90, Protein.Position <= 93)
write_xlsx(DEhDEh_motif, "roi_DEhDEh_90_93.xlsx")

# 4) hDEh motif
hDEh_motif <- combined %>% filter(Protein.Position >= 165, Protein.Position <= 167)
write_xlsx(hDEh_motif, "roi_hDEh_165_167.xlsx")

# 5) Plk4 binding domain
plk4_binding <- combined %>%
  filter((Protein.Position >= 140 & Protein.Position <= 146) |
           (Protein.Position >= 318 & Protein.Position <= 323))
write_xlsx(plk4_binding, "roi_plk4_binding.xlsx")

# 6) BCCIP-alpha binding
bccip_sites <- c(181, 188, 191, 192, 194, 158, 184, 193, 201, 203, 295, 141, 215, 190)
bccip_binding <- combined %>% filter(Protein.Position %in% bccip_sites)
write_xlsx(bccip_binding, "roi_bccip_binding.xlsx")

# 7) Poli-A activity
poli_a_sites <- c(90, 92, 166, 73, 74, 248, 282, 268)
poli_a <- combined %>% filter(Protein.Position %in% poli_a_sites)
write_xlsx(poli_a, "roi_poliA_activity.xlsx")

# 8) Conserved sites
pap_sites <- c(77, 290, 298, 72)
pap <- combined %>% filter(Protein.Position %in% pap_sites)
write_xlsx(pap, "roi_conserved_sites.xlsx")


# ======================================
# Samples with > 1 mutations
# ======================================

samples_multiple <- combined %>%
  group_by(Sample.ID) %>%
  mutate(n_mutations_in_sample = n()) %>%
  ungroup() %>%
  filter(n_mutations_in_sample > 1) %>%
  arrange(desc(n_mutations_in_sample))

samples_multiple_summary <- samples_multiple %>%
  distinct(Sample.ID, Cancer, Aggressiveness, n_mutations_in_sample) %>%
  arrange(desc(n_mutations_in_sample))

write_xlsx(samples_multiple, "samples_with_multiple_mutations_full.xlsx")

write_xlsx(samples_multiple_summary, "samples_with_multiple_mutations_summary.xlsx")
