# =========================
# Load Libraries
# =========================
library(cbioportalR)
library(dplyr)
library(purrr)
library(ggplot2)
library(stringr)
library(tidyr)
library(forcats)

# =========================
# 1. Connect to cBioPortal and Retrieve Studies
# =========================
set_cbioportal_db("public")

myeloma_studies_info <- available_studies() %>%
  filter(str_detect(name, regex("myeloma", ignore_case = TRUE))) %>%
  select(studyId, name)

myeloma_studies <- myeloma_studies_info$studyId

myeloma_studies_df <- as.data.frame(myeloma_studies)
#write_xlsx(myeloma_studies_df,
          # path = "C:/Users/isido/OneDrive/Documents/CBIOP/myeloma/myeloma_studies.xlsx")

# =========================
# 2. Get All Mutations and Add Study Names
# =========================
all_mutations <- map_dfr(myeloma_studies, function(study_id) {
  tryCatch(get_mutations_by_study(study_id), error = function(e) NULL)
}) %>%
  left_join(myeloma_studies_info, by = "studyId")

all_clinical <- map_dfr(myeloma_studies, function(study_id) {
  tryCatch(get_clinical_by_study(study_id), error = function(e) NULL)
})

sample_type_clinical <- all_clinical %>%
  filter(clinicalAttributeId == "SAMPLE_TYPE") %>%
  select(sampleId, SAMPLE_TYPE = value) %>%
  distinct()

# =========================
# 3. Classify Samples by Aggressiveness (updated with progression info)
# =========================
# Patients with progression of disease
progression_patients <- all_clinical %>%
  filter(clinicalAttributeId == "PATIENT_DEATH_REASON",
         str_detect(value, regex("Progression of disease", ignore_case = TRUE))) %>%
  pull(patientId) %>%
  unique()

# Map progression patients to sampleIds
progression_samples <- all_mutations %>%
  filter(patientId %in% progression_patients) %>%
  pull(sampleId) %>%
  unique()

# Aggressive by study name or sample type
aggressive_samples_by_study <- unique(c(
  all_mutations %>%
    filter(str_detect(name, regex("metastatic|invasive", ignore_case = TRUE))) %>%
    pull(sampleId),
  sample_type_clinical %>%
    filter(str_detect(SAMPLE_TYPE, regex("metastasis|metastatic", ignore_case = TRUE))) %>%
    pull(sampleId)
))

# Final aggressive samples = study-based + sample type + progression
aggressive_samples <- unique(c(
  aggressive_samples_by_study,
  progression_samples
))

# Classification table
sample_aggressiveness <- all_mutations %>%
  select(sampleId) %>%
  distinct() %>%
  mutate(Aggressiveness = if_else(sampleId %in% aggressive_samples,
                                  "Aggressive", "Less_aggressive"))

# helper: ensure both aggressiveness levels exist for complete()
aggr_levels <- c("Aggressive", "Less_aggressive")

summary_table <- sample_aggressiveness %>%
  count(Aggressiveness) %>%
  complete(Aggressiveness = aggr_levels, fill = list(n = 0)) %>%
  rename(Count = n) %>%
  mutate(Total_Samples = sum(Count)) %>%
  select(Total_Samples, Aggressiveness, Count)

print(summary_table)
write_xlsx(summary_table,
           path = "C:/Users/isido/OneDrive/Documents/cbioportal/summary_table.xlsx")
# =========================
# 4. Mutation Analysis by Sample (TENT5C + ACTB)
# =========================

# Define genes of interest
genes_of_interest <- c("TENT5C", "OAZ1", "SRP14", "TBP", "RPL13A", "TP53", "KRAS", "PIK3CA", "APC")

# Total number of unique samples
total_samples <- all_mutations %>% 
  pull(sampleId) %>% 
  unique() %>% 
  length()

# Count mutations per gene
mut_counts <- all_mutations %>%
  filter(hugoGeneSymbol %in% genes_of_interest) %>%
  distinct(sampleId, hugoGeneSymbol) %>%
  count(hugoGeneSymbol, name = "mut_samples")

# Ensure all genes are represented, even with 0 mutations
mut_counts <- tibble(hugoGeneSymbol = genes_of_interest) %>%
  left_join(mut_counts, by = "hugoGeneSymbol") %>%
  mutate(mut_samples = replace_na(mut_samples, 0),
         no_mut_samples = total_samples - mut_samples)

# Reshape and calculate percentages
df_overall <- mut_counts %>%
  pivot_longer(cols = c("mut_samples", "no_mut_samples"),
               names_to = "Category", values_to = "Samples") %>%
  mutate(Percentage = Samples / total_samples * 100,
         Category = recode(Category,
                           mut_samples = "Mutation",
                           no_mut_samples = "No mutation"),
         hugoGeneSymbol = ifelse(hugoGeneSymbol == "TENT5C", "FAM46C", hugoGeneSymbol))

# Create the plot
p <- ggplot(df_overall, aes(x = Category, y = Percentage, fill = Category)) +
  geom_col() +
  geom_text(aes(label = paste0(round(Percentage, 2), "%")),
            vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("Mutation" = "#4D4D4D", "No mutation" = "#56B4E9")) +
  ylim(0, 100) +
  labs(title = "Mutation percentage in myeloma samples",
       x = "", y = "Percentage (%)") +
  facet_wrap(~ hugoGeneSymbol) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold")
  )

# Display the plot
print(p)

# Save the plot
ggsave("mutation_percentage_myeloma.png", plot = p, width = 12, height = 14, dpi = 300)

## Statistics

# Fisher test per gene (global, mutation vs. no mutation)
# Create a binary mutation matrix: one row per sample, one column per gene

# Get all unique samples
all_samples <- all_mutations %>% distinct(sampleId)

# Create binary mutation matrix: one row per sample, one column per gene
mutation_matrix <- all_mutations %>%
  filter(hugoGeneSymbol %in% genes_of_interest) %>%
  distinct(sampleId, hugoGeneSymbol) %>%
  mutate(mutated = TRUE) %>%
  pivot_wider(names_from = hugoGeneSymbol, values_from = mutated, values_fill = FALSE)

# Add missing samples (those with no mutations in any gene)
mutation_matrix <- all_samples %>%
  left_join(mutation_matrix, by = "sampleId")

# Add missing gene columns (those not present in the data at all)
missing_genes <- setdiff(genes_of_interest, colnames(mutation_matrix))
mutation_matrix[missing_genes] <- FALSE

# Reorder columns: sampleId first, then genes in original order
mutation_matrix <- mutation_matrix %>%
  select(sampleId, all_of(genes_of_interest))

# Fill any remaining NAs with FALSE
mutation_matrix <- mutation_matrix %>%
  mutate(across(all_of(genes_of_interest), ~replace_na(.x, FALSE)))

# Fisher test per gene (only if table has both mutation states)
global_fisher <- map_dfr(genes_of_interest, function(gene) {
  mutated <- mutation_matrix[[gene]]
  global_mutated <- rowSums(mutation_matrix[genes_of_interest]) > 0
  contingency <- table(mutated, global_mutated)
  
  if (all(dim(contingency) == c(2, 2))) {
    test <- fisher.test(contingency)
    tibble(
      hugoGeneSymbol = gene,
      mut_sample_count = sum(mutated),
      fisher_p = test$p.value
    )
  } else {
    tibble(
      hugoGeneSymbol = gene,
      mut_sample_count = sum(mutated),
      fisher_p = NA_real_
    )
  }
}) %>%
  mutate(fisher_p_adj = p.adjust(fisher_p, method = "BH"))

print(global_fisher)

write_xlsx(global_fisher,
           "global_fisher_myeloma.xlsx")


# =========================
# 7. Mutation Percentage in all samples + FAM46 family
# =========================

genes_of_interest <- c("TENT5C", "TENT5A", "TENT5B", "TENT5D", "OAZ1", "SRP14", "TBP", "RPL13A", "TP53", "KRAS", "PIK3CA", "APC")

# Total number of unique samples
total_samples <- all_mutations %>% 
  pull(sampleId) %>% 
  unique() %>% 
  length()

# Count mutations per gene
mut_counts <- all_mutations %>%
  filter(hugoGeneSymbol %in% genes_of_interest) %>%
  distinct(sampleId, hugoGeneSymbol) %>%
  count(hugoGeneSymbol, name = "mut_samples")

# Ensure all genes are represented, even with 0 mutations
mut_counts <- tibble(hugoGeneSymbol = genes_of_interest) %>%
  left_join(mut_counts, by = "hugoGeneSymbol") %>%
  mutate(mut_samples = replace_na(mut_samples, 0),
         no_mut_samples = total_samples - mut_samples)

# Reshape and calculate percentages
df_overall <- mut_counts %>%
  pivot_longer(cols = c("mut_samples", "no_mut_samples"),
               names_to = "Category", values_to = "Samples") %>%
  mutate(Percentage = Samples / total_samples * 100,
         Category = recode(Category,
                           mut_samples = "Mutation",
                           no_mut_samples = "No mutation"))

df_overall <- df_overall %>%
  mutate(hugoGeneSymbol = recode(hugoGeneSymbol,
                                 "TENT5A" = "FAM46A",
                                 "TENT5B" = "FAM46B",
                                 "TENT5C" = "FAM46C",
                                 "TENT5D" = "FAM46D"))

p <- ggplot(df_overall, aes(x = Category, y = Percentage, fill = Category)) +
  geom_col() +
  geom_text(aes(label = paste0(round(Percentage, 2), "%")),
            vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("Mutation" = "#4D4D4D", "No mutation" = "#56B4E9")) +
  ylim(0, 100) +
  labs(title = "Mutation percentage in myeloma samples",
       x = "", y = "Percentage (%)") +
  facet_wrap(~ hugoGeneSymbol) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold")
  )

print(p)

# Save the plot
ggsave("mutation_percentage_fams_myeloma.png", plot = p, width = 14, height = 16, dpi = 300)

## Statistics FAM46A

genes_of_interest <- c("TENT5A", "OAZ1", "SRP14", "TBP", "RPL13A", "TP53", "KRAS", "PIK3CA", "APC")

# Get all unique samples
all_samples <- all_mutations %>% distinct(sampleId)

# Create binary mutation matrix: one row per sample, one column per gene
mutation_matrix <- all_mutations %>%
  filter(hugoGeneSymbol %in% genes_of_interest) %>%
  distinct(sampleId, hugoGeneSymbol) %>%
  mutate(mutated = TRUE) %>%
  pivot_wider(names_from = hugoGeneSymbol, values_from = mutated, values_fill = FALSE)

# Add missing samples (those with no mutations in any gene)
mutation_matrix <- all_samples %>%
  left_join(mutation_matrix, by = "sampleId")

# Add missing gene columns (those not present in the data at all)
missing_genes <- setdiff(genes_of_interest, colnames(mutation_matrix))
mutation_matrix[missing_genes] <- FALSE

# Reorder columns: sampleId first, then genes in original order
mutation_matrix <- mutation_matrix %>%
  select(sampleId, all_of(genes_of_interest))

# Fill any remaining NAs with FALSE
mutation_matrix <- mutation_matrix %>%
  mutate(across(all_of(genes_of_interest), ~replace_na(.x, FALSE)))

# Fisher test per gene (only if table has both mutation states)
global_fisher <- map_dfr(genes_of_interest, function(gene) {
  mutated <- mutation_matrix[[gene]]
  global_mutated <- rowSums(mutation_matrix[genes_of_interest]) > 0
  contingency <- table(mutated, global_mutated)
  
  if (all(dim(contingency) == c(2, 2))) {
    test <- fisher.test(contingency)
    tibble(
      hugoGeneSymbol = gene,
      mut_sample_count = sum(mutated),
      fisher_p = test$p.value
    )
  } else {
    tibble(
      hugoGeneSymbol = gene,
      mut_sample_count = sum(mutated),
      fisher_p = NA_real_
    )
  }
}) %>%
  mutate(fisher_p_adj = p.adjust(fisher_p, method = "BH"))

print(global_fisher)

write_xlsx(global_fisher, "famA_fisher_myeloma.xlsx")

## Statistics FAM46B

genes_of_interest <- c("TENT5B", "OAZ1", "SRP14", "TBP", "RPL13A", "TP53", "KRAS", "PIK3CA", "APC")

# Get all unique samples
all_samples <- all_mutations %>% distinct(sampleId)

# Create binary mutation matrix: one row per sample, one column per gene
mutation_matrix <- all_mutations %>%
  filter(hugoGeneSymbol %in% genes_of_interest) %>%
  distinct(sampleId, hugoGeneSymbol) %>%
  mutate(mutated = TRUE) %>%
  pivot_wider(names_from = hugoGeneSymbol, values_from = mutated, values_fill = FALSE)

# Add missing samples (those with no mutations in any gene)
mutation_matrix <- all_samples %>%
  left_join(mutation_matrix, by = "sampleId")

# Add missing gene columns (those not present in the data at all)
missing_genes <- setdiff(genes_of_interest, colnames(mutation_matrix))
mutation_matrix[missing_genes] <- FALSE

# Reorder columns: sampleId first, then genes in original order
mutation_matrix <- mutation_matrix %>%
  select(sampleId, all_of(genes_of_interest))

# Fill any remaining NAs with FALSE
mutation_matrix <- mutation_matrix %>%
  mutate(across(all_of(genes_of_interest), ~replace_na(.x, FALSE)))

# Fisher test per gene (only if table has both mutation states)
global_fisher <- map_dfr(genes_of_interest, function(gene) {
  mutated <- mutation_matrix[[gene]]
  global_mutated <- rowSums(mutation_matrix[genes_of_interest]) > 0
  contingency <- table(mutated, global_mutated)
  
  if (all(dim(contingency) == c(2, 2))) {
    test <- fisher.test(contingency)
    tibble(
      hugoGeneSymbol = gene,
      mut_sample_count = sum(mutated),
      fisher_p = test$p.value
    )
  } else {
    tibble(
      hugoGeneSymbol = gene,
      mut_sample_count = sum(mutated),
      fisher_p = NA_real_
    )
  }
}) %>%
  mutate(fisher_p_adj = p.adjust(fisher_p, method = "BH"))

print(global_fisher)

write_xlsx(global_fisher, "famB_fisher_myeloma.xlsx")

## Statistics FAM46D

genes_of_interest <- c("TENT5D", "OAZ1", "SRP14", "TBP", "RPL13A", "TP53", "KRAS", "PIK3CA", "APC")

# Get all unique samples
all_samples <- all_mutations %>% distinct(sampleId)

# Create binary mutation matrix: one row per sample, one column per gene
mutation_matrix <- all_mutations %>%
  filter(hugoGeneSymbol %in% genes_of_interest) %>%
  distinct(sampleId, hugoGeneSymbol) %>%
  mutate(mutated = TRUE) %>%
  pivot_wider(names_from = hugoGeneSymbol, values_from = mutated, values_fill = FALSE)

# Add missing samples (those with no mutations in any gene)
mutation_matrix <- all_samples %>%
  left_join(mutation_matrix, by = "sampleId")

# Add missing gene columns (those not present in the data at all)
missing_genes <- setdiff(genes_of_interest, colnames(mutation_matrix))
mutation_matrix[missing_genes] <- FALSE

# Reorder columns: sampleId first, then genes in original order
mutation_matrix <- mutation_matrix %>%
  select(sampleId, all_of(genes_of_interest))

# Fill any remaining NAs with FALSE
mutation_matrix <- mutation_matrix %>%
  mutate(across(all_of(genes_of_interest), ~replace_na(.x, FALSE)))

# Fisher test per gene (only if table has both mutation states)
global_fisher <- map_dfr(genes_of_interest, function(gene) {
  mutated <- mutation_matrix[[gene]]
  global_mutated <- rowSums(mutation_matrix[genes_of_interest]) > 0
  contingency <- table(mutated, global_mutated)
  
  if (all(dim(contingency) == c(2, 2))) {
    test <- fisher.test(contingency)
    tibble(
      hugoGeneSymbol = gene,
      mut_sample_count = sum(mutated),
      fisher_p = test$p.value
    )
  } else {
    tibble(
      hugoGeneSymbol = gene,
      mut_sample_count = sum(mutated),
      fisher_p = NA_real_
    )
  }
}) %>%
  mutate(fisher_p_adj = p.adjust(fisher_p, method = "BH"))

print(global_fisher)

write_xlsx(global_fisher, "famD_fisher_myeloma.xlsx")
