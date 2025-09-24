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
library(writexl)
library(ggpattern)
install.packages("Cairo")
library(Cairo)
install.packages("svglite")
library(svglite)

# =========================
# 1. Connect to cBioPortal and Retrieve Studies
# =========================
set_cbioportal_db("public")

melanoma_studies_info <- available_studies() %>%
  filter(str_detect(name, regex("melanoma", ignore_case = TRUE))) %>%
  select(studyId, name)

melanoma_studies <- melanoma_studies_info$studyId

#saving information about study names

melanoma_studies_df <- as.data.frame(melanoma_studies)
write_xlsx(melanoma_studies_df,
           path = "C:/Users/isido/OneDrive/Documents/cbioportal/tabele/melanoma/melanoma_studies.xlsx")

# =========================
# 2. Get All Mutations and Add Study Names
# =========================
all_mutations <- map_dfr(melanoma_studies, function(study_id) {
  tryCatch(get_mutations_by_study(study_id), error = function(e) NULL)
}) %>%
  left_join(melanoma_studies_info, by = "studyId")

all_clinical <- map_dfr(melanoma_studies, function(study_id) {
  tryCatch(get_clinical_by_study(study_id), error = function(e) NULL)
})

sample_type_clinical <- all_clinical %>%
  filter(clinicalAttributeId == "SAMPLE_TYPE") %>%
  select(sampleId, SAMPLE_TYPE = value) %>%
  distinct()

#checking names of the studies for indexing "Aggressive" and "Less aggressive" samples

unique_names <- all_mutations %>%
  distinct(name) %>%
  pull(name)

print(unique_names)

# =========================
# 3. Classify Samples by Aggressiveness
# =========================
aggressive_samples <- unique(c(
  # From study name
  all_mutations %>% filter(str_detect(name, regex("metastatic|invasive", ignore_case = TRUE))) %>% pull(sampleId),
  # From SAMPLE_TYPE clinical data
  sample_type_clinical %>% filter(str_detect(SAMPLE_TYPE, regex("metastasis|metastatic", ignore_case = TRUE))) %>% pull(sampleId)
))

sample_aggressiveness <- all_mutations %>%
  select(sampleId) %>%
  distinct() %>%
  mutate(Aggressiveness = if_else(sampleId %in% aggressive_samples, "Aggressive", "Less_aggressive"))

aggr_levels <- c("Aggressive", "Less_aggressive")

#making a summery table, total number of samples, number of aggressive samples, number of less aggressive samples

summary_table <- sample_aggressiveness %>%
  count(Aggressiveness) %>%
  complete(Aggressiveness = aggr_levels, fill = list(n = 0)) %>%
  rename(Count = n) %>%
  mutate(Total_Samples = sum(Count)) %>%
  select(Total_Samples, Aggressiveness, Count)

print(summary_table)
write_xlsx(summary_table,
           path = "C:/Users/isido/OneDrive/Documents/cbioportal/tabele/melanoma/summary_table.xlsx")

# =========================
# 4. Mutation Analysis by Sample (TENT5C + positive controls + negative controls)
# =========================
genes_of_interest <- c("TENT5C", "OAZ1", "SRP14", "TBP", "RPL13A", "TP53", "KRAS", "PIK3CA", "APC")

# Total samples
total_samples <- all_mutations %>% pull(sampleId) %>% unique() %>% length()

# Samples with mutations in genes of interest
samples_per_gene <- all_mutations %>%
  filter(hugoGeneSymbol %in% genes_of_interest) %>%
  distinct(sampleId, hugoGeneSymbol)

df_overall <- samples_per_gene %>%
  count(hugoGeneSymbol, name = "mut_samples") %>%
  mutate(no_mut_samples = total_samples - mut_samples) %>%
  pivot_longer(cols = c("mut_samples","no_mut_samples"),
               names_to = "Category", values_to = "Samples") %>%
  mutate(Percentage = Samples / total_samples * 100,
         Category = recode(Category,
                           mut_samples = "Mutation",
                           no_mut_samples = "No mutation"))

# Replace TENT5C with FAM46C in the data
df_overall <- df_overall %>%
  mutate(hugoGeneSymbol = ifelse(hugoGeneSymbol == "TENT5C", "FAM46C", hugoGeneSymbol))

unique(df_overall$Category)

# Create the plot
p <- ggplot(df_overall, aes(x = Category, y = Percentage, fill = Category)) +
  geom_col() +
  geom_text(aes(label = paste0(round(Percentage, 2), "%")),
            vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("Mutation" = "#4D4D4D", "No mutation" = "#56B4E9")) +
  ylim(0, 100) +
  labs(title = "Mutation percentage in melanoma samples",
       x = "", y = "Percentage (%)") +
  facet_wrap(~ hugoGeneSymbol) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold")
  )

print(p)

# Save the plot
ggsave("grafici/melanoma/mutation_percentage_melanoma.png", plot = p, width = 12, height = 14, dpi = 300)

## Statistika

# Fisher test per gene (global, mutation vs. no mutation)
# Create a binary mutation matrix: one row per sample, one column per gene
mutation_matrix <- all_mutations %>%
  filter(hugoGeneSymbol %in% genes_of_interest) %>%
  distinct(sampleId, hugoGeneSymbol) %>%
  mutate(mutated = TRUE) %>%
  pivot_wider(names_from = hugoGeneSymbol, values_from = mutated, values_fill = FALSE)

# Add missing samples (no mutations in any gene)
all_samples <- all_mutations %>% distinct(sampleId)
mutation_matrix <- all_samples %>%
  left_join(mutation_matrix, by = "sampleId") %>%
  mutate(across(all_of(genes_of_interest), ~replace_na(.x, FALSE)))

# Fisher test per gene
global_fisher <- map_dfr(genes_of_interest, function(gene) {
  mutated <- mutation_matrix[[gene]]
  table <- table(mutated, rowSums(mutation_matrix[genes_of_interest]) > 0)
  test <- fisher.test(table)
  
  tibble(
    hugoGeneSymbol = gene,
    mut_sample_count = sum(mutated),  # Number of samples with mutation in this gene
    fisher_p = test$p.value
  )
}) %>%
  mutate(fisher_p_adj = p.adjust(fisher_p, method = "BH"))

print(global_fisher)

write_xlsx(global_fisher,
           "global_fisher_melanoma.xlsx")

# =========================
# 5. Mutation Percentage by Aggressiveness (sample-level)
# =========================
total_samples_group <- all_mutations %>%
  select(sampleId) %>%
  distinct() %>%
  left_join(sample_aggressiveness, by = "sampleId") %>%
  count(Aggressiveness, name = "total_samples")

mut_samples_group <- all_mutations %>%
  filter(hugoGeneSymbol %in% genes_of_interest) %>%
  select(sampleId, hugoGeneSymbol) %>%
  distinct() %>%
  left_join(sample_aggressiveness, by = "sampleId") %>%
  count(Aggressiveness, hugoGeneSymbol, name = "mut_samples")

df_group_sample <- left_join(total_samples_group, mut_samples_group, 
                             by = "Aggressiveness") %>%
  mutate(mut_samples = ifelse(is.na(mut_samples), 0, mut_samples),
         perc_mut = mut_samples / total_samples * 100,
         perc_non_mut = 100 - perc_mut)

df_plot_group_sample <- df_group_sample %>%
  pivot_longer(cols = c("perc_mut","perc_non_mut"),
               names_to = "Category", values_to = "Percentage") %>%
  mutate(Category = recode(Category,
                           perc_mut = "Mutation",
                           perc_non_mut = "No mutation"))

# Rename gene from TENT5C to FAM46C
df_plot_group_sample <- df_plot_group_sample %>%
  mutate(hugoGeneSymbol = ifelse(hugoGeneSymbol == "TENT5C", "FAM46C", hugoGeneSymbol))

# Create the plot with hatch pattern for aggressive samples
p <- ggplot(df_plot_group_sample, aes(x = Aggressiveness, y = Percentage, fill = Category)) +
  geom_col_pattern(aes(pattern = Aggressiveness),
                   position = position_dodge(width = 0.9),
                   pattern_density = 0.1,
                   pattern_spacing = 0.05,
                   pattern_fill = "black",
                   pattern_colour = "black") +
  geom_text(aes(label = paste0(round(Percentage, 2), "%")),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Mutation" = "#4D4D4D", "No mutation" = "#56B4E9"),
                    name = "Mutation Status") +
  scale_pattern_manual(values = c("Aggressive" = "stripe", "Non-aggressive" = "none"),
                       name = "Aggressiveness",
                       guide = "legend") +  # keep this legend
  guides(fill = guide_legend(override.aes = list(pattern = "none"))) +  # remove hatch from fill legend
  ylim(0, 100) +
  labs(title = "Mutations by sample aggressiveness in melanoma samples",
       x = "Aggressiveness", y = "Percentage (%)") +
  facet_wrap(~ hugoGeneSymbol) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(size = 20, face = "bold"))

# Display the plot
print(p)

ggpattern::ggsave_pattern("grafici/melanoma/mutation_percentage_melanoma_agg.png",
                          plot = p, width = 12, height = 14, dpi = 300)


###########statistika##########

group_table <- df_group_sample %>%
  select(Aggressiveness, hugoGeneSymbol, mut_samples, total_samples) %>%
  mutate(no_mut_samples = total_samples - mut_samples)

group_fisher <- map_dfr(genes_of_interest, function(gene) {
  gene_data <- group_table %>% filter(hugoGeneSymbol == gene)
  
  # Ensure both levels are present
  aggr_data <- gene_data %>% complete(Aggressiveness = aggr_levels, fill = list(mut_samples = 0, no_mut_samples = 0, total_samples = 0))
  
  # Extract counts
  a <- aggr_data %>% filter(Aggressiveness == "Aggressive") %>% pull(mut_samples)
  b <- aggr_data %>% filter(Aggressiveness == "Aggressive") %>% pull(no_mut_samples)
  c <- aggr_data %>% filter(Aggressiveness == "Less_aggressive") %>% pull(mut_samples)
  d <- aggr_data %>% filter(Aggressiveness == "Less_aggressive") %>% pull(no_mut_samples)
  
  test <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
  
  tibble(
    hugoGeneSymbol = gene,
    mut_samples_aggressive = a,
    mut_samples_less_aggressive = c,
    fisher_p = test$p.value
  )
}) %>%
  mutate(fisher_p_adj = p.adjust(fisher_p, method = "BH"))

group_fisher

write_xlsx(group_fisher,
           path = "C:/Users/isido/OneDrive/Documents/CBIOP/melanoma/group_fisher_melanoma.xlsx")


######aggressive samples

aggressive_only <- sample_aggressiveness %>%
  filter(Aggressiveness == "Aggressive") %>%
  pull(sampleId)

# Total aggressive samples
total_aggressive <- all_mutations %>%
  filter(sampleId %in% aggressive_only) %>%
  distinct(sampleId) %>%
  summarise(total_samples = n())

# Mutated samples per gene
mut_samples_aggressive <- all_mutations %>%
  filter(sampleId %in% aggressive_only, hugoGeneSymbol %in% genes_of_interest) %>%
  distinct(sampleId, hugoGeneSymbol) %>%
  count(hugoGeneSymbol, name = "mut_samples") %>%
  complete(hugoGeneSymbol = genes_of_interest, fill = list(mut_samples = 0)) %>%
  mutate(no_mut_samples = total_aggressive$total_samples - mut_samples)

# Prepare for plotting
df_aggressive_plot <- mut_samples_aggressive %>%
  pivot_longer(cols = c("mut_samples", "no_mut_samples"),
               names_to = "Category", values_to = "Samples") %>%
  mutate(Percentage = Samples / total_aggressive$total_samples * 100,
         Category = recode(Category,
                           mut_samples = "Mutation",
                           no_mut_samples = "No mutation"))

# Plot
ggplot(df_aggressive_plot, aes(x = Category, y = Percentage, fill = Category)) +
  geom_col() +
  geom_text(aes(label = paste0(round(Percentage, 2), "%")),
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  ylim(0, 100) +
  labs(title = "Mutation percentage in aggressive melanoma cancer samples",
       x = "", y = "Percentage (%)") +
  facet_wrap(~ hugoGeneSymbol) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


# Define your gene list (replace with your actual list)
# genes_of_interest <- c("GAPDH", "TPB", "RPL13A", "TP53", "KRAS", "APC")

# Get aggressive sample IDs
aggressive_only <- sample_aggressiveness %>%
  filter(Aggressiveness == "Aggressive") %>%
  pull(sampleId)

# Total aggressive sample count
total_aggressive <- length(aggressive_only)

# Create binary mutation matrix for aggressive samples only
mutation_matrix <- all_mutations %>%
  filter(sampleId %in% aggressive_only, hugoGeneSymbol %in% genes_of_interest) %>%
  distinct(sampleId, hugoGeneSymbol) %>%
  mutate(mutated = TRUE) %>%
  pivot_wider(names_from = hugoGeneSymbol, values_from = mutated, values_fill = FALSE)

# Add missing aggressive samples (no mutations in any gene)
all_aggressive_samples <- tibble(sampleId = aggressive_only)
mutation_matrix <- all_aggressive_samples %>%
  left_join(mutation_matrix, by = "sampleId")

# Add missing gene columns with FALSE
missing_genes <- setdiff(genes_of_interest, colnames(mutation_matrix))
for (gene in missing_genes) {
  mutation_matrix[[gene]] <- FALSE
}

# Fill NAs with FALSE
mutation_matrix <- mutation_matrix %>%
  mutate(across(all_of(genes_of_interest), ~replace_na(.x, FALSE)))

# Run Fisher test per gene: mutated vs. non-mutated within aggressive samples
fisher_aggressive <- map_dfr(genes_of_interest, function(gene) {
  mutated <- mutation_matrix[[gene]]
  other_mutated <- rowSums(mutation_matrix[genes_of_interest]) > 0
  
  contingency_table <- table(mutated, other_mutated)
  
  # Skip invalid tables
  if (nrow(contingency_table) < 2 || ncol(contingency_table) < 2 || length(unique(c(contingency_table))) == 1) {
    return(tibble(
      hugoGeneSymbol = gene,
      mut_sample_count = sum(mutated),
      fisher_p = NA_real_
    ))
  }
  
  test <- fisher.test(contingency_table)
  
  tibble(
    hugoGeneSymbol = gene,
    mut_sample_count = sum(mutated),
    fisher_p = test$p.value
  )
}) %>%
  mutate(fisher_p_adj = p.adjust(fisher_p, method = "BH"))

# View results
print(fisher_aggressive)

write_xlsx(fisher_aggressive,
           path = "C:/Users/isido/OneDrive/Documents/CBIOP/melanoma/fisher_aggressive_melanoma.xlsx")



##########################

genes_of_interest <- c("TENT5C", "TENT5A", "TENT5B", "TENT5D", "ACTB", "GAPDH", "TBP", "RPL13A")

# Total samples
total_samples <- all_mutations %>% pull(sampleId) %>% unique() %>% length()

# Samples with mutations in genes of interest
samples_per_gene <- all_mutations %>%
  filter(hugoGeneSymbol %in% genes_of_interest) %>%
  distinct(sampleId, hugoGeneSymbol)

df_overall <- samples_per_gene %>%
  count(hugoGeneSymbol, name = "mut_samples") %>%
  mutate(no_mut_samples = total_samples - mut_samples) %>%
  pivot_longer(cols = c("mut_samples","no_mut_samples"),
               names_to = "Category", values_to = "Samples") %>%
  mutate(Percentage = Samples / total_samples * 100,
         Category = recode(Category,
                           mut_samples = "Mutation",
                           no_mut_samples = "No mutation"))

ggplot(df_overall, aes(x = Category, y = Percentage, fill = Category)) +
  geom_col() +
  geom_text(aes(label = paste0(round(Percentage, 2), "%")),
            vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  ylim(0, 100) +
  labs(title = "Mutation percentage in melanoma samples",
       x = "", y = "Percentage (%)") +
  facet_wrap(~ hugoGeneSymbol) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


#############statistika#########

mutation_matrix <- all_mutations %>%
  filter(hugoGeneSymbol %in% genes_of_interest) %>%
  distinct(sampleId, hugoGeneSymbol) %>%
  mutate(mutated = TRUE) %>%
  pivot_wider(names_from = hugoGeneSymbol, values_from = mutated, values_fill = FALSE)

# Add missing samples
all_samples <- all_mutations %>% distinct(sampleId)
mutation_matrix <- all_samples %>%
  left_join(mutation_matrix, by = "sampleId") %>%
  mutate(across(all_of(genes_of_interest), ~replace_na(.x, FALSE)))

global_fisher <- map_dfr(genes_of_interest, function(gene) {
  mutated <- mutation_matrix[[gene]]
  table <- table(mutated, rowSums(mutation_matrix[genes_of_interest]) > 0)
  test <- fisher.test(table)
  
  tibble(
    hugoGeneSymbol = gene,
    mut_sample_count = sum(mutated),
    fisher_p = test$p.value
  )
}) %>%
  mutate(fisher_p_adj = p.adjust(fisher_p, method = "BH"))

global_fisher

write_xlsx(global_fisher,
           path = "C:/Users/isido/OneDrive/Documents/CBIOP/melanoma/fams_fisher_melanoma.xlsx")



# =========================
# 5. Mutation Percentage by Aggressiveness (sample-level)
# =========================
total_samples_group <- all_mutations %>%
  select(sampleId) %>%
  distinct() %>%
  left_join(sample_aggressiveness, by = "sampleId") %>%
  count(Aggressiveness, name = "total_samples")

mut_samples_group <- all_mutations %>%
  filter(hugoGeneSymbol %in% genes_of_interest) %>%
  select(sampleId, hugoGeneSymbol) %>%
  distinct() %>%
  left_join(sample_aggressiveness, by = "sampleId") %>%
  count(Aggressiveness, hugoGeneSymbol, name = "mut_samples")

df_group_sample <- left_join(total_samples_group, mut_samples_group, 
                             by = "Aggressiveness") %>%
  mutate(mut_samples = ifelse(is.na(mut_samples), 0, mut_samples),
         perc_mut = mut_samples / total_samples * 100,
         perc_non_mut = 100 - perc_mut)

df_plot_group_sample <- df_group_sample %>%
  pivot_longer(cols = c("perc_mut","perc_non_mut"),
               names_to = "Category", values_to = "Percentage") %>%
  mutate(Category = recode(Category,
                           perc_mut = "Mutation",
                           perc_non_mut = "No mutation"))

ggplot(df_plot_group_sample, aes(x = Aggressiveness, y = Percentage, fill = Category)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = paste0(round(Percentage, 2), "%")),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  ylim(0, 100) +
  labs(title = "Mutations by sample aggressiveness",
       x = "Aggressiveness", y = "Percentage (%)") +
  facet_wrap(~ hugoGeneSymbol) +
  theme_minimal(base_size = 14)

########statistika########

group_table <- df_group_sample %>%
  select(Aggressiveness, hugoGeneSymbol, mut_samples, total_samples) %>%
  mutate(no_mut_samples = total_samples - mut_samples)

group_fisher <- map_dfr(genes_of_interest, function(gene) {
  gene_data <- group_table %>% filter(hugoGeneSymbol == gene)
  
  # Ensure both levels exist
  gene_data <- gene_data %>%
    complete(Aggressiveness = aggr_levels, fill = list(mut_samples = 0, no_mut_samples = 0, total_samples = 0))
  
  # Extract counts
  a <- gene_data %>% filter(Aggressiveness == "Aggressive") %>% pull(mut_samples)
  b <- gene_data %>% filter(Aggressiveness == "Aggressive") %>% pull(no_mut_samples)
  c <- gene_data %>% filter(Aggressiveness == "Less_aggressive") %>% pull(mut_samples)
  d <- gene_data %>% filter(Aggressiveness == "Less_aggressive") %>% pull(no_mut_samples)
  
  test <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))
  
  tibble(
    hugoGeneSymbol = gene,
    mut_samples_aggressive = a,
    mut_samples_less_aggressive = c,
    fisher_p = test$p.value
  )
}) %>%
  mutate(fisher_p_adj = p.adjust(fisher_p, method = "BH"))

write_xlsx(group_fisher,
           path = "C:/Users/isido/OneDrive/Documents/CBIOP/melanoma/fams_group_fisher_melanoma.xlsx")

######aggressive

aggressive_only <- sample_aggressiveness %>%
  filter(Aggressiveness == "Aggressive") %>%
  pull(sampleId)

# Total aggressive samples
total_aggressive <- all_mutations %>%
  filter(sampleId %in% aggressive_only) %>%
  distinct(sampleId) %>%
  summarise(total_samples = n())

# Mutated samples per gene
mut_samples_aggressive <- all_mutations %>%
  filter(sampleId %in% aggressive_only, hugoGeneSymbol %in% genes_of_interest) %>%
  distinct(sampleId, hugoGeneSymbol) %>%
  count(hugoGeneSymbol, name = "mut_samples") %>%
  complete(hugoGeneSymbol = genes_of_interest, fill = list(mut_samples = 0)) %>%
  mutate(no_mut_samples = total_aggressive$total_samples - mut_samples)

# Prepare for plotting
df_aggressive_plot <- mut_samples_aggressive %>%
  pivot_longer(cols = c("mut_samples", "no_mut_samples"),
               names_to = "Category", values_to = "Samples") %>%
  mutate(Percentage = Samples / total_aggressive$total_samples * 100,
         Category = recode(Category,
                           mut_samples = "Mutation",
                           no_mut_samples = "No mutation"))

# Plot
ggplot(df_aggressive_plot, aes(x = Category, y = Percentage, fill = Category)) +
  geom_col() +
  geom_text(aes(label = paste0(round(Percentage, 2), "%")),
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  ylim(0, 100) +
  labs(title = "Mutation percentage in aggressive melanoma cancer samples",
       x = "", y = "Percentage (%)") +
  facet_wrap(~ hugoGeneSymbol) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")


# Define your gene list (replace with your actual list)
# genes_of_interest <- c("GAPDH", "TPB", "RPL13A", "TP53", "KRAS", "APC")

# Get aggressive sample IDs
aggressive_only <- sample_aggressiveness %>%
  filter(Aggressiveness == "Aggressive") %>%
  pull(sampleId)

# Total aggressive sample count
total_aggressive <- length(aggressive_only)

# Create binary mutation matrix for aggressive samples only
mutation_matrix <- all_mutations %>%
  filter(sampleId %in% aggressive_only, hugoGeneSymbol %in% genes_of_interest) %>%
  distinct(sampleId, hugoGeneSymbol) %>%
  mutate(mutated = TRUE) %>%
  pivot_wider(names_from = hugoGeneSymbol, values_from = mutated, values_fill = FALSE)

# Add missing aggressive samples (no mutations in any gene)
all_aggressive_samples <- tibble(sampleId = aggressive_only)
mutation_matrix <- all_aggressive_samples %>%
  left_join(mutation_matrix, by = "sampleId")

# Add missing gene columns with FALSE
missing_genes <- setdiff(genes_of_interest, colnames(mutation_matrix))
for (gene in missing_genes) {
  mutation_matrix[[gene]] <- FALSE
}

# Fill NAs with FALSE
mutation_matrix <- mutation_matrix %>%
  mutate(across(all_of(genes_of_interest), ~replace_na(.x, FALSE)))

# Run Fisher test per gene: mutated vs. non-mutated within aggressive samples
fisher_aggressive <- map_dfr(genes_of_interest, function(gene) {
  mutated <- mutation_matrix[[gene]]
  other_mutated <- rowSums(mutation_matrix[genes_of_interest]) > 0
  
  contingency_table <- table(mutated, other_mutated)
  
  # Skip invalid tables
  if (nrow(contingency_table) < 2 || ncol(contingency_table) < 2 || length(unique(c(contingency_table))) == 1) {
    return(tibble(
      hugoGeneSymbol = gene,
      mut_sample_count = sum(mutated),
      fisher_p = NA_real_
    ))
  }
  
  test <- fisher.test(contingency_table)
  
  tibble(
    hugoGeneSymbol = gene,
    mut_sample_count = sum(mutated),
    fisher_p = test$p.value
  )
}) %>%
  mutate(fisher_p_adj = p.adjust(fisher_p, method = "BH"))

# View results
print(fisher_aggressive)

write_xlsx(fisher_aggressive,
           path = "C:/Users/isido/OneDrive/Documents/CBIOP/melanoma/fisher_aggressive_fams_melanoma.xlsx")
