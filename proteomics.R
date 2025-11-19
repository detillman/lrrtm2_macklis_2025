## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(biomaRt)
library(readxl)
library(VennDiagram)
library(DEP)
library(SummarizedExperiment)
library(ggrepel)
library(Mus.musculus)
library(clusterProfiler)
library(tidyverse)
go_reduce <- function (pathway_df, orgdb = "org.Hs.eg.db", threshold = 0.7,
                       scores = NULL, measure = "Wang")
{
  if (!measure %in% c("Resnik", "Lin", "Rel", "Jiang", "Wang")) {
    stop("Chosen measure is not one of the recognised measures, c(\"Resnik\", \"Lin\", \"Rel\", \"Jiang\", \"Wang\").")
  }
  if (measure == "Wang") {
    computeIC <- FALSE
  }
  else {
    computeIC <- TRUE
  }
  ont <- pathway_df %>% .[["go_type"]] %>% unique()
  if (any(!ont %in% c("BP", "CC", "MF"))) {
    stop("Column go_type does not contain the recognised sub-ontologies, c(\"BP\", \"CC\", \"MF\")")
  }
  go_similarity <- setNames(object = vector(mode = "list",
                                            length = length(ont)), nm = ont)
  for (i in 1:length(ont)) {
    print(stringr::str_c("Reducing sub-ontology: ", ont[i]))
    hsGO <- GOSemSim::godata(OrgDb = orgdb, ont = ont[i],
                             computeIC = computeIC)
    terms <- pathway_df %>% dplyr::filter(go_type ==
                                            ont[i]) %>% .[["go_id"]] %>% unique()
    sim <- GOSemSim::mgoSim(GO1 = terms, GO2 = terms, semData = hsGO,
                            measure = measure, combine = NULL)
    go_similarity[[i]] <- rrvgo::reduceSimMatrix(simMatrix = sim,
                                                 threshold = threshold, 
                                                 orgdb = orgdb, scores = scores) %>%
      tibble::as_tibble() %>% dplyr::rename(parent_id = parent,
                                            parent_term = parentTerm, 
                                            parent_sim_score = termDispensability)
  }
  go_sim_df <- go_similarity %>% qdapTools::list_df2df(col1 = "go_type")
  pathway_go_sim_df <- pathway_df %>% dplyr::inner_join(go_sim_df %>%
                                                          dplyr::select(go_type, 
                                                                        go_id = go, 
                                                                        contains("parent")),
                                                        by = c("go_type", "go_id")) %>% dplyr::arrange(go_type,
                                                                                                       parent_id, -parent_sim_score)
  return(pathway_go_sim_df)
}
cbPalette <- c("black", # white 
               "#E69F00", # orange
               "#56B4E9", # light blue
               "#009E73", # green
               "#CC79A7", # pink
               "#D55E00", # red
               "#0072B2", # dark blue
               "#F0E442") # yellow


## ----databaseLoad------------------------------------------------------------------------------------------------------------------------
ensembl <- useEnsembl(biomart = "genes")
ensembl_dataset <- useDataset(dataset = "mmusculus_gene_ensembl", 
                              mart = ensembl)
biomart_1 <- getBM(attributes = c("ensembl_gene_id", 
                                  "uniprot_gn_id", "uniprot_isoform"),
               mart = ensembl_dataset)
biomart_2 <- getBM(attributes = c("ensembl_gene_id", 
                                  "uniprotswissprot", "uniprotsptrembl"),
               mart = ensembl_dataset)
ensembl_uniprot <- biomart_1 %>%
  pivot_longer(cols = contains("uniprot"), 
               names_to = "database", values_to = "uniprot_id") %>%
  bind_rows(biomart_2 %>%
              pivot_longer(cols = contains("uniprot"), 
                           names_to = "database", values_to = "uniprot_id")) %>%
  distinct(ensembl_gene_id, uniprot_id) %>%
  filter(uniprot_id != "")
biomart_3 <- getBM(attributes = c("ensembl_gene_id", 
                                  "entrezgene_id", "external_gene_name"),
               mart = ensembl_dataset) 
final_mart <- full_join(ensembl_uniprot, biomart_3, by = "ensembl_gene_id")


## ----proteinDetection--------------------------------------------------------------------------------------------------------------------
# loading dt somata data (from "epic Dustin Tillman SAM7669 TMT12pro Uniprot Mouse report.xlsx")
dt_soma <- read_excel("input_tables.xlsx", sheet = "Input Table 1") %>%
  filter(grepl("contaminant", Description) == F) %>%
  filter(Contaminant == F) %>%
  select(-contains("Found")) %>%
  rename(uniprot_id = Accession) %>%
  filter(!(if_all(contains("Sample"), ~ is.na(.)))) %>%
  left_join(ensembl_uniprot, by = "uniprot_id") %>%
  mutate(ensembl_gene_id = case_when(is.na(ensembl_gene_id) == T ~ uniprot_id,
                                     T ~ ensembl_gene_id))

# loading pou somata data (from Poulopoulos*, Murphy* et al., Nature, 2019 - Supp Table 2)
pou_soma <- read_excel("external_tables.xlsx", sheet = "pou_somata") %>%
  mutate(external_gene_name = gsub(";.*", "", `Gene names`)) %>%
  drop_na(external_gene_name) %>%
  left_join(biomart_3 %>% distinct(ensembl_gene_id, external_gene_name), 
            by = "external_gene_name") %>%
  mutate(ensembl_gene_id = case_when(is.na(ensembl_gene_id) == T ~ `external_gene_name`,
                                     T ~ ensembl_gene_id))

# venn diagram of dt and pou somata data
venn.diagram(x = list(dt_soma$ensembl_gene_id, 
                      pou_soma$ensembl_gene_id),
             category.names = c("dt", "pou"),
             filename = "supp_fig_1e.tiff",
             output = T)
supp_fig_1e_table <- dt_soma %>%
  distinct(ensembl_gene_id) %>%
  mutate(data = "dt_soma") %>%
  bind_rows(pou_soma %>%
              distinct(ensembl_gene_id) %>%
              mutate(data = "pou_soma"))

# loading pilot gc data (from "pilot Dustin Tillman  SAM08197-8200 LFQ report.xlsx"
pilot_gcs <- read_excel("input_tables.xlsx", sheet = "Input Table 2") %>%
  filter(grepl("contaminant", Description) == F) %>%
  select(Accession, contains("Prep")) %>%
  rename(uniprot_id = Accession,
         "FASP_1" = "Prep A",
         "FASP_2" = "Prep B",
         "ISD_1" = "Prep C",
         "ISD_2" = "Prep D") %>%
  drop_na() %>%
  filter(!(FASP_1 == "NaN" & FASP_2 == "NaN" & ISD_1 == "NaN" & ISD_2 == "NaN")) %>%
  left_join(ensembl_uniprot, by = "uniprot_id") %>%
  mutate(ensembl_gene_id = case_when(is.na(ensembl_gene_id) == T ~ uniprot_id,
                                     T ~ ensembl_gene_id))

# making supplemental table 1
supp_table_1 <- pilot_gcs %>%
  distinct()
write.csv(supp_table_1, "supp_table_1.csv", quote = F, row.names = F)

# loading pou gc data (from Poulopoulos*, Murphy* et al., Nature, 2019 - Supp Table 1)
pou_gcs <- read_excel("external_tables.xlsx", sheet = "pou_gcs") %>%
  mutate(ensembl_gene_id = gsub(";.*", "", ENSG)) %>%
  mutate(ensembl_gene_id = case_when(is.na(ensembl_gene_id) == T ~ `Gene names`,
                                     T ~ ensembl_gene_id))

# venn diagram of pilot and pou gc data
venn.diagram(x = list(pilot_gcs$ensembl_gene_id, 
                      pou_gcs$ensembl_gene_id),
             category.names = c("pilot", "pou"),
             filename = "supp_fig_1f.tiff",
             output = T)
supp_fig_1f_table <- pilot_gcs %>%
  distinct(ensembl_gene_id) %>%
  mutate(data = "pilot_gc") %>%
  bind_rows(pou_gcs %>%
              distinct(ensembl_gene_id) %>%
              mutate(data = "pou_gc"))

# loading wGC data
wGC_data <- read_excel("input_tables.xlsx", sheet = "Input Table 3") %>%
  filter(is.na(Filter) == T) %>%
  filter(!if_any(everything(), ~ grepl("FGCZ", .))) %>%
  separate(Accession, into = c("source", "uniprot_id", "temp"), 
           sep = "\\|") %>%
  filter(!(if_all(contains("TotInt"), ~ . == 0))) %>%
  left_join(ensembl_uniprot, by = "uniprot_id") %>%
  mutate(ensembl_gene_id = case_when(is.na(ensembl_gene_id) == T ~ uniprot_id,
                                     T ~ ensembl_gene_id))


## ----somataProcess-----------------------------------------------------------------------------------------------------------------------
# loading somata data
data_somata <- dt_soma %>%
  rename(Sample.FLyesCre_1 = `Sample 1`,
         Ignore.WTnoCre_1 = `Sample 2`,
         Ignore.FLnoCre_1 = `Sample 3`,
         Sample.FLyesCre_2 = `Sample 4`,
         Sample.FLyesCre_3 = `Sample 5`,
         Ignore.WTnoCre_2 = `Sample 6`,
         Sample.FLyesCre_4 = `Sample 7`,
         Sample.WTyesCre_1 = `Sample 8`,
         Sample.WTyesCre_2 = `Sample 9`,
         Sample.WTyesCre_3 = `Sample 10`,
         Ignore.FLnoCre_2 = `Sample 11`,
         Ignore.FLnoCre_3 = `Sample 12`) %>%
  left_join(final_mart %>%
              distinct(external_gene_name, ensembl_gene_id), by = "ensembl_gene_id") %>%
  mutate(uniprot_id = gsub("-.*", "", uniprot_id)) %>%
  left_join(read_excel("external_tables.xlsx", sheet = "somata_unmapped") %>%
              rename(uniprot_id = From,
                     external_gene_name = To),
            by = "uniprot_id") %>%
  mutate(external_gene_name = case_when(is.na(external_gene_name.x) == T ~
                                          external_gene_name.y,
                                        is.na(external_gene_name.y) == T ~
                                        external_gene_name.x,
                                        T ~ external_gene_name.x),
        .keep = "unused")
  
# make unique names using external gene names and uniprot ids
data_unique_somata <- make_unique(data_somata, "external_gene_name", "uniprot_id", 
                               delim = ";")

# parse condition information from the column names
columns_somata <- grep("Sample", colnames(data_unique_somata))
data_se_somata <- make_se_parse(data_unique_somata, columns_somata)

# barplot of the protein identification overlap between samples
plot_frequency(data_se_somata)

# loose filtering (identified in all - n replicates of at least one condition)
data_filt_somata <- filter_missval(data_se_somata, thr = 1)

# barplot of the number of identified proteins per sample
plot_numbers(data_filt_somata)

# barplot of the protein identification overlap between samples
plot_coverage(data_filt_somata)

# normalize the data
data_norm_somata <- normalize_vsn(data_filt_somata)

# visualize normalization for all samples before and after normalization
plot_normalization(data_filt_somata, data_norm_somata)

# heatmap of proteins with missing values
plot_missval(data_filt_somata)

# intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt_somata)

# impute via quantile regression-based left-censored function (QRLIC)
set.seed(440) # setting random seed for consistent imputation draws
data_imp_somata <- impute(data_norm_somata, fun = "QRILC")

# intensity distributions before and after imputation
plot_imputation(data_filt_somata, data_norm_somata, data_imp_somata)

# identify differential genes
data_diff_manual_somata <- test_diff(data_imp_somata, type = "manual",
                              test = c("FLyesCre__vs_WTyesCre_"))

# denote significant proteins (modify "alpha" and "lfc" as needed)
dep_somata <- add_rejections(data_diff_manual_somata, alpha = 0.1, lfc = log2(0))

# pca for top 500 proteins
pca_var_somata <- apply(assay(dep_somata), 1, sd)
pca_var_df_somata <- assay(dep_somata)[order(pca_var_somata, 
                       decreasing = TRUE)[seq_len(min(c(nrow(dep_somata), 500)))],]
pca_somata <- prcomp(t(pca_var_df_somata), scale = FALSE)
pca_df_somata <- pca_somata$x %>% data.frame() %>% rownames_to_column() %>% 
  left_join(., data.frame(colData(dep_somata)), by = c(rowname = "ID"))
pca_percent_somata <- round(100 * pca_somata$sdev^2/sum(pca_somata$sdev^2), 1)

# plotting pca for top 500 proteins
pca_somata <- ggplot(pca_df_somata %>%
         mutate(Age = "P3",
                Bcl11a = case_when(condition == "FLyesCre_" ~ "Null",
                                   T ~ "WT"),
                Bcl11a = as.factor(Bcl11a),
                Bcl11a = fct_relevel(Bcl11a, "WT", "Null")), 
       aes(get(paste0("PC", 1)), 
           get(paste0("PC", 2)),
           shape = Age,
           color = Bcl11a,
           fill = Bcl11a)) + 
  labs(title = "PC Analysis", 
       x = paste0("PC", 1, ": ", pca_percent_somata[1], "%"), 
       y = paste0("PC", 2, ": ", pca_percent_somata[2], "%")) + 
  coord_fixed() +
  theme_classic(base_size = 7) +
  scale_shape_manual(values = c(21)) +
  scale_color_manual(values = c("#009E73", "#CC79A7")) +
  scale_fill_manual(values = c("#009E73", "#CC79A7")) +
  geom_point() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.3)); pca_somata
ggsave("supp_fig_1j.pdf", pca_somata, units = "in", width = 2.52, height = 1.77)

# volcano plot for specific contrast
plot_volcano(dep_somata, contrast = "FLyesCre__vs_WTyesCre_",
             label_size = 2, add_names = TRUE)

# results table
data_results_somata <- get_results(dep_somata)

# wide data.frame
df_wide_somata <- get_df_wide(dep_somata)
# long data.frame
df_long_somata <- get_df_long(dep_somata)

# making supplemental table 2
supp_table_2 <- df_wide_somata %>%
  distinct(uniprot_id, FLyesCre__vs_WTyesCre__diff,
           FLyesCre__vs_WTyesCre__p.val,
           FLyesCre__vs_WTyesCre__p.adj,
           ensembl_gene_id,
           external_gene_name) %>%
  rename(Null_over_WT_log2FC = FLyesCre__vs_WTyesCre__diff,
         Null_over_WT_pval = FLyesCre__vs_WTyesCre__p.val,
         Null_over_WT_qval = FLyesCre__vs_WTyesCre__p.adj)
write.csv(supp_table_2, "supp_table_2.csv", quote = F, row.names = F)


## ----wGCprocess--------------------------------------------------------------------------------------------------------------------------
# loading wGC data
data_wGC <- wGC_data %>%
  rename(Sample.Early1 = TotInt_126C_cpnGCs,
         Sample.Early2 = TotInt_127N_cpnGCs,
         Sample.Early3 = TotInt_127C_cpnGCs,
         Sample.Early4 = TotInt_128N_cpnGCs,
         Sample.Late1 = TotInt_128C_cpnGCs,
         Sample.Late2 = TotInt_129N_cpnGCs,
         Sample.Late3 = TotInt_129C_cpnGCs,
         Sample.Late4 = TotInt_130N_cpnGCs,
         Sample.WT1 = TotInt_130C_cpnGCs,
         Sample.WT2 = TotInt_131N_cpnGCs,
         Sample.WT3 = TotInt_131C_cpnGCs,
         Sample.WT4 = TotInt_132N_cpnGCs,
         Sample.WT5 = TotInt_132C_cpnGCs,
         Sample.Null1 = TotInt_133N_cpnGCs,
         #Sample.Null2 = TotInt_133C_cpnGCs, # outlier by pca
         Sample.Null3 = TotInt_134N_cpnGCs,
         Sample.Null4 = TotInt_134C_cpnGCs,
         Sample.Null5 = TotInt_135N_cpnGCs) %>%
  left_join(final_mart %>%
              distinct(external_gene_name, ensembl_gene_id), 
            by = "ensembl_gene_id") %>%
  mutate(uniprot_id = gsub("-.*", "", uniprot_id)) %>%
  left_join(read_excel("external_tables.xlsx", sheet = "wGCs_unmapped") %>%
              rename(uniprot_id = From,
                     external_gene_name = To),
            by = "uniprot_id") %>%
  mutate(external_gene_name = case_when(is.na(external_gene_name.x) == T ~
                                          external_gene_name.y,
                                        is.na(external_gene_name.y) == T ~
                                        external_gene_name.x,
                                        T ~ external_gene_name.x),
        .keep = "unused")
  
# make unique names using external gene names and uniprot ids
data_unique_wGC <- make_unique(data_wGC, "external_gene_name", "uniprot_id", 
                               delim = ";")

# parse condition information from the column names
columns_wGC <- grep("Sample", colnames(data_unique_wGC))
data_se_wGC <- make_se_parse(data_unique_wGC, columns_wGC)

# barplot of the protein identification overlap between samples
plot_frequency(data_se_wGC)

# loose filtering (identified in all - n replicates of at least one condition)
data_filt_wGC <- filter_missval(data_se_wGC, thr = 1)

# barplot of the number of identified proteins per sample
plot_numbers(data_filt_wGC)

# barplot of the protein identification overlap between samples
plot_coverage(data_filt_wGC)

# normalize the data
data_norm_wGC <- normalize_vsn(data_filt_wGC)

# visualize normalization for all samples before and after normalization
plot_normalization(data_filt_wGC, data_norm_wGC)

# impute via quantile regression-based left-censored function (QRLIC)
# data_imp <- impute(data_norm, fun = "QRILC")
data_imp_wGC <- data_norm_wGC

# intensity distributions before and after imputation
plot_imputation(data_filt_wGC, data_norm_wGC)#, data_imp_wGC)

# identify differential genes
data_diff_manual_wGC <- test_diff(data_imp_wGC, type = "manual",
                              test = c("Late_vs_Early",
                                       "Null_vs_WT"))

# denote significant proteins (modify "alpha" and "lfc" as needed)
dep_wGC <- add_rejections(data_diff_manual_wGC, alpha = 0.1, lfc = log2(0))

# pca for top 500 proteins
pca_var_wGC <- apply(assay(dep_wGC), 1, sd)
pca_var_df_wGC <- assay(dep_wGC)[order(pca_var_wGC, 
                       decreasing = TRUE)[seq_len(min(c(nrow(dep_wGC), 500)))],]
pca_wGC <- prcomp(t(pca_var_df_wGC), scale = FALSE)
pca_df_wGC <- pca_wGC$x %>% data.frame() %>% rownames_to_column() %>% 
  left_join(., data.frame(colData(dep_wGC)), by = c(rowname = "ID"))
pca_percent_wGC <- round(100 * pca_wGC$sdev^2/sum(pca_wGC$sdev^2), 1)

# plotting pca for top 500 proteins
pca_gcs <- ggplot(pca_df_wGC %>%
         mutate(Age = case_when(condition != "Early" ~ "P3",
                                T ~ "P1"),
                Bcl11a = case_when(condition != "Null" ~ "WT",
                                   T ~ "Null"),
                Bcl11a = as.factor(Bcl11a),
                Bcl11a = fct_relevel(Bcl11a, "WT", "Null")), 
       aes(get(paste0("PC", 1)), 
           get(paste0("PC", 2)),
           color = Bcl11a,
           fill = Bcl11a,
           shape = Age)) + 
  labs(title = "PC Analysis", 
       x = paste0("PC", 1, ": ", pca_percent_wGC[1], "%"), 
       y = paste0("PC", 2, ": ", pca_percent_wGC[2], "%")) + 
  coord_fixed() +
  theme_classic(base_size = 7) +
  scale_shape_manual(values = c(24, 21)) +
  scale_color_manual(values = c("#009E73", "#CC79A7")) +
  scale_fill_manual(values = c("#009E73", "#CC79A7")) +
  geom_point() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.3)); pca_gcs
ggsave("supp_fig_1k.pdf", pca_gcs, units = "in", width = 2.52, height = 1.77)

# volcano plot for specific contrast
plot_volcano(dep_wGC, contrast = "Late_vs_Early",
             label_size = 2, add_names = TRUE)
plot_volcano(dep_wGC, contrast = "Null_vs_WT",
             label_size = 2, add_names = TRUE)

# frequency plot of significant proteins for the different conditions
plot_cond(dep_wGC)

# results table
data_results_wGC <- get_results(dep_wGC)

# wide data.frame
df_wide_wGC <- get_df_wide(dep_wGC)
# long data.frame
df_long_wGC <- get_df_long(dep_wGC)

# making supplemental table 3
supp_table_3 <- df_wide_wGC %>%
  distinct(uniprot_id,
           Null_vs_WT_diff,
           Null_vs_WT_p.val,
           Null_vs_WT_p.adj,
           Late_vs_Early_diff,
           Late_vs_Early_p.val,
           Late_vs_Early_p.adj,
           ensembl_gene_id,
           external_gene_name) %>%
  rename(Null_over_WT_log2FC = Null_vs_WT_diff,
         Null_over_WT_pval = Null_vs_WT_p.val,
         Null_over_WT_qval = Null_vs_WT_p.adj,
         P3_over_P0_log2FC = Late_vs_Early_diff,
         P3_over_P0_pval = Late_vs_Early_p.val,
         P3_over_P0_qval = Late_vs_Early_p.adj)
write.csv(supp_table_3, "supp_table_3.csv", quote = F, row.names = F)


## ----dataAnalylsis-----------------------------------------------------------------------------------------------------------------------
# loading wGC analysis
wGC_dep <- df_wide_wGC %>%
  select(external_gene_name, uniprot_id, ensembl_gene_id, Early_1:Null_5,
         contains("diff"), contains("p."), contains("significant")) %>%
  mutate(external_gene_name = case_when(is.na(external_gene_name) == T ~
                                          uniprot_id,
                                        T ~ external_gene_name)) %>%
  rowwise() %>%
  mutate(baseMean = mean(c_across(Early_1:Null_5))) %>%
  ungroup() %>%
  mutate(Null_vs_WT_significant = case_when(Null_vs_WT_p.adj < 0.1 &
                                              Null_vs_WT_p.val < 0.1 ~ T,
                                            T ~ F),
         Late_vs_Early_significant = case_when(Late_vs_Early_p.adj < 0.1 &
                                              Late_vs_Early_p.val < 0.1 ~ T,
                                            T ~ F),
         Null_vs_WT_significant = as.factor(Null_vs_WT_significant),
         Late_vs_Early_significant = as.factor(Late_vs_Early_significant),
         Null_vs_WT_significant = fct_relevel(Null_vs_WT_significant, 
                                              "TRUE", "FALSE"),
         Late_vs_Early_significant = fct_relevel(Late_vs_Early_significant, 
                                              "TRUE", "FALSE"),
         bcl11a_label = case_when(Null_vs_WT_significant == T ~
                                     external_gene_name,
                                  external_gene_name %in% c("Lrrtm2", 
                                                               "Sema3f",
                                                               "Spire2",
                                                               "Camsap3",
                                                            "Cdk5") == T ~
                                     external_gene_name, T ~ ""),
         dev_label = case_when(Late_vs_Early_significant == T ~
                                     external_gene_name, 
                               external_gene_name %in% c("Lrrtm2", 
                                                               "Sema3f",
                                                               "Spire2",
                                                               "Camsap3",
                                                            "Cdk5") == T ~
                                     external_gene_name, T ~ ""))

# loading somata analysis
somata_dep <- df_wide_somata %>%
  select(external_gene_name, uniprot_id, ensembl_gene_id, FLyesCre__1:WTyesCre__3,
         contains("diff"), contains("p."), contains("significant")) %>%
  mutate(external_gene_name = case_when(is.na(external_gene_name) == T ~
                                          uniprot_id,
                                        T ~ external_gene_name)) %>%
  rowwise() %>%
  mutate(baseMean = mean(c_across(FLyesCre__1:WTyesCre__3))) %>%
  ungroup() %>%
  mutate(FLyesCre__vs_WTyesCre__significant = case_when(FLyesCre__vs_WTyesCre__p.adj < 0.1 &
                                                          FLyesCre__vs_WTyesCre__p.val < 0.1 ~ T,
                                                        T ~ F),
         FLyesCre__vs_WTyesCre__significant = as.factor(FLyesCre__vs_WTyesCre__significant),
         FLyesCre__vs_WTyesCre__significant = fct_relevel(FLyesCre__vs_WTyesCre__significant, 
                                                          "TRUE", "FALSE"),
         pythag = sqrt(FLyesCre__vs_WTyesCre__diff * FLyesCre__vs_WTyesCre__diff +
                         -log10(FLyesCre__vs_WTyesCre__p.val) * -log10(FLyesCre__vs_WTyesCre__p.val)),
         pythag_label = case_when(pythag > 4.7 ~ external_gene_name,
                                  external_gene_name %in% c("Lrrtm2", 
                                                            "Sema3f",
                                                            "Spire2",
                                                            "Camsap3",
                                                            "Cdk5") == T ~
                                    external_gene_name,
                                  T ~ ""),
         gc_label = case_when(external_gene_name %in% c("Lrrtm2", 
                                                            "Sema3f",
                                                            "Spire2",
                                                            "Camsap3",
                                                            "Cdk5",
                                                        "Ptpn2", "Spag9",
                                                        "Akap12", "Camk2b",
                                                        "Marcks", "Cnn3", 
                                                        "Cotl1", "Sec61b",
                                                        "Ctnnb1", "Jmy",
                                                        "Bcl11a") ~ 
                                external_gene_name,
                              T ~ ""))

# volcano plots
somata_volcano <- ggplot(somata_dep,
                             aes(x = FLyesCre__vs_WTyesCre__diff,
                                 y = -log10(FLyesCre__vs_WTyesCre__p.val),
                                 label = gc_label,
                                 color = FLyesCre__vs_WTyesCre__significant)) +
  labs(x = "log2(Fold Change)",
       y = "-log10(p-value)",
       color = "q < 0.1") +
  scale_color_manual(values = c("black", "grey50")) +
  geom_point(size = 0.000000000000000000000000000000000000000000000000000001) +
  theme_classic(base_size = 7, base_family = "sans") +
  geom_text_repel(size = 2, max.overlaps = 100) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = c(0.9, 0.85)); somata_volcano
ggsave("fig_1b.pdf", somata_volcano, units = "in", width = 2.20, height = 2.20)

wGC_volcano_bcl11a <- ggplot(wGC_dep,
                             aes(x = Null_vs_WT_diff,
                                 y = -log10(Null_vs_WT_p.val),
                                 color = Null_vs_WT_significant,
                                 label = bcl11a_label)) +
  labs(x = "log2(Fold Change)",
       y = "-log10(p-value)",
       color = "q < 0.1") +
  scale_color_manual(values = c("black", "grey50")) +
  geom_point(size = 0.000000000000000000000000000000000000000000000000000001) +
  theme_classic(base_size = 7) +
  geom_text_repel(size = 2, max.overlaps = 100) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = c(0.9, 0.85)); wGC_volcano_bcl11a
ggsave("fig_1c.pdf", wGC_volcano_bcl11a, units = "in", width = 2.20, height = 2.20)

wGC_volcano_dev <- ggplot(wGC_dep,
                             aes(x = Late_vs_Early_diff,
                                 y = -log10(Late_vs_Early_p.val),
                                 color = Late_vs_Early_significant,
                                 label = dev_label)) +
  labs(x = "log2(Fold Change)",
       y = "-log10(p-value)",
       color = "q < 0.1") +
  scale_color_manual(values = c("black", "grey50")) +
  geom_point(size = 0.000000000000000000000000000000000000000000000000000001) +
  theme_classic(base_size = 7) +
  geom_text_repel(size = 2, max.overlaps = 100) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = c(0.9, 0.85)); wGC_volcano_dev
ggsave("supp_fig_1h.pdf", wGC_volcano_dev, units = "in", width = 1.91, height = 1.85)

# protein quadrant plot dataframes
quad_bcl11a_gc_pro <- wGC_dep %>%
  select(external_gene_name, uniprot_id, ensembl_gene_id, contains("Null_vs_WT")) %>%
  rename(log2FoldChange = "Null_vs_WT_diff",
         pvalue = "Null_vs_WT_p.val",
         padj = "Null_vs_WT_p.adj",
         significant = "Null_vs_WT_significant")

quad_bcl11a_soma_pro <- somata_dep %>%
  select(external_gene_name, uniprot_id, ensembl_gene_id, contains("FLyesCre__vs_WTyesCre_")) %>%
  rename(log2FoldChange = "FLyesCre__vs_WTyesCre__diff",
         pvalue = "FLyesCre__vs_WTyesCre__p.val",
         padj = "FLyesCre__vs_WTyesCre__p.adj",
         significant = "FLyesCre__vs_WTyesCre__significant")

quad_dev_gc_pro <- wGC_dep %>%
  select(external_gene_name, uniprot_id, ensembl_gene_id, contains("Late_vs_Early")) %>%
  rename(log2FoldChange = "Late_vs_Early_diff",
         pvalue = "Late_vs_Early_p.val",
         padj = "Late_vs_Early_p.adj",
         significant = "Late_vs_Early_significant")

# loading rna data
durak_degs <- read_excel("external_tables.xlsx", sheet = "durak_degs") %>%
  rename(ensembl_gene_id = ensembl_id)
priya_degs <- read_excel("external_tables.xlsx", sheet = "priya_degs") %>%
  group_by(ensembl_gene_id) %>%
  slice(which.min(padj_P3_DivideBy_P1.soma)) %>%
  ungroup() %>%
  rename(pvalue = pvalue_P3_DivideBy_P1.soma,
         padj = padj_P3_DivideBy_P1.soma,
         log2FoldChange = log2FoldChange_P3_DivideBy_P1.soma)

# bcl11a soma protein vs bcl11a soma rna quadrant plot
quad_somaPro_somaRNA <- ggplot(quad_bcl11a_soma_pro %>%
         full_join(durak_degs,
                   by = "ensembl_gene_id") %>%
         filter(is.na(log2FoldChange.x) == F & is.na(log2FoldChange.y) == F) %>%
           mutate(first_sig_color = case_when(padj.x < 0.1 & pvalue.x < 0.1 ~ "Soma Pro",
                                              T ~ ""),
                  second_sig_color = case_when(padj.y < 0.1 & pvalue.y < 0.1 ~ "Soma RNA",
                                               T ~ ""),
                  sig_color = case_when(first_sig_color == "Soma Pro" & 
                                          second_sig_color == "Soma RNA" ~ "Both",
                                        T ~ paste0(first_sig_color, second_sig_color))) %>%
         filter(sig_color != "") %>%
         mutate(sig_color = as.factor(sig_color),
                sig_color = fct_relevel(sig_color, "Soma Pro", "Soma RNA", "Both")),
       aes(x = log2FoldChange.y, y = log2FoldChange.x, color = sig_color)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "black") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "black") +
  geom_point(size = 0.0000000000000000000000000000000000000000000000000000001) +
  theme_classic(base_size = 7) +
  labs(x = "Soma RNA: log2FC(Null/WT)",
       y = "Soma Pro: log2FC(Null/WT)",
       color = "Significance") +
  scale_color_manual(values = c("#D55E00", "#0072B2", "purple")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.1, 0.15)) +
  scale_x_continuous(sec.axis = dup_axis(name = "", 
                                           labels = NULL, breaks = 0)) +
  scale_y_continuous(sec.axis = dup_axis(name = "", 
                                         labels = NULL, breaks = 0)); quad_somaPro_somaRNA
ggsave("fig_1d.pdf", quad_somaPro_somaRNA, units = "in", width = 2.20, height = 2.20)

# bcl11a gc protein vs bcl11a soma rna quadrant plot
quad_gcPro_somaRNA <- ggplot(quad_bcl11a_gc_pro %>%
         full_join(durak_degs,
                   by = "ensembl_gene_id") %>%
         filter(is.na(log2FoldChange.x) == F & is.na(log2FoldChange.y) == F) %>%
         mutate(first_sig_color = case_when(padj.x < 0.1 & pvalue.x < 0.1 ~ "GC Pro",
                                              T ~ ""),
                  second_sig_color = case_when(padj.y < 0.1 & pvalue.y < 0.1 ~ "Soma RNA",
                                               T ~ ""),
                  sig_color = case_when(first_sig_color == "GC Pro" & 
                                          second_sig_color == "Soma RNA" ~ "Both",
                                        T ~ paste0(first_sig_color, second_sig_color))) %>%
         filter(sig_color != "") %>%
         mutate(sig_color = as.factor(sig_color),
                sig_color = fct_relevel(sig_color, "GC Pro", "Soma RNA", "Both")),
       aes(x = log2FoldChange.y, y = log2FoldChange.x, color = sig_color)) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "black") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "black") +
  geom_point(size = 0.0000000000000000000000000000000000000000000000000000001) +
  theme_classic(base_size = 7) +
  labs(x = "Soma RNA: log2FC(Null/WT)",
       y = "GC Pro: log2FC(Null/WT)",
       color = "Significance") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.1, 0.85)) +
  scale_color_manual(values = c("#D55E00", "#0072B2", "purple")) +
  scale_x_continuous(sec.axis = dup_axis(name = "", 
                                           labels = NULL, breaks = 0)) +
  scale_y_continuous(sec.axis = dup_axis(name = "", 
                                         labels = NULL, breaks = 0)); quad_gcPro_somaRNA
ggsave("fig_1e.pdf", quad_gcPro_somaRNA, units = "in", width = 2.20, height = 2.20)

# dev gc protein vs dev soma rna quadrant plot
quad_gcPro_somaRNA_dev <- ggplot(quad_dev_gc_pro %>%
         full_join(priya_degs,
                   by = "ensembl_gene_id") %>%
         filter(is.na(log2FoldChange.x) == F & is.na(log2FoldChange.y) == F) %>%
         mutate(first_sig_color = case_when(padj.x < 0.1 & pvalue.x < 0.1 ~ "GC Pro",
                                              T ~ ""),
                  second_sig_color = case_when(padj.y < 0.1 & pvalue.y < 0.1 ~ "Soma RNA",
                                               T ~ ""),
                  sig_color = case_when(first_sig_color == "GC Pro" & 
                                          second_sig_color == "Soma RNA" ~ "Both",
                                        T ~ paste0(first_sig_color, second_sig_color))) %>%
         filter(sig_color != "") %>%
         mutate(sig_color = as.factor(sig_color),
                sig_color = fct_relevel(sig_color, "GC Pro", "Soma RNA", "Both")),
       aes(x = log2FoldChange.y, y = log2FoldChange.x, color = sig_color)) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "black") +
  geom_vline(xintercept = 0, linewidth = 0.25, color = "black") +
  geom_point(size = 0.0000000000000000000000000000000000000000000000000000001) +
  theme_classic(base_size = 7) +
  labs(x = "Soma RNA: log2FC(P3/P1)",
       y = "GC Pro: log2FC(P3/P0)",
       color = "Significance") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.1, 0.85)) +
  scale_color_manual(values = c("#D55E00", "#0072B2", "purple")) +
  scale_x_continuous(sec.axis = dup_axis(name = "", 
                                           labels = NULL, breaks = 0)) +
  scale_y_continuous(sec.axis = dup_axis(name = "", 
                                         labels = NULL, breaks = 0)); quad_gcPro_somaRNA_dev
ggsave("supp_fig_1i.pdf", quad_gcPro_somaRNA_dev, units = "in", width = 1.91, height = 1.85)



## ----signalPeptidesBcl11a----------------------------------------------------------------------------------------------------------------
# loading soma go terms for rnas depleted after bcl11a depletion
bcl11a_go <- read_excel("external_tables.xlsx", sheet = "durak_go") %>%
  arrange(qvalue)

# plotting soma go terms
bcl11a_go_plot <- ggplot(bcl11a_go, aes(x = GeneRatio, y = Description, 
                      fill = qvalue, size = Num)) +
  geom_point(color = "black", pch = 21) +
  labs(y = "", title = "GO Analysis: Depleted Soma RNAs After Bcl11a Deletion") +
  theme_classic(base_size = 7, base_family = "sans") +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); bcl11a_go_plot
ggsave("supp_fig_3a.pdf", bcl11a_go_plot, units = "in",
       width = 6.72, height = 2.77)


## ----sessionInfo-------------------------------------------------------------------------------------------------------------------------
sessionInfo()

