## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(readxl)
library(ggpubr)
library(biomaRt)
library(org.Mm.eg.db)
library(clusterProfiler)
library(DEP)
library(SummarizedExperiment)
library(ggrepel)
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


## ----slFSPS------------------------------------------------------------------------------------------------------------------------------
# loading slfsps data
slfsps_data <- read_excel("input_tables.xlsx", sheet = "Input Table 10") %>%
  
  # using gates from flowjo
  mutate(P1 = case_when(FSC.PMT.A < 40000 & SSC.A < 40000 ~ T, T ~ F),
         GFP = case_when(GFP.A > 90 ~ T, T ~ F),
         RFP = case_when(tdTomato.A > 30 ~ T, T ~ F),
         Construct = case_when(Construct == "pDT016e" ~ "cytLrrtm2",
                               Construct == "pDT016es" ~ "memLrrtm2"))

# focusing on cpn/nrxn double-positive gcs
slfsps_cpn <- slfsps_data %>%
  filter(P1 == T & RFP == T & GFP == T)
slfsps_lrrtm2 <- ggplot(slfsps_cpn, aes(x = Construct, y = tdTomato.A, fill = Construct)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), 
              linewidth = 0.5, color = "black") +
  stat_compare_means(comparisons = list(c("cytLrrtm2", "memLrrtm2")),
                     method = "t.test", label = "p.signif",
                     method.args = list(alternative = "two.sided"),
                     size = 2) +
  scale_y_continuous(trans = "log10") +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#56B4E9", "#0072B2")) +
  labs(x = "", y = "Nrxn Fluoresence") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none"); slfsps_lrrtm2
ggsave("fig_4c.pdf", slfsps_lrrtm2, units = "in",
       height = 1.19, width = 1.57)

# finding representative images
test_means <- slfsps_cpn %>%
  group_by(NB, Construct) %>%
  summarize(mean_red = mean(tdTomato.A))


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


## ----inSilicoIP--------------------------------------------------------------------------------------------------------------------------
# loading alpha pulldown data
data_ap <- read_excel("input_tables.xlsx", sheet = "Input Table 11") %>%
  separate(combination, into = c("Lrrtm2", "and", "uniprot_id")) %>%
  select(uniprot_id, iptm, pDockQ) %>%
  
  # making z scores to use iptm and docking scores
  mutate(iptm_z = (iptm - mean(iptm, na.rm = T)) / sd(iptm, na.rm = T),
         dock_z = (pDockQ - mean(pDockQ, na.rm = T)) / sd(pDockQ, na.rm = T),
         z_mean = (iptm_z + dock_z) / 2) %>%
  drop_na() %>%
  left_join(final_mart %>%
              distinct(uniprot_id, external_gene_name), by = "uniprot_id") %>%
  mutate(external_gene_name = case_when(is.na(external_gene_name) == T ~ uniprot_id,
                                        T ~ external_gene_name))

# making supplemental table 4
supp_table_4 <- data_ap %>%
  distinct(uniprot_id, iptm, pDockQ, external_gene_name)
write.csv(supp_table_4, "supp_table_4.csv", quote = F, row.names = F)

# # gsea analysis of alphapulldown hits
# data_ap_uniprot <- data_ap %>%
#   select(-external_gene_name) %>%
#   distinct()
# list_iptm_dock <- data_ap_uniprot %>%
#   select(uniprot_id, z_mean) %>%
#   arrange(desc(z_mean))
# geneList_iptm_dock <- list_iptm_dock[[2]]
# names(geneList_iptm_dock) <- as.character(list_iptm_dock[[1]])
# geneList_iptm_dock <- sort(geneList_iptm_dock, decreasing = T)
# set.seed(440) # setting random seed for consistent gsea tiebreaking
# gse_iptm_dock <- gseGO(geneList = geneList_iptm_dock,
#                        ont = "ALL",
#                        OrgDb = org.Mm.eg.db,
#                        keyType = "UNIPROT",
#                        eps = 0) %>%
#   as.data.frame()

# # removing redundant gsea terms
# plot_gse_df_iptm_dock <- gse_iptm_dock %>%
#   as.data.frame() %>%
#   arrange(qvalue) %>%
#   mutate(trans_qvalues = -log10(qvalue))
# gse_df_preReduce_iptm_dock <- plot_gse_df_iptm_dock %>%
#   select(ONTOLOGY, ID) %>%
#   rename(go_type = ONTOLOGY,
#            go_id = ID)
# gse_df_scores_iptm_dock <- plot_gse_df_iptm_dock$trans_qvalue %>%
#   setNames(plot_gse_df_iptm_dock$ID)
# gse_df_reduce_iptm_dock <- go_reduce(gse_df_preReduce_iptm_dock,
#                           orgdb = "org.Mm.eg.db",
#                           threshold = 0.7,
#                           scores = gse_df_scores_iptm_dock,
#                           measure = "Wang")
# plot_gse_iptm_dock_reduce <- gse_iptm_dock %>%
#   filter(ID %in% gse_df_reduce_iptm_dock$parent_id)

# # saving locally for quicker access
# write.csv(plot_gse_iptm_dock_reduce,
#           "misc_noGitHub/figure_5/gse_alphaPulldown_reduce.csv",
#           row.names = F)

gse_df_reduce_iptm_dock <- read.csv("misc_noGitHub/figure_5/gse_alphaPulldown_reduce.csv") %>%
  mutate(gene_num = str_count(core_enrichment, "\\/") + 1,
         geneRatio = gene_num/setSize) %>%
  filter(enrichmentScore > 0)

# plotting top gsea mf terms
gsea_df_reduce_plot <- ggplot(gse_df_reduce_iptm_dock %>%
                                filter(ONTOLOGY == "MF") %>%
                                arrange(qvalue) %>%
                                slice(1:15) %>%
                                mutate(Description = case_when(Description == 
                                                                 "catalytic activity, acting on RNA" ~ 
                                                                 "catalytic activity on RNA",
                                                               Description == 
                                                                 "catalytic activity, acting on a nucleic acid" ~
                                                                 "catalytic activity on nucleic acid",
                                                               Description == 
                                                                 "hydrolase activity, acting on ester bonds" ~ 
                                                                 "hydrolase activity on ester bonds",
                                                               Description == 
                                                                 "oxidoreductase activity, acting on CH-OH group of donors" ~ 
                                                                 "oxidoreductase activity on CH-OH donors",
                                                               T ~ Description)) %>%
                                mutate(Description = fct_reorder2(Description,
                                                                  -qvalue, 
                                                                  qvalue)),
                              aes(x = geneRatio, y = Description,
                                  fill = qvalue, size = gene_num)) +
  geom_point(color = "black", pch = 21) +
  scale_size(range = c(0.25, 2.5)) +
  labs(y = "", title = "GSEA Analysis: In Silico IP",
       size = "Number") +
  theme_classic(base_size = 7, base_family = "sans") +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); gsea_df_reduce_plot
ggsave("fig_5d.pdf", gsea_df_reduce_plot, units = "in", 
       height = 2.22, width = 3.89)


## ----inVivoIP----------------------------------------------------------------------------------------------------------------------------
# loading ip data
data_ip <- read_excel("input_tables.xlsx", sheet = "Input Table 12") %>%
  
  rename(Sample.IgG1 = Corrected_060725A_Dustin_Strap1,
         Sample.IgG2 = Corrected_060725A_Dustin_Strap2,
         Sample.IgG3 = Corrected_060725A_Dustin_Strap3,

         Sample.IgG4 = Corrected_060725A_Dustin_Strap4,
         Sample.IgG5 = Corrected_060725A_Dustin_Strap5,
         Sample.IgG6 = Corrected_060725A_Dustin_Strap6,

         Sample.MemIP1 = Corrected_060725A_Dustin_Strap7,
         Sample.MemIP2 = Corrected_060725A_Dustin_Strap8,
         Sample.MemIP3 = Corrected_060725A_Dustin_Strap9,

         Sample.CytIP1 = Corrected_060725A_Dustin_Strap10,
         Sample.CytIP2 = Corrected_060725A_Dustin_Strap11,
         Sample.CytIP3 = Corrected_060725A_Dustin_Strap12) %>%
  
  filter(is.na(Filter) == T) %>%
  filter(!if_any(everything(), ~ grepl("FGCZ", .))) %>%
  separate(Accession, into = c("source", "Protein.IDs", "temp"), 
           sep = "\\|") %>%
  filter(!(if_all(contains("Sample"), ~ . == 0))) %>%
  separate(temp, into = c("Gene.Names", "dump"), sep = "_") %>%
  select(-dump, -Filter) %>%
  mutate(across(contains("Sample"), ~ na_if(., 0)))

# make unique names using external gene names and uniprot ids
data_unique_ip <- make_unique(data_ip, "Gene.Names", "Protein.IDs", 
                               delim = ";")

# parse condition information from the column names
columns_ip <- grep("Sample", colnames(data_unique_ip))
data_se_ip <- make_se_parse(data_unique_ip, columns_ip)

# barplot of the protein identification overlap between samples
plot_frequency(data_se_ip)

# loose filtering (identified in all - n replicates of at least one condition)
data_filt_ip <- filter_missval(data_se_ip, thr = 2)

# barplot of the number of identified proteins per sample
plot_numbers(data_filt_ip)

# barplot of the protein identification overlap between samples
plot_coverage(data_filt_ip)

# normalize the data
data_norm_ip <- normalize_vsn(data_filt_ip)

# visualize normalization for all samples before and after normalization
plot_normalization(data_filt_ip, data_norm_ip)

# heatmap of proteins with missing values
plot_missval(data_filt_ip)

# intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt_ip)

# impute via quantile regression-based left-censored function (QRLIC)
set.seed(311) # setting random seed for consistent imputation draws
data_imp_ip <- impute(data_norm_ip, fun = "QRILC")

# intensity distributions before and after imputation
plot_imputation(data_filt_ip, data_norm_ip, data_imp_ip)

# identify differential genes
data_diff_manual_ip <- test_diff(data_imp_ip, type = "manual",
                                 test = c("IgG_vs_MemIP",
                                          "IgG_vs_CytIP"))

# denote significant proteins (modify "alpha" and "lfc" as needed)
dep_ip <- add_rejections(data_diff_manual_ip, alpha = 0.1, lfc = log2(0))

# pca for top 500 proteins
pca_var_ip <- apply(assay(dep_ip), 1, sd)
pca_var_df_ip <- assay(dep_ip)[order(pca_var_ip,
                       decreasing = TRUE)[seq_len(min(c(nrow(dep_ip), 500)))],]
pca_ip <- prcomp(t(pca_var_df_ip), scale = FALSE)
pca_df_ip <- pca_ip$x %>% data.frame() %>% rownames_to_column() %>%
  left_join(., data.frame(colData(dep_ip)), by = c(rowname = "ID"))
pca_percent_ip <- round(100 * pca_ip$sdev^2/sum(pca_ip$sdev^2), 1)

# plotting pca for top n proteins
pca_ip <- ggplot(pca_df_ip %>%
                   mutate(condition = case_when(condition == "CytIP" ~ 
                                                  "cytLrrtm2",
                                                condition == "MemIP" ~ 
                                                  "memLrrtm2",
                                                condition == "IgG" ~ 
                                                  "IgG"),
                          condition = fct_relevel(as.factor(condition), 
                                                  "cytLrrtm2", 
                                                  "memLrrtm2", 
                                                  "IgG")),
       aes(get(paste0("PC", 1)),
           get(paste0("PC", 2)),
           color = condition,
           fill = condition)) +
  labs(title = "In Vivo IP: PC Analysis",
       x = paste0("PC", 1, ": ", pca_percent_ip[1], "%"),
       y = paste0("PC", 2, ": ", pca_percent_ip[2], "%"),
       color = "", fill = "") +
  coord_fixed() +
  theme_classic(base_size = 7) +
  scale_color_manual(values = c("#56b4e9", "#0072b2", "violet")) +
  scale_fill_manual(values = c("#56b4e9", "#0072b2", "violet")) +
  geom_point(size = 0.25) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.8, 0.3)); pca_ip
ggsave("supp_fig_7r.pdf", pca_ip, units = "in", width = 2, height = 1.66)

# volcano plot for specific contrast
plot_volcano(dep_ip, contrast = "IgG_vs_MemIP",
             label_size = 2, add_names = TRUE)
plot_volcano(dep_ip, contrast = "IgG_vs_CytIP",
             label_size = 2, add_names = TRUE)

# frequency plot of significant proteins for the different conditions
plot_cond(dep_ip)

# results table
data_results_ip <- get_results(dep_ip)

# wide data.frame
df_wide_ip <- get_df_wide(dep_ip)

# long data.frame
df_long_ip <- get_df_long(dep_ip)

# making data.frame for volcano plots
df_volcano <- df_wide_ip %>%
  select(Protein.IDs, contains("IgG_vs_")) %>%
  select(-contains("CI.")) %>%
  mutate(IgG_vs_CytIP_significant = case_when(IgG_vs_CytIP_p.val < 0.1 &
                                                IgG_vs_CytIP_p.adj < 0.1 ~ T,
                                              T ~ F),
         IgG_vs_MemIP_significant = case_when(IgG_vs_MemIP_p.val < 0.1 &
                                                IgG_vs_MemIP_p.adj < 0.1 ~ T,
                                              T ~ F),
         
         CytIP_vs_IgG_diff = -IgG_vs_CytIP_diff,
         MemIP_vs_IgG_diff = -IgG_vs_MemIP_diff,
         
         IgG_vs_CytIP_significant = as.factor(IgG_vs_CytIP_significant),
         IgG_vs_MemIP_significant = as.factor(IgG_vs_MemIP_significant),
         IgG_vs_CytIP_significant = fct_relevel(IgG_vs_CytIP_significant, 
                                              "TRUE", "FALSE"),
         IgG_vs_MemIP_significant = fct_relevel(IgG_vs_MemIP_significant, 
                                              "TRUE", "FALSE")) %>%
  rename(uniprot_id = Protein.IDs) %>%
  left_join(final_mart %>% distinct(uniprot_id, external_gene_name),
            by = "uniprot_id") %>%
  mutate(uniprot_id = gsub("-.*", "", uniprot_id)) %>%
  left_join(read_excel("external_tables.xlsx", sheet = "ip_unmapped") %>%
              rename(uniprot_id = From,
                     external_gene_name = To),
            by = "uniprot_id") %>%
  mutate(external_gene_name = case_when(is.na(external_gene_name.x) == T ~
                                          external_gene_name.y,
                                        is.na(external_gene_name.y) == T ~
                                        external_gene_name.x,
                                        T ~ external_gene_name.x),
        .keep = "unused") %>%
  mutate(cyt_label = case_when(external_gene_name %in% c("Pcdh15", "Cyth2",
                                                         "Dock3", "Dhx9",
                                                         "Hnrnpa0", "Hnrnpl",
                                                         "Cbx1", "Hmgn2",
                                                         "Mybbp1a", "Prdm14",
                                                         "Safb", "Sltm",
                                                         "Igkv4-57", "Nrxn2",
                                                         "Sec31b") ~ 
                                 external_gene_name, T ~ ""),
         mem_label = case_when(external_gene_name %in% c("Pcdhac1", "Dock3",
                                                         "Sec31b", "Dppa3",
                                                         "Hmgn2", "Ints4",
                                                         "Igkv4-57", "Anxa1",
                                                         "Nrxn2") ~ 
                                 external_gene_name, T ~ ""))

# cytLrrmtm2 IP volcano plot
cyt_volcano <- ggplot(df_volcano,
                             aes(x = CytIP_vs_IgG_diff,
                                 y = -log10(IgG_vs_CytIP_p.val),
                                 color = IgG_vs_CytIP_significant,
                                 label = cyt_label)) +
  labs(x = "cytLrrtm2-V5 In Vivo IP: log2(V5 - IgG)",
       y = "-log10(p-value)",
       color = "q < 0.1") +
  scale_color_manual(values = c("black", "grey50")) +
  geom_point(size = 0.000000000000000000000000000000000000000000000000000001) +
  theme_classic(base_size = 7) +
  geom_text_repel(size = 2, max.overlaps = 100) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = c(0.9, 0.85)); cyt_volcano
ggsave("fig_5i.pdf", cyt_volcano, units = "in", width = 2.47, height = 1.91)

# memLrrtm2 IP volcano plot
mem_volcano <- ggplot(df_volcano,
                             aes(x = MemIP_vs_IgG_diff,
                                 y = -log10(IgG_vs_MemIP_p.val),
                                 color = IgG_vs_MemIP_significant,
                                 label = mem_label)) +
  labs(x = "memLrrtm2-V5 In Vivo IP: log2(V5 - IgG)",
       y = "-log10(p-value)",
       color = "q < 0.1") +
  scale_color_manual(values = c("black", "grey50")) +
  geom_point(size = 0.000000000000000000000000000000000000000000000000000001) +
  theme_classic(base_size = 7) +
  geom_text_repel(size = 2, max.overlaps = 100) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = c(0.9, 0.85)); mem_volcano
ggsave("fig_5j.pdf", mem_volcano, units = "in", width = 2.47, height = 1.91)

# making supplemental table 5
supp_table_5 <- df_volcano %>%
  distinct(uniprot_id,
           CytIP_vs_IgG_diff,
           IgG_vs_CytIP_p.val,
           IgG_vs_CytIP_p.adj,
           MemIP_vs_IgG_diff,
           IgG_vs_MemIP_p.val,
           IgG_vs_MemIP_p.adj,
           external_gene_name) %>%
  rename(Cyt_over_IgG_log2FC = CytIP_vs_IgG_diff,
         Cyt_over_IgG_pval = IgG_vs_CytIP_p.val,
         Cyt_over_IgG_qval = IgG_vs_CytIP_p.adj,
         Mem_over_IgG_log2FC = MemIP_vs_IgG_diff,
         Mem_over_IgG_pval = IgG_vs_MemIP_p.val,
         Mem_over_IgG_qval = IgG_vs_MemIP_p.adj)
write.csv(supp_table_5, "supp_table_5.csv", quote = F, row.names = F)


## ----sessionInfo-------------------------------------------------------------------------------------------------------------------------
sessionInfo()

