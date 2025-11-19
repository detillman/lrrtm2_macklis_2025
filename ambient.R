## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(readxl)
library(ggpubr)
library(biomaRt)
library(Mus.musculus)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(UpSetR)
library(clusterProfiler)
library(ggridges)
library(tidyverse)
cbPalette <- c("black", # white 
               "#E69F00", # orange
               "#56B4E9", # light blue
               "#009E73", # green
               "#CC79A7", # pink
               "#D55E00", # red
               "#0072B2", # dark blue
               "#F0E442") # yellow

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


## ----loadingBcl11aSNPs-------------------------------------------------------------------------------------------------------------------
# making metadataframe for bcl11a snp data
metadata_samples <- data.frame(sample = c(12, 13, 14, 
                                              15, 16, 17, 18, 
                                              19, 20, 21, 22),
                                   genotype = c("wt", "wt", "wt", 
                                                "het", "het","het", "het", 
                                                "null", "null", "null", "null"),
                                   replicate = c(1, 2, 3, 
                                                 1, 2, 3, 4, 
                                                 1, 2, 3, 4),
                                   percentLabeledMass = c(0.6470588, 
                                                          0.6470588, 
                                                          0.4285714,
                                                          0.4, 0.4, 
                                                          0.4230769, 
                                                          0.3773585,
                                                          0.7692308, 
                                                          0.4705882, 
                                                          0.4705882, 
                                                          0.4375),
                                   percentCD1Mass = c(0.3529412, 
                                                      0.3529412, 
                                                      0, 0, 0, 0, 0, 
                                                      0.2307692, 
                                                      0.5294118, 
                                                      0.5294118, 
                                                      0.5625))
  
# loading data for bcl11a snps
samples <- read_excel("input_tables.xlsx", sheet = "Input Table 13",
                               col_types = c("text", "numeric", "text", "text",
                                             "text", "text", "numeric")) %>%
  separate("AD", c("RC", "AC_1", "AC_2"), sep = ",", 
           convert = T, fill = "right") %>%
  mutate(AC_1_0 = case_when(is.na(AC_1) == T ~ 0, T ~ AC_1),
         AC_2_0 = case_when(is.na(AC_2) == T ~ 0, T ~ AC_2), 
         .keep = "unused") %>%
  mutate(AC = AC_1_0 + AC_2_0, .keep = "unused") %>%
  mutate(TC = RC + AC, percentAC = AC/TC) %>%
  rename(sample = SAMPLE) %>%
  unite("location", "CHROM", "POS", sep = ".", remove = F) %>%
  left_join(metadata_samples, by = "sample") %>%
  filter(CHROM != "JH584304.1", CHROM != "JH584295.1") %>%
  mutate(location = gsub("MT", "M", location))


## ----vizAllSNPs--------------------------------------------------------------------------------------------------------------------------
# distribution of total read counts for SNPs
total_snp_dist <- ggplot(samples, aes(x = TC)) + 
  geom_density(color = "black", linewidth = 1) + 
  theme_classic(base_size = 10) + scale_x_continuous(trans = "log10") + 
  labs(x = "Total Read Counts") + 
  geom_vline(xintercept = 8.58, color = "black", linetype = "dashed"); total_snp_dist
# total_snp_dist_values <- ggplot_build(total_snp_dist)
# write.csv(total_snp_dist_values$data[[1]], 
#           file = "total_snp_dist_values.csv")

# setting a read filter of > 8 total counts
high_samples <- samples %>%  
  filter(TC > 8)

# gathering SNP information
snp_number <- high_samples %>%
  group_by(sample) %>%
  summarize(snp_number = n(),
            genotype = unique(genotype),
            percentCD1Mass = unique(percentCD1Mass),
            percentLabeledMass = unique(percentLabeledMass))

# number of SNPs across samples (higher in CD1-supplemented samples)
ggplot(snp_number, 
       aes(x = sample, y = snp_number)) + 
  geom_point() + theme_classic(base_size = 10) + 
  labs(x = "Sample Number", y = "SNP Number") +
  theme(legend.position = c(0.15, 0.85))

# number of SNPs directly correlates with ratio of CD1 supplementation
ggplot(snp_number, aes(x = 100*percentCD1Mass, y = snp_number)) +
  scale_color_manual(values = cbPalette) +
  geom_point() + theme_classic(base_size = 10) +
  geom_smooth(method='lm', formula = 'y ~ x', color = "black") +
  stat_regline_equation(label.y = 19000, color = "black", 
                        aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 18000, color = "black",
                        aes(label = ..rr.label..)) +
  labs(x = "% CD1 Mass", y = "SNP Number") +
  theme(legend.position = c(0.15, 0.85))


## ----B6SNPremoval------------------------------------------------------------------------------------------------------------------------
# grouping samples by SNP location and finding average penetrance
pre_adj_samples <- high_samples %>%
  pivot_wider(names_from = sample, values_from = percentAC) %>%
  group_by(location) %>%
  summarize(mean_noSupp = mean(c(`14`, `15`, `16`, `17`, `18`), na.rm = T),
            mean_Supp = mean(c(`12`, `13`, `19`, `20`, `21`, `22`),
                             na.rm = T),
            presence = case_when(is.na(mean_noSupp) == F & 
                                   is.na(mean_Supp) == F ~ "both",
                                 is.na(mean_noSupp) == T & 
                                   is.na(mean_Supp) == F ~ "supp_only",
                                 is.na(mean_noSupp) == F &
                                   is.na(mean_Supp) == T ~ "noSupp_only"))

# finding effect of supplementation on SNP penetrance
pre_adj_samples[is.na(pre_adj_samples) == T] <- 0
pre_adj_samples_noZeroes <- pre_adj_samples %>%
  mutate(supp_diff = mean_Supp - mean_noSupp)
ggplot(pre_adj_samples_noZeroes, aes(x = supp_diff, color = presence)) + 
  geom_density(size = 1) + theme_classic(base_size = 10) + 
  labs(x = expression(bold(Delta*"SNP Penetrance (CD1 Supp - No CD1 Supp)")))

# summarizing SNP counts based on supplement status
pre_count_adj_samples <- pre_adj_samples_noZeroes %>%
  group_by(presence) %>%
  summarize(number = n())
pre_count_adj_samples

# only keeping SNPs that are not identified without CD1 supplement
keep_snps <- pre_adj_samples_noZeroes %>%
  filter(presence == "supp_only")
adj_samples <- high_samples %>%
  filter(percentCD1Mass > 0, location %in% keep_snps$location)


## ----vizCD1SNPs--------------------------------------------------------------------------------------------------------------------------
# gathering information for true SNPs
adj_snp_number <- adj_samples %>%
  group_by(sample) %>%
  summarize(adj_snp_number = n(),
            genotype = unique(genotype),
            percentCD1Mass = unique(percentCD1Mass),
            percentLabeledMass = unique(percentLabeledMass),
            mean_percentAC = mean(percentAC),
            med_percentAC = median(percentAC))

# including zeroes for non-snp samples
adj_snp_number_all <- data.frame(sample = c(14, 15, 16, 17, 18),
                                 adj_snp_number = c(0, 0, 0, 0, 0),
                                 genotype = c("wt", "het", "het", "het", "het"),
                                 percentCD1Mass = c(0, 0, 0, 0, 0),
                                 percentLabeledMass = c(1, 1, 1, 1, 1),
                                 mean_percentAC = c(0, 0, 0, 0, 0),
                                 med_percentAC = c(0, 0, 0, 0, 0)) %>%
  bind_rows(adj_snp_number)


# SNPs are now only present in samples with CD1 supplement
ggplot(adj_snp_number_all, 
       aes(x = sample, y = adj_snp_number, 
           color = genotype)) + 
  geom_point() + theme_classic(base_size = 10) + 
  labs(x = "Sample Number", y = "CD1 SNP Number")

# SNP number still strongly correlates with percentage of CD1 supplement mass
ggplot(adj_snp_number_all, aes(x = percentCD1Mass, y = adj_snp_number)) +
  geom_point(color = "black") + theme_classic(base_size = 10) +
  geom_smooth(method='lm', formula = 'y ~ x', color = "black") +
  stat_regline_equation(label.y = 20000, color = "black", 
                        aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 19000, color = "black",
                        aes(label = ..rr.label..)) +
  labs(x = "CD1 Mass Percentage", y = "CD1 SNP Number")

# CD1 SNPs have varying penetrance, partially due to different mass percentages
ggplot(adj_samples, aes(x = percentAC, color = as.factor(percentCD1Mass))) +
  stat_ecdf(size = 1) + theme_classic(base_size = 10) +
  labs(x = "CD1 SNP Penetrance")

# CD1 SNPs have varying penetrance, partially due to different mass percentages
ggplot(adj_samples, aes(x = 100*percentAC, 
                         color = as.factor(round(100*percentCD1Mass)))) +
  geom_density(size = 1) + theme_classic(base_size = 10) + 
  guides(color=guide_legend(title="% CD1 Mass")) +
  labs(x = "CD1 SNP Penetrance")


## ----snpsTOgenes, include = F------------------------------------------------------------------------------------------------------------
# loading genes and info for mouse
mouse_genes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
mouse_ensembl <- as.data.frame(org.Mm.egENSEMBL)
mouse_chr <- as.data.frame(org.Mm.egCHR)
mouse_pseudo <- as.data.frame(org.Mm.egGENETYPE)

# mouse gene dataframe and grange
mouse_symbols <- as.data.frame(org.Mm.egSYMBOL) %>%
  full_join(mouse_ensembl, by = "gene_id") %>%
  full_join(mouse_pseudo, by = "gene_id")
grange_samples <- adj_samples %>%
  select(-CHROM) %>%
  separate(location, c("chrom", "start"), convert = T) %>%
  mutate(end = start, chrom = paste0('chr', chrom)) %>%
  ungroup() %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) 

# converting variants from position info to gene info
overlap_samples <- mergeByOverlaps(grange_samples, mouse_genes) %>%
  as.list() %>%
  as.data.frame() %>%
  select(-contains("mouse"), -contains("strand"),
         -contains("end"), -contains("width"), -contains("start")) %>%
  select(contains("grange"), contains("gene_id")) %>%
  rename(grange_samples.gene_id = gene_id)
names(overlap_samples) <- substring(names(overlap_samples), 16)

gene_samples <- overlap_samples %>%
  group_by(sample, gene_id) %>%
  summarize(CHR = unique(seqnames), gene_RC = sum(RC),
            gene_AC = sum(AC), gene_TC = sum(TC), 
            replicate = unique(replicate),genotype = unique(genotype),
            percentCD1Mass = unique(percentCD1Mass),
            percentLabeledMass = unique(percentLabeledMass)) %>%
  mutate(geneAdjPercentAC = gene_AC/gene_TC) %>%
  left_join(mouse_ensembl, by = "gene_id") %>%
  group_by(sample)


## ----vizGenes----------------------------------------------------------------------------------------------------------------------------
# gathering information for SNP-containing genes
gene_number <- gene_samples %>%
  group_by(sample) %>%
  summarize(gene_number = n(),
            genotype = unique(genotype),
            percentCD1Mass = unique(percentCD1Mass),
            percentLabeledMass = unique(percentLabeledMass))

# including zeroes for non-snp samples
gene_number_all <- data.frame(sample = c(14, 15, 16, 17, 18),
                                 gene_number = c(0, 0, 0, 0, 0),
                                 genotype = c("wt", "het", "het", "het", "het"),
                                 percentCD1Mass = c(0, 0, 0, 0, 0),
                                 percentLabeledMass = c(1, 1, 1, 1, 1)) %>%
  bind_rows(gene_number)

# number of SNP-containing genes correlates with mass of CD1 supplement
ggplot(gene_number_all, 
       aes(x = sample, y = gene_number, color = genotype)) + 
  geom_point() + theme_classic(base_size = 10) + 
  labs(y = "CD1 Gene Number", x = "Sample Number")
snp_mass <- ggplot(gene_number_all, aes(x = percentCD1Mass * 100,
                                        y = gene_number)) + 
  geom_point(color = "black", size = 0.25) + 
  theme_classic(base_size = 7) +
  geom_smooth(method='lm', formula = 'y ~ x', color = "black",
              linewidth = 0.5) +
  stat_regline_equation(label.y = 3700, color = "black", 
                        aes(label = ..eq.label..),
                        size = 2) +
  stat_regline_equation(label.y = 3400, color = "black",
                        aes(label = ..rr.label..),
                        size = 2) +
  theme(axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          axis.ticks = element_line(color = "black")) +
  labs(y = "# of Genes with CD1 SNPs", 
       x = "Relative CD1 Mass (%)"); snp_mass
ggsave("supp_fig_8h.pdf", snp_mass, units = "in",
       height = 1.72, width = 1.72)

# CD1 genes have varying penetrance, partially due to different mass percentages
snp_density <- ggplot(gene_samples, aes(x = 100*geneAdjPercentAC, 
                         color = as.factor(round(100*percentCD1Mass)))) +
  geom_density(size = 1) + theme_classic(base_size = 10) + 
  guides(color=guide_legend(title="% CD1 Mass")) +
  scale_color_manual(values = c(cbPalette[2],
                                cbPalette[4],
                                cbPalette[3],
                                cbPalette[5])) +
  theme(axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          axis.ticks = element_line(color = "black"),
        legend.position = c(0.8, 0.8)) +
  labs(x = "CD1 SNP Penetrance (%)", y = "Density"); snp_density

# creating and using a function for rowwise binomial testing of variant penetrance
binom.p <- function(x, n, p){binom.test(x, n, p, alternative="less")$p.value}
gene_samples$binom.pValue <- mapply(binom.p,
                                    gene_samples$gene_AC,
                                    gene_samples$gene_TC,
                                    gene_samples$percentCD1Mass)
gene_samples$p.adjust <- p.adjust(gene_samples$binom.pValue, method = "fdr")
gene_samples_final <- gene_samples %>%
  mutate(ambience = case_when(binom.pValue < 0.1 ~ "real",
                              TRUE ~ "ambient"),
         adj_ambience = case_when(p.adjust < 0.1 ~ "real",
                                  TRUE ~ "ambient"),
         cd1_snps_over_expected = geneAdjPercentAC / percentCD1Mass)

# making supplemental table 9
supp_table_9 <- gene_samples_final %>%
  ungroup() %>%
  distinct(gene_id, genotype, replicate, cd1_snps_over_expected,
           binom.pValue, p.adjust, ensembl_id) %>%
  rename(entrez_id = gene_id,
         cd1_snps_over_expected_ratio = cd1_snps_over_expected,
         cd1_snps_over_expected_pval = binom.pValue,
         cd1_snps_over_expected_qval = p.adjust,
         ensembl_gene_id = ensembl_id)
write.csv(supp_table_9, "supp_table_9.csv", quote = F, row.names = F)

# exclusion list for genes with "bad" snps in > 2/6 samples
ensembl_percents <- gene_samples %>%
  pivot_wider(names_from = sample, values_from = p.adjust,
              names_prefix = "number_") %>%
  group_by(ensembl_id) %>%
  summarize(number_12 = mean(number_12, na.rm = T),
            number_13 = mean(number_13, na.rm = T),
            number_19 = mean(number_19, na.rm = T),
            number_20 = mean(number_20, na.rm = T),
            number_21 = mean(number_21, na.rm = T),
            number_22 = mean(number_22, na.rm = T))
ensembl_percents[ensembl_percents < 0.1] <- NaN
ensembl_percents$number_good <- rowSums(is.na(ensembl_percents))
ensembl_percents$number_bad <- 6 - ensembl_percents$number_good
ensembl_percents <- filter(ensembl_percents, number_bad > 2)


## ----goAnalysis--------------------------------------------------------------------------------------------------------------------------
# loading bcl11a gc transcriptome data
gc_universe <- read_excel("external_tables.xlsx", sheet = "durak_dags") %>%
  distinct(ensembl_id)

# identifying relavent genes
gc_ambient <- gc_universe %>%
  filter(ensembl_id %in% ensembl_percents$ensembl_id)
gc_real <- gc_universe %>%
  filter(!(ensembl_id %in% ensembl_percents$ensembl_id))
gc_snp_real <- gc_universe %>%
  filter(ensembl_id %in% gene_samples$ensembl_id) %>%
  filter(!(ensembl_id %in% ensembl_percents$ensembl_id))

# performing GO analysis
go_ambient <- enrichGO(gc_ambient$ensembl_id,
                             OrgDb = org.Mm.eg.db,
                             keyType = "ENSEMBL",
                             ont = "ALL",
                             pAdjustMethod = "BH",
                             readable = T,
                             pool = F,
                       universe = gc_universe$ensembl_id) %>%
  as.data.frame() %>%
  mutate(comp = "ambient")
go_snp_real <- enrichGO(gc_snp_real$ensembl_id,
                             OrgDb = org.Mm.eg.db,
                             keyType = "ENSEMBL",
                             ont = "ALL",
                             pAdjustMethod = "BH",
                             readable = T,
                             pool = F,
                        universe = gc_universe$ensembl_id) %>%
  as.data.frame() %>%
  mutate(comp = "snp_real")

go_df_realVSambient <- bind_rows(go_ambient, go_snp_real) %>%
  separate(GeneRatio, c("Num", "Denom")) %>%
  mutate(GeneRatio = as.numeric(Num) / as.numeric(Denom))

# reducing GO terms
go_df_realVSambient_comps <- go_df_realVSambient %>%
  distinct(comp)
go_df_realVSambient_reduced <- data.frame()
for (i in 1:nrow(go_df_realVSambient_comps)) {

  plot_comp <- go_df_realVSambient_comps[i, ]

  plot_go_df <- go_df_realVSambient %>%
    filter(comp == plot_comp) %>%
    arrange(qvalue) %>%
    mutate(trans_qvalue = -log10(qvalue))

  go_df_reduceTest <- plot_go_df %>%
    group_by(ONTOLOGY) %>%
    summarize(ont_count = n()) %>%
    filter(ont_count == 1)

  rows_to_bind <- plot_go_df %>%
    filter(ONTOLOGY %in% go_df_reduceTest$ONTOLOGY)

  plot_go_df <- plot_go_df %>%
    filter(!(ONTOLOGY %in% go_df_reduceTest$ONTOLOGY))

  go_df_preReduce <- plot_go_df %>%
    select(ONTOLOGY, ID) %>%
    rename(go_type = ONTOLOGY,
           go_id = ID)

  go_df_scores <- plot_go_df$trans_qvalue %>%
    setNames(plot_go_df$ID)

  go_df_reduce <- go_reduce(go_df_preReduce,
                            orgdb = "org.Mm.eg.db",
                            threshold = 0.7,
                            scores = go_df_scores,
                            measure = "Wang")

  plot_go_df_reduce <- plot_go_df %>%
    filter(ID %in% go_df_reduce$parent_id) %>%
    bind_rows(rows_to_bind)

  go_df_realVSambient_reduced <- bind_rows(go_df_realVSambient_reduced, plot_go_df_reduce)
}

# plotting reduced GO terms
go_df_realVSambient_comps <- go_df_realVSambient_reduced %>%
  distinct(comp)
for (i in 1:nrow(go_df_realVSambient_comps)) {
  
  plot_comp <- go_df_realVSambient_comps[i, ]
  
  plot_go_df <- go_df_realVSambient_reduced %>%
    filter(comp == plot_comp) %>%
    arrange(qvalue) %>%
    group_by(ONTOLOGY) %>%
    slice(1:10) %>%
    ungroup() %>%
    mutate(Description = paste0(ONTOLOGY,
                                ": ",
                                Description),
           Description = fct_reorder2(Description,
                                         -qvalue,
                                         qvalue))
  
  print(ggplot(plot_go_df, aes(x = GeneRatio, y = Description, 
                         fill = qvalue, size = Count)) +
    geom_point(color = "black", pch = 21) +
    facet_wrap(~ ONTOLOGY) + 
    labs(y = "", title = paste0("GO Enrichment: ", plot_comp)) +
    theme_classic(base_size = 10, base_family = "sans") +
    scale_fill_distiller(palette = "RdBu", direction = 1) +
    theme(axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          axis.ticks = element_line(color = "black")))
}

# final go plot
bcl11a_snp_go <- ggplot(go_df_realVSambient_reduced %>%
                          #filter(ONTOLOGY == "MF") %>%
                          arrange(qvalue) %>%
                          slice(1:15) %>%
                          mutate(Description = case_when(Description == 
                                                           "regulation of postsynaptic membrane neurotransmitter receptor levels" ~ 
                                                           "regulation of postsynaptic receptor levels",
                                                         Description == 
                                                           "phosphotransferase activity, alcohol group as acceptor" ~
                                                           "phosphotransferase activity",
                                                         T ~ Description)) %>%
                          mutate(Description = fct_reorder2(Description,
                                                            -qvalue, 
                                                            qvalue)), 
                        aes(x = GeneRatio, y = Description, 
                            fill = qvalue, size = Count)) +
  geom_point(color = "black", pch = 21) +
  scale_size(range = c(0.25, 2.5)) +
  labs(y = "", title = "GO Analysis: Non-Ambient GC RNAs", size = "Number") +
  theme_classic(base_size = 7, base_family = "sans") +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); bcl11a_snp_go
ggsave("supp_fig_8i.pdf", bcl11a_snp_go, units = "in",
       height = 1.75, width = 3.5)


## ----spinWB------------------------------------------------------------------------------------------------------------------------------
# load western blot validation data (equal mass per sample by bradford assay)
mass <- read_excel("input_tables.xlsx", sheet = "Input Table 4") %>%
  filter(Marker == "Gap43" | Marker == "Map2")

# calculating significance values for wGCF signal compared to gcf signal
mass_sig <- mass %>%
  filter(Condition != "PNH") %>%
  select(-contains("Percent")) %>%
  pivot_wider(names_from = Condition, values_from = Signal) %>%
  group_by(Marker) %>%
  summarize(pval = t.test(wGCF, GCF, alternative = "less")$p.value) %>%
  mutate(p_symbol = case_when(pval < 0.0001 ~ "****", 
                              pval < 0.001 ~ "***",
                              pval < 0.01 ~ "**",
                              pval < 0.05 ~ "*",
                              TRUE ~ "ns"))

# plot western blot validation data
wb_validate <- ggplot(mass %>%
                        
                        # adding significance values
                        left_join(mass_sig, by = "Marker") %>%
                        
                        # plotting "washed" GCFs                
                        filter(Condition == "wGCF"), 
                      aes(x = as.factor(Marker), y = as.numeric(percentOfGCF),
                          fill = Condition)) +
  geom_bar(stat = "summary", position = "dodge") +
  geom_errorbar(stat='summary', width = 0.25, 
                fun.data = mean_cl_boot, linewidth = 0.5) +
  geom_point(size = 0.5) +
  geom_text(aes(x = Marker, y = 175, label = p_symbol), size = 2) +
  theme_classic(base_size = 7) +
  labs(x = "", y = "wGCF WB Signal (% of GCF)", fill = "") +
  scale_fill_manual(values = c("#CA6368")) +
  geom_hline(yintercept = 100, linetype = "dashed", 
             color = "black", linewidth = 0.25) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none"); wb_validate
ggsave("supp_fig_8l.pdf", wb_validate, units = "in", width = 2, height = 1.55)


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


## ----loadingWashSNPs---------------------------------------------------------------------------------------------------------------------
# loading gc wash snp data
allele_gcs <- read_excel("input_tables.xlsx", sheet = "Input Table 14",
                         col_types = c("text", "numeric", "text", "text", "text",
                                       "text", "numeric", "text", "numeric"))

# loading soma wash snp data
allele_somata <- read_excel("input_tables.xlsx", sheet = "Input Table 15",
                         col_types = c("text", "numeric", "text", "text", "text",
                                       "text", "numeric", "text", "numeric"))

# combining gc and soma wash snp data
allele_all <- bind_rows(allele_gcs, allele_somata) %>%
  
  separate("counts", c("ref", "alt"), sep = ",", 
           convert = T, fill = "right") %>%
  mutate(alt = case_when(is.na(alt) == T ~ 0, T ~ alt)) %>%
  mutate(bi_depth = ref + alt) %>%
  filter(FILTER == "PASS") %>%
  unite("position", CHROM:POS, sep = ".", remove = F) %>%
  mutate(position = gsub("MT", "M", position))

# focusing on relavent chromosomes
allele_all_chrom <- allele_all %>%
  distinct(CHROM) %>%
  mutate(CHROM = gsub("MT", "M", CHROM))
allele_all_chrom_filter <- allele_all_chrom %>%
  slice(1:22)
allele_all_chrom_keep <- allele_all %>%
  filter(CHROM %in% allele_all_chrom_filter$CHROM)


## ----allSNPsWash-------------------------------------------------------------------------------------------------------------------------
# distribution of total read counts for SNPs
total_snp_dist_spin <- ggplot(allele_all, aes(x = bi_depth)) + 
  geom_density(color = "black", linewidth = 1) + 
  theme_classic(base_size = 10) + scale_x_continuous(trans = "log10") + 
  labs(x = "Bialleleic Depth", y = "Density") + 
  geom_vline(xintercept = 9.66, color = "black", linetype = "dashed"); total_snp_dist_spin
total_snp_dist_values_spin <- ggplot_build(total_snp_dist_spin)
# write.csv(total_snp_dist_values_spin$data[[1]], 
#           file = "total_snp_dist_values_spin.csv")

# setting a read filter of > 9 biallelic counts
high_alleles <- allele_all %>%  
  filter(bi_depth > 9)

# desnity plots of non-reference reads per sample
ggplot(high_alleles %>% filter(condition != "cd1_somata",
                               condition != "b6_somata"), 
       aes(x = alt/bi_depth*100, color = condition)) +
  geom_density(size = 1.5) + theme_classic(base_size = 10) +
  xlab("Percent of Non-Reference Reads per SNP") +
  ylab("Density")


## ----idCD1Wash---------------------------------------------------------------------------------------------------------------------------
snp_matching <- high_alleles %>%
  group_by(position, condition)%>%
  summarize(percent = mean(alt/bi_depth*100, na.rm = T)) %>%
  pivot_wider(names_from = condition, values_from = percent)

# removing SNPs that aren't in either soma sample
snp_simple_first <- snp_matching[rowSums(is.na(snp_matching[c("cd1_somata", "b6_somata")])) != 2, ]

# removing SNPs that are in both soma samples
snp_simple_second <- snp_simple_first[rowSums(is.na(snp_simple_first[c("cd1_somata", "b6_somata")])) != 0, ]

# deciding to use 85.37997% as cutoff
soma_snp_dist <- ggplot(snp_simple_second,
       aes(x = cd1_somata)) + geom_density(color = "black", size = 1.5) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  geom_vline(xintercept = 85.37997, linetype = "dashed", linewidth = 0.5) +
  theme_classic(base_size = 12) +
  xlab("Percent of Non-Reference Reads per SNP") + 
  ylab("Density"); soma_snp_dist
soma_snp_dist_values <- ggplot_build(soma_snp_dist)
# write.csv(soma_snp_dist_values$data[[1]], 
#           file = "soma_snp_dist_values.csv")

# filter on CD1 SNPs
snp_spectrum_cd1 <- snp_simple_second %>%
  drop_na(cd1_somata) %>%
  filter(cd1_somata > 85.37997) %>%
  pivot_longer(contains("_") | contains("supers"), 
               names_to = "condition", values_to = "snp_spectrum_cd1") %>%
  drop_na(snp_spectrum_cd1)

# remove duplicate SNP positions
cd1_snps <- unique(snp_spectrum_cd1$position)

# filter for cd1 snps
allConditions_onlyCD1 <- high_alleles %>%
  filter(position %in% cd1_snps)

# density plots for CD1 SNPs
ggplot(allConditions_onlyCD1 %>% filter(condition != "cd1_somata"), 
       aes(x = alt/bi_depth*100, color = condition)) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  geom_density(size = 1.5) + theme_classic(base_size = 10) +
  xlab("Percent of Non-Reference Reads per CD1 SNP") +
  ylab("Density")


## ----cd1GenesWash------------------------------------------------------------------------------------------------------------------------
# loading genes and info for mouse
mouse_genes <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
mouse_ensembl <- as.data.frame(org.Mm.egENSEMBL)
mouse_chr <- as.data.frame(org.Mm.egCHR)
mouse_pseudo <- as.data.frame(org.Mm.egGENETYPE)

# mouse gene dataframe and grange
mouse_symbols <- as.data.frame(org.Mm.egSYMBOL) %>%
  full_join(mouse_ensembl, by = "gene_id") %>%
  full_join(mouse_pseudo, by = "gene_id")
cd1_grange <- allConditions_onlyCD1 %>%
  select(-CHROM) %>%
  separate(position, c("chrom", "start"), convert = T) %>%
  mutate(end = start, chrom = paste0('chr', chrom)) %>%
  ungroup() %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) 

# converting variants from position info to gene info
cd1_overlap <- mergeByOverlaps(cd1_grange, mouse_genes) %>%
  as.list() %>%
  as.data.frame() %>%
  select(-contains("mouse"), -contains("strand"),
         -contains("end"), -contains("width"), -contains("start")) %>%
  select(contains("grange"), contains("gene_id")) %>%
  rename(cd1_grange.gene_id = gene_id)
names(cd1_overlap) <- substring(names(cd1_overlap), 12)

allConditions_snpGenes <- cd1_overlap %>%
  group_by(condition, replicate, gene_id) %>%
  summarize(CHR = unique(seqnames), gene_ref = sum(ref), gene_alt = sum(alt),
            gene_bi_depth = sum(bi_depth)) %>%
  mutate(alt_per = gene_alt / gene_bi_depth * 100) %>%
  left_join(mouse_ensembl, by = "gene_id") %>%
  ungroup()

# density plots for combined CD1 SNP Genes
cd1_snp_dist <- ggplot(allConditions_snpGenes %>% filter(condition != "cd1_somata"), 
       aes(x = alt_per)) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  geom_density(size = 1.5) + theme_classic(base_size = 10) +
  geom_vline(xintercept = 90.16302, color = "black", 
             linewidth = 0.5, linetype = "dashed") +
  xlab("Percent of Non-Reference Reads per CD1 Gene") +
  ylab("Density"); cd1_snp_dist
cd1_snp_dist_values <- ggplot_build(cd1_snp_dist)
# write.csv(cd1_snp_dist_values$data[[1]], 
#           file = "cd1_snp_dist_values.csv")

# making supplemental table 10
supp_table_10 <- allConditions_snpGenes %>%
  distinct(gene_id, condition, replicate, alt_per, ensembl_id) %>%
  rename(entrez_id = gene_id,
         cd1_snp_penetrance = alt_per,
         ensembl_gene_id = ensembl_id)
write.csv(supp_table_10, "supp_table_10.csv", quote = F, row.names = F)

# density plots for CD1 SNP Genes
spin_snp_dens <- ggplot(allConditions_snpGenes %>% 
                          filter(condition != "cd1_somata") %>%
                          mutate(condition = case_when(condition == "spin_GCs" ~ 
                                                         "wGCs",
                                                       condition == "noSpin_GCs" ~ 
                                                         "GCs",
                                                       condition == "supers" ~ 
                                                         "Supers")) %>%
                          mutate(condition = as.factor(condition)) %>%
                          mutate(condition = fct_relevel(condition, 
                                                        "wGCs", 
                                                        "GCs", 
                                                        "Supers")), 
                        aes(x = alt_per, color = condition)) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  geom_vline(xintercept = 90.16302, color = "black", 
             linewidth = 0.25, linetype = "dashed") +
  geom_vline(xintercept = 50, color = "black", 
             linewidth = 0.25, linetype = "dashed") +
  geom_density(linewidth = 0.5) + theme_classic(base_size = 7) +
  scale_color_manual(values = c("#CA6368", "#653234", "#d1d3d4")) +
  labs(x = "CD1 SNP Penetrance (%)", y = "Density", color = "") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = c(0.1, .9)); spin_snp_dens
ggsave("supp_fig_9c.pdf", spin_snp_dens, units = "in",
       height = 2.16, width = 2.02)

# summarizing confidence by condition
cd1_snp_confidence <- allConditions_snpGenes %>%
  mutate(high_confidence = case_when(alt_per > 90.16302 ~ T,
                                     T ~ F),
         low_confidence = case_when(alt_per > 50 ~ T,
                                    T ~ F),
         total = T) %>%
  group_by(condition) %>%
  summarize(high_sum = sum(high_confidence),
            low_sum = sum(low_confidence),
            total_sum = sum(total)) %>%
  mutate(high_percent = high_sum / total_sum * 100,
         low_percent = low_sum / total_sum * 100)
print(cd1_snp_confidence)


## ----deseq2run---------------------------------------------------------------------------------------------------------------------------
# loading featureCounts
spin_counts <- read_excel("input_tables.xlsx", sheet = "Input Table 16") %>%
  pivot_longer(cols = contains("rep"), 
               names_to = "sample", values_to = "counts") %>%
  select(Geneid, sample, counts) %>%
  mutate(sample = gsub("dt_", "", sample)) %>%
  mutate(sample = gsub("_gc", "", sample)) %>%
  separate("sample", into = c("condition", "replicate"), remove = F)

# making metadata frame
spin_metadata <- data.frame(sample = c("noSpin_rep1",
                                       "noSpin_rep2",
                                       "noSpin_rep3",
                                       "spin_rep1",
                                       "spin_rep2",
                                       "spin_rep3",
                                       "super_rep1",
                                       "super_rep2",
                                       "super_rep3"),
                            condition = c("noSpin", "noSpin", "noSpin",
                                          "spin", "spin", "spin",
                                          "super", "super", "super"),
                            replicate = c(1, 2, 3, 1, 2, 3, 1, 2, 3))

# making and filtering count matrix
all_na_counts_spin <- spin_counts %>%
  select(Geneid, sample, counts) %>%
  pivot_wider(names_from = "sample", values_from = "counts")
all_na_counts_spin[all_na_counts_spin == 0] <- NA
all_filtered_counts_spin <- all_na_counts_spin %>%
  filter(!if_all(.cols = noSpin_rep1:noSpin_rep3, .fns = is.na) |
           !if_all(.cols = spin_rep1:spin_rep3, .fns = is.na) |
           !if_all(.cols = super_rep1:super_rep3, .fns = is.na)) %>%
  column_to_rownames(var = "Geneid")
all_filtered_counts_spin[is.na(all_filtered_counts_spin)] <- 0

# checking row and column names
rownames(spin_metadata) <- colnames(all_filtered_counts_spin)
all(rownames(spin_metadata) %in% colnames(all_filtered_counts_spin))
all(rownames(spin_metadata) == colnames(all_filtered_counts_spin))

# running DESeq2
spin_dds <- DESeqDataSetFromMatrix(countData = all_filtered_counts_spin,
                                  colData = spin_metadata,
                                  design = ~ condition)
spin_dds <- DESeq(spin_dds)

# loose filter to remove junk
spin_dds <- spin_dds[rowSums(counts(spin_dds)) >= 10,]

# transforming values with rlog
spin_rld_blind <- rlog(spin_dds, blind = T)

# pca for top 500 genes
spin_pca_var <- apply(assay(spin_rld_blind), 1, sd)
spin_pca_var_df <- assay(spin_rld_blind)[order(spin_pca_var, 
                       decreasing = TRUE)[seq_len(500)],]
spin_pca <- prcomp(t(spin_pca_var_df), scale = FALSE)
spin_pca_df <- spin_pca$x %>% data.frame() %>% rownames_to_column() %>% 
  left_join(., data.frame(colData(spin_rld_blind)), 
            by = c(rowname = "sample"))
spin_pca_percent <- round(100 * spin_pca$sdev^2/sum(spin_pca$sdev^2), 1)

# plotting PC1 and PC2 for top 500 genes
spin_pca_plot <- ggplot(spin_pca_df %>%
                     mutate(condition = case_when(condition == "spin" ~ 
                                                         "wGCs",
                                                       condition == "noSpin" ~ 
                                                         "GCs",
                                                       condition == "super" ~ 
                                                         "Supers")) %>%
                          mutate(condition = as.factor(condition)) %>%
                          mutate(condition = fct_relevel(condition, 
                                                        "wGCs", 
                                                        "GCs", 
                                                        "Supers")),
                   aes(get(paste0("PC", 1)), 
                       get(paste0("PC", 2)),
                       col = condition)) + 
  scale_color_manual(values = c("#CA6368", "#653234", "#d1d3d4")) +
  labs(title = "GC Wash: PC Analysis", 
       x = paste0("PC", 1, ": ", spin_pca_percent[1], "%"), 
       y = paste0("PC", 2, ": ", spin_pca_percent[2], "%"),
       color = "") + 
  coord_fixed() +
  theme_classic(base_size = 7) +
  geom_point(size = 0.25) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = c(0.2, .9)); spin_pca_plot
ggsave("supp_fig_9d.pdf", spin_pca_plot, units = "in",
       height = 1.39, width = 2.88)

# making deseq2 dataframes
dds_spin_super <- spin_dds
dds_spin_super$condition <- relevel(dds_spin_super$condition, ref = "super")
deseq_spin_super <- DESeq(dds_spin_super)
resultsNames(deseq_spin_super)
res_spin_super <- results(deseq_spin_super, name = "condition_spin_vs_super")
summary(res_spin_super)
plotMA(res_spin_super)
res_shrink_spin_super <- lfcShrink(deseq_spin_super, coef = "condition_spin_vs_super")
summary(res_shrink_spin_super)
plotMA(res_shrink_spin_super)
deseq_df_spin_super <- as.data.frame(res_shrink_spin_super) %>%
  mutate(comp = "spin_over_super") %>%
  rownames_to_column("Geneid")

dds_spin_noSpin <- spin_dds
dds_spin_noSpin$condition <- relevel(dds_spin_noSpin$condition, ref = "noSpin")
deseq_spin_noSpin <- DESeq(dds_spin_noSpin)
resultsNames(deseq_spin_noSpin)
res_spin_noSpin <- results(deseq_spin_noSpin, name = "condition_spin_vs_noSpin")
summary(res_spin_noSpin)
plotMA(res_spin_noSpin)
res_shrink_spin_noSpin <- lfcShrink(deseq_spin_noSpin, coef = "condition_spin_vs_noSpin")
summary(res_shrink_spin_noSpin)
plotMA(res_shrink_spin_noSpin)
deseq_df_spin_noSpin <- as.data.frame(res_shrink_spin_noSpin) %>%
  mutate(comp = "spin_over_noSpin") %>%
  rownames_to_column("Geneid")

dds_noSpin_super <- spin_dds
dds_noSpin_super$condition <- relevel(dds_noSpin_super$condition, ref = "super")
deseq_noSpin_super <- DESeq(dds_noSpin_super)
resultsNames(deseq_noSpin_super)
res_noSpin_super <- results(deseq_noSpin_super, name = "condition_noSpin_vs_super")
summary(res_noSpin_super)
plotMA(res_noSpin_super)
res_shrink_noSpin_super <- lfcShrink(deseq_noSpin_super, coef = "condition_noSpin_vs_super")
summary(res_shrink_noSpin_super)
plotMA(res_shrink_noSpin_super)
deseq_df_noSpin_super <- as.data.frame(res_shrink_noSpin_super) %>%
  mutate(comp = "noSpin_over_super") %>%
  rownames_to_column("Geneid")

# saving deseq dataframe
spin_deseq_df <- bind_rows(deseq_df_spin_super,
                           deseq_df_spin_noSpin,
                           deseq_df_noSpin_super) %>%
  mutate(sig = case_when(padj < 0.01 ~ T, T ~ F),
         direction = case_when(log2FoldChange > 0 ~ "up",
                               log2FoldChange < 0 ~ "down")) %>%
  unite("full_comp", "comp", "direction", remove = F) %>%
  filter(is.na(padj) == F)

# making supplemental table 11
supp_table_11 <- spin_deseq_df %>%
  distinct(Geneid, log2FoldChange, pvalue, padj, comp) %>%
  rename(ensembl_gene_id = Geneid,
         log2FC = log2FoldChange,
         pval = pvalue,
         qval = padj) %>%
  pivot_wider(names_from = comp, values_from = c("log2FC", "pval", "qval")) %>%
  relocate(contains("spin_over_noSpin"), .after = "ensembl_gene_id") %>%
  relocate(contains("spin_over_super", ignore.case = F), 
           .after = "qval_spin_over_noSpin")
write.csv(supp_table_11, "supp_table_11.csv", quote = F, row.names = F)


## ----deseq2analyze-----------------------------------------------------------------------------------------------------------------------
# ma plots
spin_ma_plot <- ggplot(spin_deseq_df %>%
                         mutate(comp = case_when(comp == "spin_over_super" ~ 
                                                   "wGCs_over_Supers",
                                                 comp == "spin_over_noSpin" ~
                                                   "wGCs_over_GCs",
                                                 comp == "noSpin_over_super" ~
                                                   "GCs_over_Supers")) %>%
                         mutate(comp = as.factor(comp),
                                sig = as.factor(sig)) %>%
                         mutate(comp = fct_relevel(comp, "wGCs_over_GCs",
                                                   "wGCs_over_Supers",
                                                   "GCs_over_Supers"),
                                sig = fct_relevel(sig, "TRUE", "FALSE")),
                       aes(x = log10(baseMean), y = log2FoldChange,
                          color = sig)) +
  labs(x = "log10(Mean of Normalized Counts)",
       y = "log2(Fold Change)",
       color = "q < 0.01") +
  scale_color_manual(values = c("black", "grey50")) +
  geom_point(size = 0.000000000000000000000000000000000000000000000000000001) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = c(0.9, 0.85)) +
  facet_wrap(~ comp); spin_ma_plot

# volcano plots
spin_volcano <- ggplot(spin_deseq_df %>%
                         mutate(comp = case_when(comp == "spin_over_super" ~ 
                                                   "wGCs_over_Supers",
                                                 comp == "spin_over_noSpin" ~
                                                   "wGCs_over_GCs",
                                                 comp == "noSpin_over_super" ~
                                                   "GCs_over_Supers")) %>%
                         mutate(comp = as.factor(comp),
                                sig = as.factor(sig)) %>%
                         mutate(comp = fct_relevel(comp, "wGCs_over_GCs",
                                                   "wGCs_over_Supers",
                                                   "GCs_over_Supers"),
                                sig = fct_relevel(sig, "TRUE", "FALSE")),
       aes(x = log2FoldChange, y = -log10(pvalue),
                          color = sig)) +
  labs(x = "log2(Fold Change)",
       y = "-log10(p-value)",
       color = "q < 0.01") +
  scale_color_manual(values = c("black", "grey50")) +
  geom_point(size = 0.000000000000000000000000000000000000000000000000000001) +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = c(0.9, 0.85)) +
  facet_wrap(~ comp); spin_volcano
ggsave("supp_fig_9e.pdf", spin_volcano, units = "in",
       height = 2.16, width = 4.08)

# upset plot
spin_deseq_upset_dir <- spin_deseq_df %>%
  mutate(comp = case_when(comp == "spin_over_super" ~ 
                            "wGCs_over_Supers",
                          comp == "spin_over_noSpin" ~
                            "wGCs_over_GCs",
                          comp == "noSpin_over_super" ~
                            "GCs_over_Supers")) %>%
  mutate(upset_value = case_when(sig == T ~ 1,
                                 TRUE ~ 0),
         direction = case_when(log2FoldChange > 0 ~ "up",
                               log2FoldChange < 0 ~ "down")) %>%
  unite(full_comp, comp, direction) %>%
  select(Geneid, full_comp, upset_value) %>%
  pivot_wider(names_from = "full_comp", values_from = "upset_value") %>%
  as.data.frame()
spin_deseq_upset_dir[is.na(spin_deseq_upset_dir)] <- 0
plot_spin_deseq_upset_dir <- upset(spin_deseq_upset_dir,
                                   nsets = 6,
                                   order.by = "freq",
                                   sets = c("GCs_over_Supers_down", 
                                            "GCs_over_Supers_up",
                                            "wGCs_over_Supers_down",
                                            "wGCs_over_Supers_up",
                                            "wGCs_over_GCs_down",
                                            "wGCs_over_GCs_up"),
                                   keep.order = T,
                                   
                                   matrix.color = "black", main.bar.color = "black",
                                   sets.bar.color = "black", shade.color = "gray75",
                                   text.scale = 1.25); plot_spin_deseq_upset_dir

# high-confidence genes in deseq2
deseq_snpConfidence <- spin_deseq_df %>%
  left_join(allConditions_snpGenes %>%
              filter(condition != "cd1_somata") %>%
              rename(Geneid = ensembl_id) %>%
              mutate(high_confidence = case_when(alt_per > 90.16302 ~ T,
                                                 T ~ F),
                     low_confidence = case_when(alt_per > 50 ~ T,
                                                T ~ F)) %>%
              group_by(Geneid) %>%
              summarize(high_confident = max(high_confidence),
                        low_confident = max(low_confidence)),
            by = "Geneid") %>%
  filter(sig == T) %>%
  group_by(full_comp) %>%
  summarize(num_high = sum(high_confident, na.rm = T),
            num_low = sum(low_confident, na.rm = T),
            total = n(),
            percent_high = num_high / total * 100,
            percent_low = num_low / total * 100)
print(deseq_snpConfidence)

# snp confidence
snp_confidence <- allConditions_snpGenes %>%
  filter(condition != "cd1_somata", alt_per > 90.16302) %>%
  distinct(ensembl_id) %>%
  drop_na()

# differential expression by condition
spin_deseq_sig <- spin_deseq_df %>%
  mutate(isSNP = case_when(Geneid %in% snp_confidence$ensembl_id ~ T,
         T ~ F)) %>%
  group_by(full_comp, isSNP) %>%
  count(sig) %>%
  ungroup()

# SNP odds ratios
spin_snpOdds <- data.frame()
for (i in 1:(nrow(spin_deseq_sig)/4)) {
  
  plot_start <- spin_deseq_sig %>%
    slice((4*i-3):(4*i)) %>%
    pivot_wider(names_from = "sig", values_from = "n") %>%
    arrange(-isSNP) %>%
    relocate("TRUE", .before = "FALSE")
  
  plot_odds <- fisher.test(plot_start %>% select(3, 4))
  
  plot_results <- data.frame(full_comp = plot_start[[1,1]],
                             odds_ratio = plot_odds$estimate[[1]],
                             lower_limit = plot_odds$conf.int[[1]],
                             upper_limit = plot_odds$conf.int[[2]],
                             p_value = plot_odds$p.value[[1]])
  
  spin_snpOdds <- bind_rows(spin_snpOdds, plot_results)
  
}

# plotting snp odds ratios
spin_snpOdds_plot <- spin_snpOdds %>%
  mutate(Significant = case_when(p_value < 0.05 ~ T, T ~ F)) %>%
  mutate(Significant = as.factor(Significant),
         full_comp = case_when(full_comp == "spin_over_super_up" ~ 
                                 "wGCs_over_Supers",
                               full_comp == "spin_over_super_down" ~
                                 "Supers_over_wGCs",
                               full_comp == "spin_over_noSpin_up" ~ 
                                 "wGCs_over_GCs",
                               full_comp == "spin_over_noSpin_down" ~
                                 "GCs_over_wGCs",
                               full_comp == "noSpin_over_super_up" ~ 
                                 "GCs_over_Supers",
                               full_comp == "noSpin_over_super_down" ~
                                 "Supers_over_GCs")) %>%
  mutate(Significant = fct_relevel(Significant, "TRUE", "FALSE"),
         full_comp = fct_relevel(full_comp, "wGCs_over_GCs", 
                                 "GCs_over_wGCs",
                                 "wGCs_over_Supers",
                                 "Supers_over_wGCs",
                                 "GCs_over_Supers",
                                 "Supers_over_GCs"))

spin_odds <- ggplot(spin_snpOdds_plot, aes(x = odds_ratio, y = as.factor(full_comp), 
                                         color = Significant)) +
  geom_point(size = 0.5) +
  geom_errorbarh(aes(xmin = lower_limit, xmax = upper_limit), height = 0.25) +
  geom_vline(xintercept = 1, color = "black", linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c("black", "gray50")) +
  labs(x = "Odds Ratio of SNP DAGs",
       y = "",
       color = "p < 0.05") +
  scale_x_continuous(trans = "log10") +
  scale_y_discrete(limits = rev) +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = c(0.9, 0.9)); spin_odds
ggsave("supp_fig_9f.pdf", spin_odds, units = "in",
       height = 0.89, width = 2.00)


## ----deseq2GO----------------------------------------------------------------------------------------------------------------------------
# loading spin universe
spin_universe <- spin_deseq_df %>%
  distinct(Geneid)

# identifying significant genes
spin_over_super <- spin_deseq_df %>%
  filter(comp == "spin_over_super", sig == T, log2FoldChange > 0)
super_over_spin <- spin_deseq_df %>%
  filter(comp == "spin_over_super", sig == T, log2FoldChange < 0)
spin_over_noSpin <- spin_deseq_df %>%
  filter(comp == "spin_over_noSpin", sig == T, log2FoldChange > 0)
noSpin_over_spin <- spin_deseq_df %>%
  filter(comp == "spin_over_noSpin", sig == T, log2FoldChange < 0)
noSpin_over_super <- spin_deseq_df %>%
  filter(comp == "noSpin_over_super", sig == T, log2FoldChange > 0)
super_over_noSpin <- spin_deseq_df %>%
  filter(comp == "noSpin_over_super", sig == T, log2FoldChange < 0)

# running go analysis
go_spin_over_super <- enrichGO(spin_over_super$Geneid,
                             OrgDb = org.Mm.eg.db,
                             keyType = "ENSEMBL",
                             ont = "ALL",
                             pAdjustMethod = "BH",
                             readable = T,
                             pool = F,
                             universe = spin_universe$Geneid) %>%
  as.data.frame() %>%
  mutate(comp = "spin_over_super", direction = "spin")
go_super_over_spin <- enrichGO(super_over_spin$Geneid,
                             OrgDb = org.Mm.eg.db,
                             keyType = "ENSEMBL",
                             ont = "ALL",
                             pAdjustMethod = "BH",
                             readable = T,
                             pool = F,
                             universe = spin_universe$Geneid) %>%
  as.data.frame() %>%
  mutate(comp = "super_over_spin", direction = "super")
go_spin_over_noSpin <- enrichGO(spin_over_noSpin$Geneid,
                             OrgDb = org.Mm.eg.db,
                             keyType = "ENSEMBL",
                             ont = "ALL",
                             pAdjustMethod = "BH",
                             readable = T,
                             pool = F,
                             universe = spin_universe$Geneid) %>%
  as.data.frame() %>%
  mutate(comp = "spin_over_noSpin", direction = "spin")
go_noSpin_over_spin <- enrichGO(noSpin_over_spin$Geneid,
                             OrgDb = org.Mm.eg.db,
                             keyType = "ENSEMBL",
                             ont = "ALL",
                             pAdjustMethod = "BH",
                             readable = T,
                             pool = F,
                             universe = spin_universe$Geneid) %>%
  as.data.frame() %>%
  mutate(comp = "noSpin_over_spin", direction = "noSpin")
go_noSpin_over_super <- enrichGO(noSpin_over_super$Geneid,
                             OrgDb = org.Mm.eg.db,
                             keyType = "ENSEMBL",
                             ont = "ALL",
                             pAdjustMethod = "BH",
                             readable = T,
                             pool = F,
                             universe = spin_universe$Geneid) %>%
  as.data.frame() %>%
  mutate(comp = "noSpin_over_super", direction = "noSpin")
go_super_over_noSpin <- enrichGO(super_over_noSpin$Geneid,
                             OrgDb = org.Mm.eg.db,
                             keyType = "ENSEMBL",
                             ont = "ALL",
                             pAdjustMethod = "BH",
                             readable = T,
                             pool = F,
                             universe = spin_universe$Geneid) %>%
  as.data.frame() %>%
  mutate(comp = "super_over_noSpin", direction = "super")

go_df_spin <- bind_rows(go_spin_over_super,
                       go_super_over_spin,
                       go_spin_over_noSpin,
                       go_noSpin_over_spin,
                       go_noSpin_over_super,
                       go_super_over_noSpin) %>%
  separate(GeneRatio, c("Num", "Denom")) %>%
  mutate(GeneRatio = as.numeric(Num) / as.numeric(Denom)) %>%
  unite("full_comp", "comp", "direction", remove = F)

# reducing go dataframe
go_df_spin_comps <- go_df_spin %>%
  distinct(full_comp)
go_df_spin_reduced <- data.frame()
for (i in 1:nrow(go_df_spin_comps)) {

  plot_comp <- go_df_spin_comps[i, ]

  plot_go_df <- go_df_spin %>%
    filter(full_comp == plot_comp) %>%
    arrange(qvalue) %>%
    mutate(trans_qvalue = -log10(qvalue))

  if (nrow(plot_go_df) > 1) {

  go_df_reduceTest <- plot_go_df %>%
    group_by(ONTOLOGY) %>%
    summarize(ont_count = n()) %>%
    filter(ont_count == 1)

  rows_to_bind <- plot_go_df %>%
    filter(ONTOLOGY %in% go_df_reduceTest$ONTOLOGY)

  plot_go_df <- plot_go_df %>%
    filter(!(ONTOLOGY %in% go_df_reduceTest$ONTOLOGY))

  go_df_preReduce <- plot_go_df %>%
    select(ONTOLOGY, ID) %>%
    rename(go_type = ONTOLOGY,
           go_id = ID)

  go_df_scores <- plot_go_df$trans_qvalue %>%
    setNames(plot_go_df$ID)

  go_df_reduce <- go_reduce(go_df_preReduce,
                            orgdb = "org.Mm.eg.db",
                            threshold = 0.7,
                            scores = go_df_scores,
                            measure = "Wang")

  plot_go_df_reduce <- plot_go_df %>%
    filter(ID %in% go_df_reduce$parent_id) %>%
    bind_rows(rows_to_bind)

  go_df_spin_reduced <- bind_rows(go_df_spin_reduced, plot_go_df_reduce)

  }

  if (nrow(plot_go_df) < 2) {
    go_df_spin_reduced <- bind_rows(go_df_spin_reduced, plot_go_df)
  }

}

go_df_spin_reduced <- go_df_spin_reduced %>%
  mutate(comp = case_when(comp == "spin_over_super" ~ 
                                 "wGCs_over_Supers",
                               comp == "super_over_spin" ~ 
                                 "Supers_over_wGCs",
                               comp == "spin_over_noSpin" ~ 
                                 "wGCs_over_GCs",
                               comp == "noSpin_over_spin" ~
                                 "GCs_over_wGCs",
                               comp == "noSpin_over_super" ~
                                 "GCs_over_Supers",
                               comp == "super_over_noSpin" ~
                                 "Supers_over_GCs"))
go_df_spin_comps <- go_df_spin_reduced %>%
  distinct(comp)
for (i in 1:nrow(go_df_spin_comps)) {
  
  plot_comp <- go_df_spin_comps[i, ]
  
  plot_go_df <- go_df_spin_reduced %>%
    filter(comp == plot_comp) %>%
    arrange(qvalue) %>%
    group_by(ONTOLOGY) %>%
    slice(1:10) %>%
    ungroup() %>%
    mutate(Description = paste0(ONTOLOGY,
                                ": ",
                                Description),
           Description = fct_reorder2(Description,
                                         -qvalue,
                                         qvalue))
  
  print(ggplot(plot_go_df, aes(x = GeneRatio, y = Description, 
                         fill = qvalue, size = Count)) +
    geom_point(color = "black", pch = 21) +
    facet_wrap(~ ONTOLOGY) + 
    labs(y = "", title = paste0("GO Analysis: ", plot_comp)) +
    theme_classic(base_size = 10, base_family = "sans") +
    scale_fill_distiller(palette = "RdBu", direction = 1) +
    theme(axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          axis.ticks = element_line(color = "black")))
}

# final wGC over GC go plot
wGC_over_gc_go <- ggplot(go_df_spin_reduced %>%
                           filter(comp == "wGCs_over_GCs") %>%
                          arrange(qvalue) %>%
                          slice(1:15) %>%
                          mutate(Description = case_when(Description ==
                                                           "oxidoreduction-driven active transmembrane transporter activity" ~
                                                           "transmembrane transporter activity",
                                                         Description ==
                                                           "inner mitochondrial membrane protein complex" ~
                                                           "inner mitochondrial membrane protein",
                                                         Description ==
                                                           "proton motive force-driven ATP synthesis" ~
                                                           "proton motive ATP synthesis",
                                                         Description ==
                                                           "ribose phosphate biosynthetic process" ~
                                                           "ribose phosphate biosynthesis",
                                                         Description ==
                                                           "NADH dehydrogenase complex assembly" ~
                                                           "NADH dehydrogenase complex",
                                                         T ~ Description)) %>%
                          mutate(Description = fct_reorder2(Description,
                                                            -qvalue, 
                                                            qvalue)), 
                        aes(x = GeneRatio, y = Description, 
                            fill = qvalue, size = Count)) +
  geom_point(color = "black", pch = 21) +
  scale_size(range = c(0.25, 2.5)) +
  labs(y = "", title = "GO Analysis: wGCs Over GCs", size = "Number") +
  theme_classic(base_size = 7, base_family = "sans") +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); wGC_over_gc_go
ggsave("supp_fig_9g.pdf", wGC_over_gc_go, units = "in",
       height = 1.75, width = 3.5)

# final wGC over super go plot
wGC_over_super_go <- ggplot(go_df_spin_reduced %>%
                           filter(comp == "wGCs_over_Supers") %>%
                          arrange(qvalue) %>%
                          slice(1:15) %>%
                          mutate(Description = case_when(Description ==
                                                           "oxidoreduction-driven active transmembrane transporter activity" ~
                                                           "transmembrane transporter activity",
                                                         Description ==
                                                           "inner mitochondrial membrane protein complex" ~
                                                           "inner mitochondrial membrane protein",
                                                         Description ==
                                                           "proton motive force-driven mitochondrial ATP synthesis" ~
                                                           "proton motive mitochondrial ATP synthesis",
                                                         Description ==
                                                           "ATP synthesis coupled electron transport" ~
                                                           "ATP synthesis electron transport",
                                                         Description ==
                                                           "plus-end-directed microtubule motor activity" ~
                                                           "plus-end-directed microtubule motor",
                                                         T ~ Description)) %>%
                          mutate(Description = fct_reorder2(Description,
                                                            -qvalue, 
                                                            qvalue)), 
                        aes(x = GeneRatio, y = Description, 
                            fill = qvalue, size = Count)) +
  geom_point(color = "black", pch = 21) +
  scale_size(range = c(0.25, 2.5)) +
  labs(y = "", title = "GO Analysis: wGCs Over Supers", size = "Number") +
  theme_classic(base_size = 7, base_family = "sans") +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); wGC_over_super_go
ggsave("supp_fig_9h.pdf", wGC_over_super_go, units = "in",
       height = 1.75, width = 3.5)


## ----sessionInfo-------------------------------------------------------------------------------------------------------------------------
sessionInfo()

