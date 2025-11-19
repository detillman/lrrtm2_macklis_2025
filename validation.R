## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(readxl)
library(ggpubr)
library(tidyverse)


## ----validateWB--------------------------------------------------------------------------------------------------------------------------
# load western blot validation data (equal mass per sample by bradford assay)
mass <- read_excel("input_tables.xlsx", sheet = "Input Table 4") %>%
  filter(is.na(Marker) == F) %>%
  mutate(Marker = case_when(Marker == "Nrxn1" ~ "Nrxn", T ~ Marker),
         Marker = as.factor(Marker)) %>%
  mutate(Marker = fct_relevel(Marker, c("Gap43", "Map2", "Lrrtm2", "Nrxn",
                                        "Sema3f", "Nrp2", "PlexA3", "Spire2",
                                        "Camsap3", "Cdk5")))

# calculating significance values for wGCF signal compared to pnh signal
mass_sig <- mass %>%
  filter(Condition != "GCF") %>%
  select(-contains("Percent")) %>%
  pivot_wider(names_from = Condition, values_from = Signal) %>%
  group_by(Marker) %>%
  summarize(pval = t.test(wGCF, PNH, alternative = "greater")$p.value) %>%
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
                      aes(x = as.factor(Marker), y = as.numeric(percentOfPNH),
                          fill = Condition)) +
  geom_bar(stat = "summary", position = "dodge") +
  geom_errorbar(stat='summary', width = 0.25, 
                fun.data = mean_cl_boot, linewidth = 0.5) +
  geom_point(size = 0.5) +
  geom_text(aes(x = Marker, y = 350, label = p_symbol), size = 2) +
  theme_classic(base_size = 7) +
  labs(x = "", y = "GCF WB (% of PNH)", fill = "") +
  scale_fill_manual(values = c("#CA6368")) +
  geom_hline(yintercept = 100, linetype = "dashed", 
             color = "black", linewidth = 0.25) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none"); wb_validate
ggsave("fig_2b.pdf", wb_validate, units = "in", width = 4.41, height = 1.54)

# get means for each candidate GC protein
mass_means <- mass %>%
  filter(Condition != "PNH") %>%
  group_by(Marker, Condition) %>%
  summarize(mean_percent_pnh = mean(percentOfPNH))

# getting representative images for each candidate GC protein
mass_rep <- mass %>%
  filter(Condition != "PNH") %>%
  left_join(mass_means, by = c("Marker", "Condition")) %>%
  mutate(pnh_diff = abs(percentOfPNH - mean_percent_pnh))


## ----validateWgcf------------------------------------------------------------------------------------------------------------------------
# load wgcf coverslip data
data_wgcf <- read_excel("input_tables.xlsx", sheet = "Input Table 5") %>%
  mutate(marker = case_when(marker == "Nrxn1" ~ "Nrxn", T ~ marker),
         marker = case_when(marker == "GfpGfp" ~ "Gfp",
                            T ~ marker),
         marker = fct_relevel(marker, "Gfp", "Gap43", "Map2", "Lrrtm2", "Nrxn",
                              "Sema3f", "Nrp2", "PlexA3",  "Spire2", "Camsap3", 
                              "Cdk5"),
         cpn_with_marker_percent = Q2_count / (Q2_count + Q3_count) * 100,
         marker_with_cpn_percent = Q2_count / (Q2_count + Q1_count) * 100,
         colocalization_score = (cpn_with_marker_percent + marker_with_cpn_percent) / 2)

# plotting gcf coverslip colocalization scores
wgcf_validate <- ggplot(data_wgcf, 
       aes(x = marker, y = colocalization_score)) +
  geom_bar(stat = "summary", fill = "#A96635") +
  geom_errorbar(stat='summary', width = 0.25, fun.data = mean_cl_boot) +
  geom_point(size = 0.5) +
  theme_classic(base_size = 7) +
  stat_compare_means(method = "t.test", ref.group = "Map2", 
                     label = "p.signif", 
                     method.args = list(alternative = "greater"), size = 2) +
  scale_y_continuous(trans = "log10") +
  labs(x = "", y = "CPN Colocalization Score") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); wgcf_validate
ggsave("fig_2f.pdf", wgcf_validate, units = "in", width = 5.76, height = 1.54)


## ----sessionInfo-------------------------------------------------------------------------------------------------------------------------
sessionInfo()

