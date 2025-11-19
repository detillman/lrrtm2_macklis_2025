## ----setup, include=FALSE----------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(readxl)
library(PerformanceAnalytics)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)


## ----qpcrShrna---------------------------------------------------------------------------------------------------------------------------
# loading lrrtm2 shrna data
lrrtm2_shrna <- read_excel("input_tables.xlsx", sheet = "Input Table 6") %>%
  distinct(`Sample Name`, KnockDown) %>%
  separate(`Sample Name`, into = c("Sample", "Replicate")) %>%
  filter(Sample == 2 & is.na(KnockDown) == F) %>%
  mutate(KnockDown = 100 * KnockDown,
         Sample = "shLrrtm2")

# plotting lrrtm2 shrna knockdown
lrrtm2_kd <- ggplot(lrrtm2_shrna, 
       aes(x = Sample, y = KnockDown)) +
  geom_bar(stat = "summary", position = "dodge", fill = "#7C86BB") +
  geom_errorbar(stat = 'summary', width = 0.25, 
                fun.data = mean_cl_boot,
                position = position_dodge(width = 0.9)) +
  geom_point(position = position_dodge(width = 0.9), size = 0.5) +
  theme_classic(base_size = 7) +
  labs(x = "", y = "Lrrtm2 KD Efficiency (%)") +
  # geom_hline(yintercept = 100, linetype = "dashed", 
  #            color = "black", linewidth = 0.25) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none"); lrrtm2_kd
ggsave("supp_fig_3d.pdf", lrrtm2_kd, units = "in",
       height = 1.5, width = 0.87)


## ----loadBlaData-------------------------------------------------------------------------------------------------------------------------
# loading bla data
bla_original <- read_excel("input_tables.xlsx", sheet = "Input Table 7") %>%
  
  # fixing data entry error
  mutate(file = gsub("pTA052_pDT", "pDT", file)) %>%
  mutate(file = gsub("_split", "", file)) %>%
  
  # changing construct names to gene names
  mutate(file = gsub("pTA052", "tdTomato", file)) %>%
  mutate(file = gsub("pDT016e_", "cytLrrtm2_", file)) %>%
  mutate(file = gsub("pDT017me", "Nrxn3aMS4", file)) %>%
  mutate(file = gsub("pDT016es", "memLrrtm2", file)) %>%
  mutate(file = gsub("pDT021pt5", "shLrrtm2", file)) %>%
  mutate(file = gsub("SAC002", "shControl", file)) %>%
  
  # replacing sexed brains with generic brains
  mutate(file = gsub("215_memLrrtm2_p21_female1", "215_memLrrtm2_p21_brain1", file)) %>%
  mutate(file = gsub("215_memLrrtm2_p21_female2", "215_memLrrtm2_p21_brain2", file)) %>%
  mutate(file = gsub("215_memLrrtm2_p21_male1", "215_memLrrtm2_p21_brain3", file)) %>%
  mutate(file = gsub("215_memLrrtm2_p21_male2", "215_memLrrtm2_p21_brain4", file)) %>%
  mutate(file = gsub("219_memLrrtm2_p21_male1", "219_memLrrtm2_p21_brain1", file)) %>%
  
  mutate(file = gsub("242_shLrrtm2_p21_male1", "242_shLrrtm2_p21_brain1", file)) %>%
  mutate(file = gsub("242_shLrrtm2_p21_male2", "242_shLrrtm2_p21_brain2", file)) %>%
  mutate(file = gsub("244_shLrrtm2_p21_female1", "244_shLrrtm2_p21_brain1", file)) %>%
  mutate(file = gsub("244_shLrrtm2_p21_male1", "244_shLrrtm2_p21_brain2", file)) %>%
  mutate(file = gsub("253_shLrrtm2_p21_male1", "253_shLrrtm2_p21_brain1", file)) %>%
  
  mutate(file = gsub("242_shControl_p21_female1", "242_shControl_p21_brain1", file)) %>%
  mutate(file = gsub("242_shControl_p21_female2", "242_shControl_p21_brain2", file)) %>%
  mutate(file = gsub("242_shControl_p21_female3", "242_shControl_p21_brain3", file)) %>%
  mutate(file = gsub("244_shControl_p21_female1", "244_shControl_p21_brain1", file)) %>%
  mutate(file = gsub("244_shControl_p21_male1", "244_shControl_p21_brain2", file)) %>%
  mutate(file = gsub("244_shControl_p21_male2", "244_shControl_p21_brain3", file)) %>%
  mutate(file = gsub("253_shControl_p21_female1", "253_shControl_p21_brain1", file)) %>%
  
  # getting image info
  separate(file, into = c("nb", "exp", "construct", "age", 
                          "brain", "section", "zoom"), sep = "_") %>%
  select(-nb, -zoom) %>%
  
  # making sample columns
  unite(col = "sample", "exp", "construct", "age", "brain", "section",
        remove = F) %>%
  unite(col = "combo", "exp", "construct", "age", "brain", remove = F)


# calculating bla metrics
bla <- bla_original %>%
  mutate(total_intensity = mean * thresholdArea,
         exp = as.numeric(exp)) %>%
  
  # replace NaN with zero
  mutate(mean = case_when(is.na(mean) ~ 0, 
                          T ~ mean),
         thresholdMean = case_when(is.na(thresholdMean) ~ 0, 
                                   T ~ thresholdMean),
         thresholdFraction = case_when(is.na(thresholdFraction) ~ 0, 
                                       T ~ thresholdFraction),
         thresholdNumber = case_when(is.na(thresholdNumber) ~ 0, 
                                     T ~ thresholdNumber),
         total_intensity = case_when(is.na(total_intensity) ~ 0, 
                                     T ~ total_intensity)) %>%
  mutate(totalBLAtotal = contraBLAtotal + ipsiBLAtotal,
         totalBLAdensity = contraBLAdensity + ipsiBLAdensity,
         totalBLAmean = contraBLAmean + ipsiBLAmean,
         totalBLAtotal_thresh = contraBLAtotal_thresh + ipsiBLAtotal_thresh,
         totalBLAdensity_thresh = contraBLAdensity_thresh + ipsiBLAdensity_thresh,
         totalBLAmean_thresh = contraBLAmean_thresh + ipsiBLAmean_thresh) %>%
  mutate(construct = fct_relevel(construct, "tdTomato", "memLrrtm2",
                                 "Nrxn3aMS4", "cytLrrtm2", "shControl",
                                 "shLrrtm2"))

# plotting number of brains per construct
bla_metadata <- bla %>%
  distinct(exp, construct, brain) %>%
  group_by(construct) %>%
  summarize(number_brains = n())
ggplot(bla_metadata, aes(x = construct, y = number_brains)) +
  geom_col(pos = "dodge", fill = "black", color = "black") +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black", 
                                   angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank())


## ----blaInnervateQuantification----------------------------------------------------------------------------------------------------------
# calculating measurement sums for each ipsi section
bla_iue_sums <- bla %>%
  
  # restrict to less than 9000 microns
  filter(distance_microns < 9000) %>%
  group_by(sample) %>%
  summarize(total_mean = sum(mean),
            total_totalIntensity = sum(total_intensity),
            total_threshFraction = sum(thresholdFraction),
            total_threshNumber = sum(thresholdNumber))

# calculating measurement percentages of each ipsi box per section
bla_iue <- bla %>%
  
  # restrict to less than 9000 microns
  filter(distance_microns < 9000) %>%
  left_join(bla_iue_sums, by = "sample") %>%
  mutate(percent_mean = mean / total_mean * 100,
         percent_totalIntensity = total_intensity / total_totalIntensity * 100,
         percent_fraction = thresholdFraction / total_threshFraction * 100,
         percent_number = thresholdNumber / total_threshNumber * 100)

# correlating intensity metrics
bla_cor <- bla_iue %>%
  select(sample, contains("BLA"), contains("total_")) %>%
  select(-total_intensity, -contains("contra"), -contains("ipsi")) %>%
  distinct() %>%
  column_to_rownames("sample")
chart.Correlation(bla_cor)

# determining bla innervation metric by comparing tdtomato OE and Bcl11a-Null
bla_metric <- bla %>%
  filter(construct == "tdTomato" | construct == "shControl") %>%
  select(sample, construct, contains("totalBLA")) %>%
  distinct() %>%
  pivot_longer(cols = contains("totalBLA"), 
               names_to = "bla_metric", values_to = "bla_value")
ggplot(bla_metric, aes(x = bla_metric, y = bla_value, fill = construct)) +
  geom_bar(stat = "summary", position = "dodge") +
  geom_errorbar(stat = 'summary', width = 0.25, 
                fun.data = mean_cl_boot,
                position = position_dodge(width = 0.9)) +
  geom_point(position = position_dodge(width = 0.9)) +
  stat_compare_means(label = "p.format", size = 2) +
  theme_classic() +
  facet_wrap(~ bla_metric, scales = "free") +
  scale_fill_manual(values = c("#E69F00", "#D55E00")) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none")
## using mean because most stark difference between tdTomato OE and Bcl11a-Null

# plotting all ipsi brains
ggplot(bla_iue, 
       aes(x = distance_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              span = 0.1,
              aes(group = combo, color = combo, fill = combo)) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  labs(x = "Distance from midline (microns)",
       y = "Mean Intensity (%)") +
  facet_wrap(~ construct, scales = "free")

# plotting tdTomato iues
ggplot(bla_iue %>%
         filter(construct == "tdTomato"), 
       aes(x = distance_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              span = 0.1,
              aes(group = combo, color = combo, fill = combo)) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) +
  labs(x = "Distance from midline (microns)",
       y = "Mean Intensity (%)")

# removing lateral tdTomato iues
bla_iue <- bla_iue %>%
  filter(combo != "176_tdTomato_p21_brain1") %>% # lateral iue
  filter(combo != "181_tdTomato_p21_brain2") %>% # lateral iue
  filter(combo != "185_tdTomato_p21_brain1") # lateral iue

# distinct bla samples
bla_constructs <- bla_iue %>%
  distinct(construct)

# plotting all constructs with tdTomato average to identify medial IUEs
for (i in 1:nrow(bla_constructs)) {
  
  # select appropriate construct
  plot_construct <- bla_constructs %>%
    slice(i)
  
  # filter data for appropriate construct
  plot_bla <- bla_iue %>%
    filter(construct %in% plot_construct$construct)
  
  # plot electroporation site
  print(ggplot(plot_bla, aes(x = distance_microns, 
                                 y = percent_mean)) +
          stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
                      span = 0.1,
                      aes(group = combo, color = combo, fill = combo)) +
          stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
                      
                      # plotting tdTomato data
                      data = bla_iue %>%
                        filter(construct == "tdTomato") %>%  
                        transform(combo = NULL),
                      span = 0.1,
                      aes(group = construct, color = construct, fill = construct)) +
          facet_wrap(~ combo, scales = "free") +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black"),
                axis.text.y = element_text(color = "black")) +
          labs(x = "Distance from midline (microns)",
               y = "Mean Intensity (%)"))
}

# identifying brains with medial and/or lateral iues
bla_iue_medial <- bla_iue %>%
  filter(combo != "176_Nrxn3aMS4_p21_brain2") %>% # lateral iue

  filter(combo != "176_cytLrrtm2_p21_brain2") %>% # lateral iue
  
  filter(combo != "242_shControl_p21_brain1") %>% # lateral iue
  filter(combo != "253_shControl_p21_brain1") %>% # lateral peaks
  
  filter(combo != "242_shLrrtm2_p21_brain2") %>% # lateral peaks
  filter(combo != "244_shLrrtm2_p21_brain1") # lateral iue

# plotting iues with medial electroporations
bla_iue <- ggplot(bla_iue_medial, 
                        aes(x = distance_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              span = 0.1,
              aes(group = construct, color = construct, fill = construct)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#009E73", "#56B4E9", "#cc79a7", "#7C86BB")) +
  scale_color_manual(values = c("#D55E00", "#0072B2", "#009E73", "#56B4E9", "#cc79a7", "#7C86BB")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.8)) +
  labs(x = "Electroporation: um From Midline",
       y = "% Mean Intensity",
       fill = "", color = ""); bla_iue
ggsave("supp_fig_3e.pdf", bla_iue, units = "in",
       height = 1.5, width = 2.90)

# calculating bla innervation intensity based on mean
bla_innervate <- bla_iue_medial %>%
  distinct(sample, combo, construct, ipsiBLAmean, contraBLAmean,
           totalBLAmean, total_mean) %>%
  group_by(combo, construct) %>%
  summarize(ipsi_final = sum(ipsiBLAmean),
            contra_final = sum(contraBLAmean),
            total_final = sum(totalBLAmean),
            iue_final = sum(total_mean)) %>%
  pivot_longer(cols = contains("final"), 
               names_to = "metric", values_to = "value") %>%
  ungroup()
ggplot(bla_innervate %>%
         filter(metric != "iue_final"), 
       aes(x = metric, y = value, fill = construct)) +
  geom_bar(stat = "summary", position = "dodge") +
  geom_errorbar(stat = 'summary', width = 0.25, 
                fun.data = mean_cl_boot,
                position = position_dodge(width = 0.9)) +
  geom_point(position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#009E73", "#56B4E9", "#cc79a7", "#7C86BB")) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))
ggplot(bla_innervate %>%
         filter(metric == "total_final"), 
       aes(x = construct, y = value, fill = construct)) +
  geom_bar(stat = "summary", position = "dodge") +
  geom_errorbar(stat = 'summary', width = 0.25, 
                fun.data = mean_cl_boot,
                position = position_dodge(width = 0.9)) +
  geom_point(position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#009E73", "#56B4E9", "#cc79a7", "#7C86BB")) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))

# normalizing by residuals of mean iue intensity
bla_resid_plot <- bla_innervate %>%
  pivot_wider(names_from = "metric", values_from = "value")
ggplot(bla_resid_plot, aes(x = iue_final, y = ipsi_final)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method='lm', formula = 'y ~ x', color = "black") +
  stat_regline_equation(label.y = 5*10^3, color = "black",
                        aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 6*10^3, color = "black",
                        aes(label = ..rr.label..))
ggplot(bla_resid_plot %>%
         filter(construct == "tdTomato"), 
       aes(x = iue_final, y = ipsi_final)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method='lm', formula = 'y ~ x', color = "black") +
  stat_regline_equation(label.y = 5*10^3, color = "black",
                        aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 6*10^3, color = "black",
                        aes(label = ..rr.label..))
ggplot(bla_resid_plot %>%
         filter(construct == "tdTomato"), 
       aes(x = iue_final, y = contra_final)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method='lm', formula = 'y ~ x', color = "black") +
  stat_regline_equation(label.y = 5*10^3, color = "black",
                        aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 6*10^3, color = "black",
                        aes(label = ..rr.label..))
ggplot(bla_resid_plot %>%
         filter(construct == "tdTomato"), 
       aes(x = iue_final, y = total_final)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method='lm', formula = 'y ~ x', color = "black") +
  stat_regline_equation(label.y = 5*10^3, color = "black",
                        aes(label = ..eq.label..)) +
  stat_regline_equation(label.y = 6*10^3, color = "black",
                        aes(label = ..rr.label..))

# calculating bla innervation intensity based on tdTomato residuals
bla_resid_final <- bla_resid_plot %>%
  mutate(ipsi_expect = 700 + 0.00024 * iue_final,
         ipsi_resid = ipsi_final / ipsi_expect * 100,
         contra_expect = 570 + 3.4*10^-5 * iue_final,
         contra_resid = contra_final / contra_expect * 100,
         total_expect = 1300 + 0.00028 * iue_final,
         total_resid = total_final / total_expect * 100) %>%
  select(-contains("final"), -contains("expect")) %>%
  pivot_longer(cols = contains("resid"), 
               names_to = "metric", values_to = "value")
ggplot(bla_resid_final, 
       aes(x = metric, y = value, fill = construct)) +
  geom_bar(stat = "summary", position = "dodge") +
  geom_errorbar(stat = 'summary', width = 0.25, 
                fun.data = mean_cl_boot,
                position = position_dodge(width = 0.9)) +
  geom_point(position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#009E73", 
                               "#56B4E9", "#cc79a7", "#7C86BB")) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))
bla_innervation <- ggplot(bla_resid_final %>% filter(metric == "total_resid"), 
       aes(x = construct, y = value, fill = construct)) +
  geom_bar(stat = "summary", position = "dodge") +
  geom_errorbar(stat = 'summary', width = 0.25, 
                fun.data = mean_cl_boot,
                position = position_dodge(width = 0.9)) +
  geom_point(position = position_dodge(width = 0.9), size = 0.5) +
  scale_fill_manual(values = c("#D55E00", "#0072B2", "#009E73", 
                               "#56B4E9", "#cc79a7", "#7C86BB")) +
  stat_compare_means(method = "t.test", ref.group = "tdTomato", 
                     label = "p.signif", 
                     method.args = list(alternative = "two.sided"), size = 2) +
  theme_classic(base_size = 7) +
  labs(x = "", y = "Normalized BLA Intensity (%)") +
  geom_hline(yintercept = 100, linetype = "dashed", 
             color = "black", linewidth = 0.25) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none"); bla_innervation
ggsave("fig_3c.pdf", bla_innervation, units = "in", height = 2.16, width = 3.44)


## ----loadCortexData----------------------------------------------------------------------------------------------------------------------
# loading cortex data
cortex <- read_excel("input_tables.xlsx", sheet = "Input Table 8") %>%
  
  # fixing data entry error
  mutate(file = gsub("pTA052_pDT", "pDT", file)) %>%
  
  # changing construct names to gene names
  mutate(file = gsub("pTA052", "tdTomato", file)) %>%
  mutate(file = gsub("pDT016e", "Lrrtm2", file)) %>%
  mutate(file = gsub("pDT017me", "Nrxn3aMS4", file)) %>%
  mutate(file = gsub("pDT018e", "Sema3f", file)) %>%
  mutate(file = gsub("pDT019e", "Nrp2", file)) %>%
  mutate(file = gsub("pDT020e", "PlexA3", file)) %>%
  
  # getting image info
  separate(file, into = c("nb", "exp", "construct", "age", 
                          "brain", "section", "zoom"), sep = "_") %>%
  select(-nb, -zoom) %>%
  
  # making sample columns
  unite(col = "sample", "exp", "construct", "age", "brain", "section",
        remove = F) %>%
  unite(col = "combo", "exp", "construct", "age", "brain", remove = F) %>%
  mutate(construct = fct_relevel(construct, "tdTomato", "Lrrtm2", "Nrxn3aMS4",
                                 "Sema3f", "Nrp2", "PlexA3"),
         exp = as.numeric(exp)) %>%
  mutate(total_intensity = mean * thresholdArea) %>%
  replace(is.na(.), 0)

# plotting number of brains per construct
cortex_metadata <- cortex %>%
  distinct(exp, construct, brain) %>%
  group_by(construct) %>%
  summarize(number_brains = n())
ggplot(cortex_metadata, aes(x = construct, y = number_brains)) +
  geom_col(pos = "dodge", fill = "black", color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", 
                                   angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank())


## ----cortexInnervateLrrtm2---------------------------------------------------------------------------------------------------------------
# filtering for lrrtm2 brains (all cortical layers)
cortex_lrrtm2 <- cortex %>%
  filter(exp < 178)

# plotting number of lrrtm2 brains per construct
cortex_lrrtm2_metadata <- cortex_lrrtm2 %>%
  distinct(exp, construct, brain) %>%
  group_by(construct) %>%
  summarize(number_brains = n())
ggplot(cortex_lrrtm2_metadata, aes(x = construct, y = number_brains)) +
  geom_col(pos = "dodge", fill = "black", color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", 
                                   angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank())

# calculating measurement sums for each ipsi section
cortex_iue_lrrtm2_sums <- cortex_lrrtm2 %>%
  
  # restrict to less than 6000 microns
  filter(distance_microns < 6000) %>%
  filter(feature == "ipsi") %>%
  group_by(sample) %>%
  summarize(total_mean = sum(mean),
            total_totalIntensity = sum(total_intensity),
            total_threshFraction = sum(thresholdFraction),
            total_threshNumber = sum(thresholdNumber))

# calculating measurement percentages of each ipsi box per section
cortex_iue_lrrtm2 <- cortex_lrrtm2 %>%
  
  # restrict to less than 6000 microns
  filter(distance_microns < 6000) %>%
  filter(feature == "ipsi") %>%
  left_join(cortex_iue_lrrtm2_sums, by = "sample") %>%
  mutate(percent_mean = mean / total_mean * 100,
         percent_totalIntensity = total_intensity / total_totalIntensity * 100,
         percent_fraction = thresholdFraction / total_threshFraction * 100,
         percent_number = thresholdNumber / total_threshNumber * 100)

# plotting all ipsi brains
ggplot(cortex_iue_lrrtm2, 
       aes(x = distance_microns, y = percent_number)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              span = 0.1,
              aes(group = combo, color = combo, fill = combo)) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  labs(x = "Distance from midline (microns)",
       y = "Pixels Above Threshold (%)") +
  facet_wrap(~ construct, scales = "free")

# calculating measurement sums for each contra section
cortex_cortex_lrrtm2_sums <- cortex_lrrtm2 %>%
  
  # restrict to less than 6000 microns
  filter(distance_microns < 6000) %>%
  filter(feature == "contra") %>%
  group_by(sample) %>%
  summarize(total_mean = sum(mean),
            total_totalIntensity = sum(total_intensity),
            total_threshFraction = sum(thresholdFraction),
            total_threshNumber = sum(thresholdNumber))

# calculating measurement percentages of each contra box per section
cortex_cortex_lrrtm2 <- cortex_lrrtm2 %>%
  
  # restrict to less than 6000 microns
  filter(distance_microns < 6000) %>%
  filter(feature == "contra") %>%
  left_join(cortex_cortex_lrrtm2_sums, by = "sample") %>%
  mutate(percent_mean = mean / total_mean * 100,
         percent_totalIntensity = total_intensity / total_totalIntensity * 100,
         percent_fraction = thresholdFraction / total_threshFraction * 100,
         percent_number = thresholdNumber / total_threshNumber * 100)

# plotting all contra brains
ggplot(cortex_cortex_lrrtm2, 
       aes(x = distance_microns, y = percent_number)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              span = 0.05,
              aes(group = combo, color = combo, fill = combo)) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  labs(x = "Distance from midline (microns)",
       y = "Pixels Above Threshold (%)") +
  facet_wrap(~ construct, scales = "free")

# distinct cortex samples
cortex_lrrtm2_constructs <- cortex_iue_lrrtm2 %>%
  distinct(construct)

# plotting all constructs to identify bilateral IUEs
for (i in 1:nrow(cortex_lrrtm2_constructs)) {

  # select appropriate construct
  plot_construct <- cortex_lrrtm2_constructs %>%
    slice(i)

  # filter data for appropriate construct
  plot_cortex_lrrtm2 <- cortex_cortex_lrrtm2 %>%
    filter(construct %in% plot_construct$construct)

  # plot cortical innervation
  print(ggplot(plot_cortex_lrrtm2, aes(x = distance_microns, y = percent_number)) +
    stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
                span = 0.05,
                aes(group = combo, color = combo, fill = combo)) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black")) +
    labs(x = "Distance from midline (microns)",
         y = "Pixels Above Threshold (%)"))
}

# removing bilateral iues
cortex_cortex_lrrtm2 <- cortex_cortex_lrrtm2 %>%
  filter(combo != "177_Lrrtm2_p7_brain2",
         combo != "174_Lrrtm2_p7_brain1",
         combo != "175_Lrrtm2_p7_brain2",
         
         combo != "172_Nrxn3aMS4_p7_brain1",
         combo != "171_Nrxn3aMS4_p7_brain2",
         combo != "171_Nrxn3aMS4_p7_brain3",
         
         combo != "176_tdTomato_p7_brain3",
         combo != "176_tdTomato_p7_brain1")

# plotting all constructs with tdTomato average to identify medial IUEs
for (i in 1:nrow(cortex_lrrtm2_constructs)) {
  
  # select appropriate construct
  plot_construct <- cortex_lrrtm2_constructs %>%
    slice(i)
  
  # filter data for appropriate construct
  plot_cortex_lrrtm2 <- cortex_iue_lrrtm2 %>%
    filter(construct %in% plot_construct$construct)
  
  # plot electoroporation site
  print(ggplot(plot_cortex_lrrtm2, aes(x = distance_microns, y = percent_number)) +
          stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
                      span = 0.1,
                      aes(group = combo, color = combo, fill = combo)) +
          stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
                      
                      # plotting tdTomato data
                      data = cortex_iue_lrrtm2 %>%
                        filter(construct == "tdTomato") %>%  
                        transform(combo = NULL),
                      span = 0.1,
                      aes(group = construct, color = construct, fill = construct)) +
          facet_wrap(~ combo, scales = "free") +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black"),
                axis.text.y = element_text(color = "black")) +
          labs(x = "Distance from midline (microns)",
               y = "Pixels Above Threshold (%)"))
}

# identifying brains with medial iues
cortex_iue_lrrtm2_medial <- cortex_iue_lrrtm2 %>%
  mutate(medial = case_when(construct == "tdTomato" ~ T,
                            
                            combo == "170_Lrrtm2_p7_brain2" ~ T,
                            combo == "172_Lrrtm2_p7_brain1" ~ T,
                            combo == "174_Lrrtm2_p7_brain1" ~ T,
                            combo == "175_Lrrtm2_p7_brain1" ~ T,
                            combo == "175_Lrrtm2_p7_brain2" ~ T,
                            combo == "176_Lrrtm2_p7_brain2" ~ T,
                            
                            combo == "171_Nrxn3aMS4_p7_brain2" ~ T,
                            combo == "171_Nrxn3aMS4_p7_brain3" ~ T,
                            combo == "172_Nrxn3aMS4_p7_brain1" ~ T,
                            combo == "175_Nrxn3aMS4_p7_brain2" ~ T,
                            combo == "175_Nrxn3aMS4_p7_brain3" ~ T,
                            combo == "176_Nrxn3aMS4_p7_brain1" ~ T,
                            
                            T ~ F)) %>%
  filter(medial == T)

cortex_cortex_lrrtm2_medial <- cortex_cortex_lrrtm2 %>%
  filter(combo %in% cortex_iue_lrrtm2_medial$combo)

# plotting iues with medial electroporations
lrrtm2_lrrtm2_iue <- ggplot(cortex_iue_lrrtm2_medial %>%
                          filter(construct != "Nrxn3aPS4") %>%
                          mutate(construct = case_when(construct == "Lrrtm2" ~ 
                                                         "cytLrrtm2",
                                                       T ~ construct),
                                 construct = fct_relevel(construct, "tdTomato",
                                                         "cytLrrtm2", "Nrxn3aMS4")), 
                        aes(x = distance_microns, y = percent_number)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              span = 0.1,
              aes(group = construct, color = construct, fill = construct)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#56B4E9", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.8)) +
   labs(x = "Electroporation: um From Midline",
       y = "% Thresholded",
       fill = "", color = ""); lrrtm2_lrrtm2_iue
ggsave("supp_fig_3f.pdf", lrrtm2_lrrtm2_iue, units = "in",
       height = 1.5, width = 2.90)

# plotting iues with medial electroporations
lrrtm2_lrrtm2_cortex <- ggplot(cortex_cortex_lrrtm2_medial %>%
                          filter(construct != "Nrxn3aPS4") %>%
                          mutate(construct = case_when(construct == "Lrrtm2" ~ 
                                                         "cytLrrtm2",
                                                       T ~ construct),
                                 construct = fct_relevel(construct, "tdTomato",
                                                         "cytLrrtm2", "Nrxn3aMS4")), 
                        aes(x = distance_microns, y = percent_number)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              span = 0.05,
              aes(group = construct, color = construct, fill = construct)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#56B4E9", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#56B4E9", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.1, 0.8)) +
   labs(x = "Innervation: um from midline",
       y = "% Thresholded",
       fill = "", color = ""); lrrtm2_lrrtm2_cortex
ggsave("fig_3e.pdf", lrrtm2_lrrtm2_cortex, units = "in",
       height = 1.33, width = 1.69)

# finding representative lrrtm2 images
test_means_lrrtm2 <- cortex_cortex_lrrtm2_medial %>%
  group_by(construct, box) %>%
  summarize(mean_percent = mean(percent_number, na.rm = T))

test_rep_lrrtm2 <- cortex_cortex_lrrtm2_medial %>%
  left_join(test_means_lrrtm2, by = c("construct", "box")) %>%
  mutate(mean_diff = abs(mean_percent - percent_number)) %>%
  group_by(sample) %>%
  summarize(sum_diff = sum(mean_diff, na.rm = T))


## ----cortexInnervateSema3f---------------------------------------------------------------------------------------------------------------
# filtering for sema3f brains (superficial cortical layers)
cortex_sema3f <- cortex %>%
  filter(exp > 178)

# plotting number of sema3f brains per construct
cortex_sema3f_metadata <- cortex_sema3f %>%
  distinct(exp, construct, brain) %>%
  group_by(construct) %>%
  summarize(number_brains = n())
ggplot(cortex_sema3f_metadata, aes(x = construct, y = number_brains)) +
  geom_col(pos = "dodge", fill = "black", color = "black") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", 
                                   angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank())

# calculating measurement sums for each ipsi section
cortex_iue_sema3f_sums <- cortex_sema3f %>%
  
  # restrict to less than 6000 microns
  filter(distance_microns < 6000) %>%
  filter(feature == "ipsi") %>%
  group_by(sample) %>%
  summarize(total_mean = sum(mean),
            total_totalIntensity = sum(total_intensity),
            total_threshFraction = sum(thresholdFraction),
            total_threshNumber = sum(thresholdNumber))

# calculating measurement percentages of each ipsi box per section
cortex_iue_sema3f <- cortex_sema3f %>%
  
  # restrict to less than 6000 microns
  filter(distance_microns < 6000) %>%
  filter(feature == "ipsi") %>%
  left_join(cortex_iue_sema3f_sums, by = "sample") %>%
  mutate(percent_mean = mean / total_mean * 100,
         percent_totalIntensity = total_intensity / total_totalIntensity * 100,
         percent_fraction = thresholdFraction / total_threshFraction * 100,
         percent_number = thresholdNumber / total_threshNumber * 100)

# plotting all ipsi brains
ggplot(cortex_iue_sema3f, 
       aes(x = distance_microns, y = percent_number)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              span = 0.1,
              aes(group = combo, color = combo, fill = combo)) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  labs(x = "Distance from midline (microns)",
       y = "Pixels Above Threshold (%)") +
  facet_wrap(~ construct, scales = "free")

# calculating measurement sums for each contra section
cortex_cortex_sema3f_sums <- cortex_sema3f %>%
  
  # restrict to less than 6000 microns
  filter(distance_microns < 6000) %>%
  filter(feature == "contra") %>%
  group_by(sample) %>%
  summarize(total_mean = sum(mean),
            total_totalIntensity = sum(total_intensity),
            total_threshFraction = sum(thresholdFraction),
            total_threshNumber = sum(thresholdNumber))

# calculating measurement percentages of each contra box per section
cortex_cortex_sema3f <- cortex_sema3f %>%
  
  # restrict to less than 6000 microns
  filter(distance_microns < 6000) %>%
  filter(feature == "contra") %>%
  left_join(cortex_cortex_sema3f_sums, by = "sample") %>%
  mutate(percent_mean = mean / total_mean * 100,
         percent_totalIntensity = total_intensity / total_totalIntensity * 100,
         percent_fraction = thresholdFraction / total_threshFraction * 100,
         percent_number = thresholdNumber / total_threshNumber * 100)

# plotting all contra brains
ggplot(cortex_cortex_sema3f, 
       aes(x = distance_microns, y = percent_number)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              span = 0.05,
              aes(group = combo, color = combo, fill = combo)) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  labs(x = "Distance from midline (microns)",
       y = "Pixels Above Threshold (%)") +
  facet_wrap(~ construct, scales = "free")

# distinct cortex samples
cortex_sema3f_constructs <- cortex_iue_sema3f %>%
  distinct(construct)

# plotting all constructs to identify bilateral IUEs
for (i in 1:nrow(cortex_sema3f_constructs)) {

  # select appropriate construct
  plot_construct <- cortex_sema3f_constructs %>%
    slice(i)

  # filter data for appropriate construct
  plot_cortex_sema3f <- cortex_cortex_sema3f %>%
    filter(construct %in% plot_construct$construct)

  # plot cortical innervation
  print(ggplot(plot_cortex_sema3f, aes(x = distance_microns, y = percent_number)) +
    stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
                span = 0.05,
                aes(group = combo, color = combo, fill = combo)) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black")) +
    labs(x = "Distance from midline (microns)",
         y = "Pixels Above Threshold (%)"))
}

# removing bilateral iues
cortex_cortex_sema3f <- cortex_cortex_sema3f %>%
  filter(combo != "180_tdTomato_p7_brain3",
         combo != "187_tdTomato_p7_brain1",
         combo != "184_Nrp2_p7_brain4",
         combo != "181_Nrp2_p7_brain2",
         combo != "184_Nrp2_p7_brain6",
         combo != "184_PlexA3_p7_brain1",
         combo != "185_PlexA3_p7_brain1")
cortex_iue_sema3f <- cortex_iue_sema3f %>%
  filter(combo != "180_tdTomato_p7_brain3",
         combo != "187_tdTomato_p7_brain1",
         combo != "184_Nrp2_p7_brain4",
         combo != "181_Nrp2_p7_brain2",
         combo != "184_Nrp2_p7_brain6",
         combo != "184_PlexA3_p7_brain1",
         combo != "185_PlexA3_p7_brain1")

# plotting all constructs with tdTomato average to identify medial IUEs
for (i in 1:nrow(cortex_sema3f_constructs)) {
  
  # select appropriate construct
  plot_construct <- cortex_sema3f_constructs %>%
    slice(i)
  
  # filter data for appropriate construct
  plot_cortex_sema3f <- cortex_iue_sema3f %>%
    filter(construct %in% plot_construct$construct)
  
  # plot electoroporation site
  print(ggplot(plot_cortex_sema3f, aes(x = distance_microns, y = percent_number)) +
          stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
                      span = 0.1,
                      aes(group = combo, color = combo, fill = combo)) +
          stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
                      
                      # plotting tdTomato data
                      data = cortex_iue_sema3f %>%
                        filter(construct == "tdTomato") %>%  
                        transform(combo = NULL),
                      span = 0.1,
                      aes(group = construct, color = construct, fill = construct)) +
          facet_wrap(~ combo, scales = "free") +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black"),
                axis.text.y = element_text(color = "black")) +
          labs(x = "Distance from midline (microns)",
               y = "Pixels Above Threshold (%)"))
}

# identifying brains with medial iues
cortex_iue_sema3f_medial <- cortex_iue_sema3f %>%
  mutate(medial = case_when(construct == "tdTomato" ~ T,
                            
                            combo == "185_Sema3f_p7_brain2" ~ T,
                            combo == "187_Sema3f_p7_brain1" ~ T,
                            combo == "187_Sema3f_p7_brain3" ~ T,
                            
                            combo == "181_Nrp2_p7_brain1" ~ T,
                            combo == "184_Nrp2_p7_brain1" ~ T,
                            combo == "185_Nrp2_p7_brain1" ~ T,
                            
                            combo == "187_PlexA3_p7_brain1" ~ T,
                            combo == "187_PlexA3_p7_brain2" ~ T,
                            combo == "188_PlexA3_p7_brain2" ~ T,
                            
                            T ~ F)) %>%
  filter(medial == T)

cortex_cortex_sema3f_medial <- cortex_cortex_sema3f %>%
  filter(combo %in% cortex_iue_sema3f_medial$combo)

# plotting iues with medial electroporations
sema3f_sema3f_iue <- ggplot(cortex_iue_sema3f_medial, 
                        aes(x = distance_microns, y = percent_number)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              span = 0.1,
              aes(group = construct, color = construct, fill = construct)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#2B93CE", "#80CFB9", "#004F3A")) +
  scale_color_manual(values = c("#D55E00", "#2B93CE", "#80CFB9", "#004F3A")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.8)) +
  labs(x = "Electroporation: um From Midline",
       y = "% Thresholded",
       fill = "", color = ""); sema3f_sema3f_iue
ggsave("supp_fig_4a_ii.pdf", sema3f_sema3f_iue,
       units = "in", height = 1.04, width = 2.94)

# plotting iues with medial electroporations
sema3f_sema3f_cortex <- ggplot(cortex_cortex_sema3f_medial, 
                        aes(x = distance_microns, y = percent_number)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              span = 0.05,
              aes(group = construct, color = construct, fill = construct)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#2B93CE", "#80CFB9", "#004F3A")) +
  scale_color_manual(values = c("#D55E00", "#2B93CE", "#80CFB9", "#004F3A")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.8)) +
  labs(x = "Innervation: um From Midline",
       y = "% Thresholded",
       fill = "", color = ""); sema3f_sema3f_cortex
ggsave("supp_fig_4a_i.pdf", sema3f_sema3f_cortex,
       units = "in", height = 0.97, width = 2.94)

test_means_sema3f <- cortex_cortex_sema3f_medial %>%
  group_by(construct, box) %>%
  summarize(mean_percent = mean(percent_number, na.rm = T))

test_rep_sema3f <- cortex_cortex_sema3f_medial %>%
  left_join(test_means_sema3f, by = c("construct", "box")) %>%
  mutate(mean_diff = abs(mean_percent - percent_number)) %>%
  group_by(sample) %>%
  summarize(sum_diff = sum(mean_diff, na.rm = T))


## ----loadExtendData----------------------------------------------------------------------------------------------------------------------
# loading axon data
axon_original <- read_excel("input_tables.xlsx", sheet = "Input Table 9") %>%
  
  # fixing data entry error
  mutate(file = gsub("pTA052_pDT", "pDT", file),
         file = gsub("_split", "", file)) %>%
  
  # changing construct names to gene names
  mutate(file = gsub("pTA052", "tdTomato", file)) %>%
  mutate(file = gsub("pDT018e", "Sema3f", file)) %>%
  mutate(file = gsub("pDT019e", "Nrp2", file)) %>%
  mutate(file = gsub("pDT020e", "PlexA3", file)) %>%
  mutate(file = gsub("pDT018eB", "Sema3f", file)) %>%
  mutate(file = gsub("pDT019eB", "Nrp2", file)) %>%
  mutate(file = gsub("pDT020eB", "PlexA3", file)) %>%
  
  # getting image info
  separate(file, into = c("nb", "exp", "construct", "age", 
                          "brain", "section", "zoom"), sep = "_") %>%
  select(-nb, -zoom) %>%
  
  
  # making beam color column
  separate("feature", into = c("color", "comp"), fill = "left", remove = F) %>%
  mutate(color = case_when(is.na(color) == T ~ "red", T ~ color)) %>%
  
  # making sample columns
  unite(col = "sample", "exp", "construct", "age", "brain", "section",
        remove = F) %>%
  unite(col = "combo", "exp", "construct", "age", "brain", remove = F) %>%
  mutate(construct = fct_relevel(construct, "tdTomato", "Sema3f", 
                                 "Nrp2", "PlexA3", "Sema3fB", "Nrp2B",
                                 "PlexA3B"),
         exp = as.numeric(exp),
         mean = as.numeric(mean)) %>%
  mutate(total_intensity = mean * thresholdArea)

# remove midline rows with nans
axon_samples <- axon_original %>%
  distinct(sample)
axon_betterROIs <- data.frame()

for (i in 1:nrow(axon_samples)) {
  
  # selecting appropriate sample
  plot_sample <- axon_samples %>%
    slice(i)
  plot_df <- axon_original %>%
    filter(sample %in% plot_sample$sample)
  
  # processing distal
  plot_distal <- plot_df %>%
    filter(comp == "distal") %>%
    arrange(box)
  plot_distal_first <- plot_distal %>%
    slice(1)
  
  # only improving if first row is NaN
  if (is.na(plot_distal_first$mean) == T) {
    
    # identifying all NaN positions
    na_indices <- which(is.na(plot_distal$mean))
    initial_na_count <- 0
    
    # identifying initial NaN positions
    for (j in seq_along(na_indices)) {
      if (na_indices[j] == j) {
        initial_na_count <- initial_na_count + 1
      } else { break }
    }
    
    # removing initial NaN positions
    distal_remove <- na_indices[1:initial_na_count]
    distal_final <- plot_distal[-distal_remove, ]
    
  } else { distal_final <- plot_distal }
  
  
  # processing proximal
  plot_proximal <- plot_df %>%
    filter(comp == "proximal") %>%
    arrange(-box)
  plot_proximal_first <- plot_proximal %>%
    slice(1)
  
  # only improving if first row is NaN
  if (is.na(plot_proximal_first$mean) == T) {
    
    # identifying all NaN positions
    na_indices <- which(is.na(plot_proximal$mean))
    initial_na_count <- 0
    
    # identifying initial NaN positions
    for (j in seq_along(na_indices)) {
      if (na_indices[j] == j) {
        initial_na_count <- initial_na_count + 1
      } else { break }
    }
    
    # removing initial NaN positions
    proximal_remove <- na_indices[1:initial_na_count]
    proximal_final <- plot_proximal[-proximal_remove, ]
    
  } else { proximal_final <- plot_proximal }
  
  # processing somata
  plot_somata <- plot_df %>%
    filter(comp == "somata") %>%
    arrange(box)
  
  # processing somata layers
  plot_somataLayer <- plot_df %>%
    filter(comp == "somataLayer") %>%
    arrange(box)
  
  # final dataframe
  axon_betterROIs <- bind_rows(axon_betterROIs,
                               distal_final, 
                               proximal_final,
                               plot_somata,
                               plot_somataLayer)
}
  
# identifying proximal midline values
axon_maxProximal <- axon_betterROIs %>%
  filter(comp == "proximal") %>%
  group_by(sample) %>%
  summarize(proximal_pixels = max(distance_pixels),
            proximal_microns = max(distance_microns))
axon <- axon_betterROIs %>%
  left_join(axon_maxProximal, by = "sample") %>%
  mutate(full_comp = case_when(comp == "proximal" | 
                                    comp == "distal" ~ "axon",
                                  T ~ comp),
         full_box = case_when(comp == "proximal" ~ box - 500, 
                              T ~ box),
         full_pixels = case_when(comp == "proximal" ~ distance_pixels -
                                   proximal_pixels, 
                                 T ~ distance_pixels),
         full_microns = case_when(comp == "proximal" ~ distance_microns -
                                    proximal_microns, 
                                  T ~ distance_microns),
         thresholdFraction = as.numeric(thresholdFraction),
         thresholdNumber = as.numeric(thresholdNumber)) %>%
  replace(is.na(.), 0)

# plotting number of brains per construct and age
axon_metadata <- axon %>%
  distinct(exp, construct, brain, age) %>%
  group_by(construct, age) %>%
  summarize(number_brains = n())
ggplot(axon_metadata, aes(x = construct, y = number_brains)) +
  geom_col(pos = "dodge", fill = "black", color = "black") +
  facet_wrap(~ age) +
  theme_classic(base_size = 15, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black", 
                                   angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.title.x = element_blank())


## ----p3AxonExtendNoBEAM------------------------------------------------------------------------------------------------------------------
# filtering for p3 non-BEAM brains
axon_p3_noBEAM <- axon %>%
  filter(age == "p3") %>%
  filter(construct == "tdTomato" |
           construct == "Sema3f" |
           construct == "Nrp2" |
           construct == "PlexA3")

# calculating measurement sums for each ipsi section
axon_iue_p3_noBEAM_sums <- axon_p3_noBEAM %>%
  
  # restrict to less than 4000 microns
  filter(distance_microns < 4000) %>%
  filter(comp == "somata") %>%
  group_by(sample) %>%
  summarize(total_mean = sum(mean),
            total_totalIntensity = sum(total_intensity),
            total_threshFraction = sum(thresholdFraction),
            total_threshNumber = sum(thresholdNumber))

# calculating measurement percentages of each ipsi box per section
axon_iue_p3_noBEAM <- axon_p3_noBEAM %>%
  
  # restrict to less than 4000 microns
  filter(distance_microns < 4000) %>%
  filter(comp == "somata") %>%
  left_join(axon_iue_p3_noBEAM_sums, by = "sample") %>%
  mutate(percent_mean = mean / total_mean * 100,
         percent_totalIntensity = total_intensity / total_totalIntensity * 100,
         percent_fraction = thresholdFraction / total_threshFraction * 100,
         percent_number = thresholdNumber / total_threshNumber * 100)

# plotting all ipsi brains
ggplot(axon_iue_p3_noBEAM, 
       aes(x = distance_microns, y = percent_fraction)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = combo, color = combo, fill = combo)) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  labs(x = "Distance from midline (microns)",
       y = "Fraction of Pixels Above Threshold (%)") +
  facet_wrap(~ construct, scales = "free")

# calculating measurement sums for each contra section
axon_axon_p3_noBEAM_sums <- axon_p3_noBEAM %>%
  
  # restrict to less than 4000 microns and greater than -2000 microns
  filter(full_microns < 4000 & full_microns > -2000) %>%
  filter(full_comp == "axon") %>%
  group_by(sample) %>%
  summarize(total_mean = sum(mean),
            total_totalIntensity = sum(total_intensity),
            total_threshFraction = sum(thresholdFraction),
            total_threshNumber = sum(thresholdNumber))

# calculating measurement percentages of each contra box per section
axon_axon_p3_noBEAM <- axon_p3_noBEAM %>%
  
  # restrict to less than 4000 microns and greater than -2000 microns
  filter(full_microns < 4000 & full_microns > -2000) %>%
  filter(full_comp == "axon") %>%
  left_join(axon_axon_p3_noBEAM_sums, by = "sample") %>%
  mutate(percent_mean = mean / total_mean * 100,
         percent_totalIntensity = total_intensity / total_totalIntensity * 100,
         percent_fraction = thresholdFraction / total_threshFraction * 100,
         percent_number = thresholdNumber / total_threshNumber * 100)

# plotting all contra brains
ggplot(axon_axon_p3_noBEAM, 
       aes(x = full_microns, y = percent_fraction)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = combo, color = combo, fill = combo)) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  labs(x = "Distance from midline (microns)",
       y = "Fraction of Pixels Above Threshold (%)") +
  facet_wrap(~ construct, scales = "free")

# distinct axon samples
axon_p3_noBEAM_constructs <- axon_iue_p3_noBEAM %>%
  distinct(construct)

# plotting all constructs to identify bilateral IUEs
for (i in 1:nrow(axon_p3_noBEAM_constructs)) {

  # select appropriate construct
  plot_construct <- axon_p3_noBEAM_constructs %>%
    slice(i)

  # filter data for appropriate construct
  plot_axon_p3_noBEAM <- axon_iue_p3_noBEAM %>%
    filter(construct %in% plot_construct$construct)

  # plot cortical extension
  print(ggplot(plot_axon_p3_noBEAM, aes(x = distance_microns, y = percent_fraction)) +
    stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",

                # can play with span setting for different smoothness
                span = 0.1,
                aes(group = combo, color = combo, fill = combo)) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black")) +
    labs(x = "Distance from midline (microns)",
         y = "Fraction of Pixels Above Threshold (%)"))
}

# removing bilateral iues
axon_axon_p3_noBEAM <- axon_axon_p3_noBEAM
axon_iue_p3_noBEAM <- axon_iue_p3_noBEAM

# plotting all constructs with tdTomato average to identify medial IUEs
for (i in 1:nrow(axon_p3_noBEAM_constructs)) {
  
  # select appropriate construct
  plot_construct <- axon_p3_noBEAM_constructs %>%
    slice(i)
  
  # filter data for appropriate construct
  plot_axon_p3_noBEAM <- axon_iue_p3_noBEAM %>%
    filter(construct %in% plot_construct$construct)
  
  # plot electoroporation site
  print(ggplot(plot_axon_p3_noBEAM, aes(x = distance_microns, y = percent_fraction)) +
          stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
                      
                      # can play with span setting for different smoothness
                      span = 0.1,
                      aes(group = combo, color = combo, fill = combo)) +
          stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
                      
                      # plotting tdTomato data
                      data = axon_iue_p3_noBEAM %>%
                        filter(construct == "tdTomato") %>%  
                        transform(combo = NULL),
                      
                      # can play with span setting for different smoothness
                      span = 0.1,
                      aes(group = construct, color = construct, fill = construct)) +
          facet_wrap(~ combo, scales = "free") +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black"),
                axis.text.y = element_text(color = "black")) +
          labs(x = "Distance from midline (microns)",
               y = "Fraction of Pixels Above Threshold (%)"))
}

# identifying brains without medial iues
axon_iue_p3_noBEAM_medial <- axon_iue_p3_noBEAM %>%
  filter(combo != "186_tdTomato_p3_brain6") %>%
  filter(combo != "186_Sema3f_p3_brain1") %>%
  filter(combo != "186_Sema3f_p3_brain6") %>%
  filter(combo != "182_Nrp2_p3_brain1")
axon_axon_p3_noBEAM_medial <- axon_axon_p3_noBEAM %>%
  filter(combo %in% axon_iue_p3_noBEAM_medial$combo)

# plotting iues with medial electroporations
sema3f_p3_noBEAM_iue <- ggplot(axon_iue_p3_noBEAM_medial,
                        aes(x = distance_microns, y = percent_fraction)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",

              linewidth = 0.25,

              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = construct, color = construct, fill = construct)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#2B93CE", "#80CFB9", "#004F3A")) +
  scale_color_manual(values = c("#D55E00", "#2B93CE", "#80CFB9", "#004F3A")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.8)) +
  labs(x = "Electroporation: um From Midline",
       y = "% Thresholded",
       fill = "", color = ""); sema3f_p3_noBEAM_iue
ggsave("supp_fig_4b_ii.pdf", sema3f_p3_noBEAM_iue,
       units = "in", height = 1.00, width = 2.94)

# plotting iues with medial electroporations
sema3f_p3_noBEAM_axon <- ggplot(axon_axon_p3_noBEAM_medial,
                        aes(x = full_microns, y = percent_fraction)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",

              linewidth = 0.25,

              # can play with span setting for different smoothness
              span = 0.05,
              aes(group = construct, color = construct, fill = construct)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#2B93CE", "#80CFB9", "#004F3A")) +
  scale_color_manual(values = c("#D55E00", "#2B93CE", "#80CFB9", "#004F3A")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.8)) +
  labs(x = "Extension: um From Midline",
       y = "% Thresholded",
       fill = "", color = ""); sema3f_p3_noBEAM_axon
ggsave("supp_fig_4b_i.pdf", sema3f_p3_noBEAM_axon,
       units = "in", height = 1.08, width = 2.94)

test_means <- axon_axon_p3_noBEAM_medial %>%
  group_by(construct, full_box) %>%
  summarize(mean_percent = mean(percent_number, na.rm = T))

test_rep <- axon_axon_p3_noBEAM_medial %>%
  left_join(test_means, by = c("construct", "full_box")) %>%
  mutate(mean_diff = abs(mean_percent - percent_number)) %>%
  group_by(sample) %>%
  summarize(sum_diff = sum(mean_diff, na.rm = T))


## ----p0AxonExtendBEAM--------------------------------------------------------------------------------------------------------------------
# filtering for p0 BEAM brains
axon_p0_BEAM <- axon %>%
  filter(age == "p0") %>%
  filter(construct == "Sema3fB" |
           construct == "Nrp2B" |
           construct == "PlexA3B")

# calculating measurement sums for each ipsi section
axon_iue_p0_BEAM_sums <- axon_p0_BEAM %>%
  
  # restrict to less than 4000 microns
  filter(distance_microns < 4000) %>%
  filter(comp == "somata") %>%
  group_by(sample, color) %>%
  summarize(total_mean = sum(mean),
            total_totalIntensity = sum(total_intensity),
            total_threshFraction = sum(thresholdFraction),
            total_threshNumber = sum(thresholdNumber))

# calculating measurement percentages of each ipsi box per section
axon_iue_p0_BEAM <- axon_p0_BEAM %>%
  
  # restrict to less than 4000 microns
  filter(distance_microns < 4000) %>%
  filter(comp == "somata") %>%
  left_join(axon_iue_p0_BEAM_sums, by = c("sample", "color")) %>%
  mutate(percent_mean = mean / total_mean * 100,
         percent_totalIntensity = total_intensity / total_totalIntensity * 100,
         percent_fraction = thresholdFraction / total_threshFraction * 100,
         percent_number = thresholdNumber / total_threshNumber * 100)

# plotting all ipsi brains
ggplot(axon_iue_p0_BEAM, 
       aes(x = distance_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = combo, color = combo, fill = combo)) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  labs(x = "Distance from midline (microns)",
       y = "Mean Intensity (%)") +
  facet_wrap(~ construct + color, scales = "free")

# calculating measurement sums for each contra section
axon_axon_p0_BEAM_sums <- axon_p0_BEAM %>%
  
  # restrict to less than 4000 microns and greater than -2000 microns
  filter(full_microns < 4000 & full_microns > -2000) %>%
  filter(full_comp == "axon") %>%
  group_by(sample, color) %>%
  summarize(total_mean = sum(mean),
            total_totalIntensity = sum(total_intensity),
            total_threshFraction = sum(thresholdFraction),
            total_threshNumber = sum(thresholdNumber))

# calculating measurement percentages of each contra box per section
axon_axon_p0_BEAM <- axon_p0_BEAM %>%
  
  # restrict to less than 4000 microns and greater than -2000 microns
  filter(full_microns < 4000 & full_microns > -2000) %>%
  filter(full_comp == "axon") %>%
  left_join(axon_axon_p0_BEAM_sums, by = c("sample", "color")) %>%
  mutate(percent_mean = mean / total_mean * 100,
         percent_totalIntensity = total_intensity / total_totalIntensity * 100,
         percent_fraction = thresholdFraction / total_threshFraction * 100,
         percent_number = thresholdNumber / total_threshNumber * 100)

# plotting all contra brains
ggplot(axon_axon_p0_BEAM, 
       aes(x = full_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = combo, color = combo, fill = combo)) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  labs(x = "Distance from midline (microns)",
       y = "Mean Intensity (%)") +
  facet_wrap(~ construct + color, scales = "free")

# distinct axon samples
axon_p0_BEAM_constructs <- axon_iue_p0_BEAM %>%
  distinct(construct)

# plotting all constructs to identify bilateral IUEs
for (i in 1:nrow(axon_p0_BEAM_constructs)) {

  # select appropriate construct
  plot_construct <- axon_p0_BEAM_constructs %>%
    slice(i)

  # filter data for appropriate construct
  plot_axon_p0_BEAM <- axon_iue_p0_BEAM %>%
    filter(construct %in% plot_construct$construct)

  # plot cortical extension
  print(ggplot(plot_axon_p0_BEAM, aes(x = distance_microns, y = percent_mean)) +
    stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",

                # can play with span setting for different smoothness
                span = 0.1,
                aes(group = combo, color = combo, fill = combo)) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black")) +
    labs(x = "Distance from midline (microns)",
         y = "Mean Intensity (%)") +
      facet_wrap(~ color, scales = "free"))
}

# removing bilateral iues
axon_axon_p0_BEAM <- axon_axon_p0_BEAM %>%
  filter(combo != "204_PlexA3B_p0_brain1") %>%
  filter(combo != "204_PlexA3B_p0_brain2")
axon_iue_p0_BEAM <- axon_iue_p0_BEAM %>%
  filter(combo != "204_PlexA3B_p0_brain1") %>%
  filter(combo != "204_PlexA3B_p0_brain2")

# plotting all constructs with tdTomato average to identify medial IUEs
for (i in 1:nrow(axon_p0_BEAM_constructs)) {
  
  # select appropriate construct
  plot_construct <- axon_p0_BEAM_constructs %>%
    slice(i)
  
  # filter data for appropriate construct
  plot_axon_p0_BEAM <- axon_iue_p0_BEAM %>%
    filter(construct %in% plot_construct$construct)
  
  # plot electoroporation site
  print(ggplot(plot_axon_p0_BEAM, aes(x = distance_microns, y = percent_mean)) +
          stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
                      
                      # can play with span setting for different smoothness
                      span = 0.1,
                      aes(group = combo, color = combo, fill = combo)) +
          stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
                      
                      # plotting average for all constructs
                      data = axon_iue_p0_BEAM %>%  
                        transform(combo = NULL),
                      
                      # can play with span setting for different smoothness
                      span = 0.1,
                      aes(group = color, color = color, fill = color)) +
          facet_wrap(~ combo + color, scales = "free") +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black"),
                axis.text.y = element_text(color = "black")) +
          labs(x = "Distance from midline (microns)",
               y = "Mean Intensity (%)"))
}

# identifying brains without medial iues
axon_iue_p0_BEAM_medial <- axon_iue_p0_BEAM %>%
  filter(combo != "204_Nrp2B_p0_brain1") %>%
  filter(combo != "209_PlexA3B_p0_brain2") %>%
  filter(combo != "209_PlexA3B_p0_brain3") %>%
  filter(combo != "209_PlexA3B_p0_brain4") %>%
  filter(combo != "209_PlexA3B_p0_brain5") %>%
  filter(combo != "206_Sema3fB_p0_brain1") %>%
  mutate(Condition = case_when(color == "green" ~ "WT",
                            color == "red" ~ "OE"),
         Condition = fct_relevel(Condition, "OE", "WT"))
axon_axon_p0_BEAM_medial <- axon_axon_p0_BEAM %>%
  filter(combo %in% axon_iue_p0_BEAM_medial$combo) %>%
  mutate(Condition = case_when(color == "green" ~ "WT",
                            color == "red" ~ "OE"),
         Condition = fct_relevel(Condition, "OE", "WT"))

# plotting iues with medial electroporations
sema3f_p0_BEAM_iue <- ggplot(axon_iue_p0_BEAM_medial %>%
         filter(construct == "Sema3fB"), aes(x = distance_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = Condition, color = Condition, fill = Condition)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)) +
  labs(x = "um From Midline",
       y = "% Intensity",
       color = "", fill = ""); sema3f_p0_BEAM_iue
ggsave("supp_fig_6a_i.pdf", sema3f_p0_BEAM_iue, units = "in",
       width = 1.44, height = 1.27)

nrp2_p0_BEAM_iue <- ggplot(axon_iue_p0_BEAM_medial %>%
         filter(construct == "Nrp2B"), aes(x = distance_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = Condition, color = Condition, fill = Condition)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)) +
  labs(x = "um From Midline",
       y = "% Intensity",
       color = "", fill = ""); nrp2_p0_BEAM_iue
ggsave("supp_fig_6b_i.pdf", nrp2_p0_BEAM_iue, units = "in",
       width = 1.44, height = 1.09)

plexa3_p0_BEAM_iue <- ggplot(axon_iue_p0_BEAM_medial %>%
         filter(construct == "PlexA3B"), aes(x = distance_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = Condition, color = Condition, fill = Condition)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)) +
  labs(x = "um From Midline",
       y = "% Intensity",
       color = "", fill = ""); plexa3_p0_BEAM_iue
ggsave("supp_fig_6c_i.pdf", plexa3_p0_BEAM_iue, units = "in",
       width = 1.44, height = 1.14)

# plotting iues with medial electroporations
ggplot(axon_axon_p0_BEAM_medial, aes(x = full_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              # can play with span setting for different smoothness
              span = 0.05,
              aes(group = Condition, color = Condition, fill = Condition)) +
  facet_wrap(~ construct, scales = "free", nrow = 3) +
  theme_classic() +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Distance from midline (um)",
       y = "% Mean Intensity",
       color = "", fill = "")

sema3f_p0_BEAM_axon <- ggplot(axon_axon_p0_BEAM_medial %>%
         filter(construct == "Sema3fB"), aes(x = full_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              
              # can play with span setting for different smoothness
              span = 0.05,
              aes(group = Condition, color = Condition, fill = Condition)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)) +
  labs(x = "um From Midline",
       y = "% Intensity",
       color = "", fill = ""); sema3f_p0_BEAM_axon
ggsave("supp_fig_5a.pdf", sema3f_p0_BEAM_axon, units = "in",
       width = 1.47, height = 0.89)

nrp2_p0_BEAM_axon <- ggplot(axon_axon_p0_BEAM_medial %>%
         filter(construct == "Nrp2B"), aes(x = full_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              
              # can play with span setting for different smoothness
              span = 0.05,
              aes(group = Condition, color = Condition, fill = Condition)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)) +
  labs(x = "um From Midline",
       y = "% Intensity",
       color = "", fill = ""); nrp2_p0_BEAM_axon
ggsave("supp_fig_5b.pdf", nrp2_p0_BEAM_axon, units = "in",
       width = 1.47, height = 0.95)

plexa3_p0_BEAM_axon <- ggplot(axon_axon_p0_BEAM_medial %>%
         filter(construct == "PlexA3B"), aes(x = full_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              
              # can play with span setting for different smoothness
              span = 0.05,
              aes(group = Condition, color = Condition, fill = Condition)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)) +
  labs(x = "um From Midline",
       y = "% Intensity",
       color = "", fill = ""); plexa3_p0_BEAM_axon
ggsave("supp_fig_5c.pdf", plexa3_p0_BEAM_axon, units = "in",
       width = 1.47, height = 1.12)

# finding representative section
test_means <- axon_axon_p0_BEAM_medial %>%
  group_by(construct, color, full_box) %>%
  summarize(mean_percent = mean(percent_mean, na.rm = T))

test_rep <- axon_axon_p0_BEAM_medial %>%
  left_join(test_means, by = c("construct", "color", "full_box")) %>%
  mutate(mean_diff = abs(mean_percent - percent_mean)) %>%
  group_by(sample) %>%
  summarize(sum_diff = sum(mean_diff, na.rm = T))

# calculating measurement sums for each layer section
axon_layer_p0_BEAM_sums <- axon_p0_BEAM %>%
  
  filter(comp == "somataLayer") %>%
  group_by(sample, color) %>%
  summarize(total_mean = sum(mean),
            total_totalIntensity = sum(total_intensity),
            total_threshFraction = sum(thresholdFraction),
            total_threshNumber = sum(thresholdNumber),
            max_distance = max(distance_microns))

# calculating measurement percentages of each layer box per section
axon_layer_p0_BEAM <- axon_p0_BEAM %>%
  
  filter(comp == "somataLayer") %>%
  left_join(axon_layer_p0_BEAM_sums, by = c("sample", "color")) %>%
  mutate(percent_mean = mean / total_mean * 100,
         percent_totalIntensity = total_intensity / total_totalIntensity * 100,
         percent_fraction = thresholdFraction / total_threshFraction * 100,
         percent_number = thresholdNumber / total_threshNumber * 100,
         surface_box = -(box-499),
         surface_microns = -(distance_microns - max_distance))

# plotting all layer brains
ggplot(axon_layer_p0_BEAM, 
       aes(x = surface_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = combo, color = combo, fill = combo)) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  labs(x = "Distance from surface (um)",
       y = "Mean Intensity (%)") +
  facet_wrap(~ construct + color, scales = "free")

axon_layer_p0_BEAM_medial <- axon_layer_p0_BEAM %>%
  filter(combo %in% axon_iue_p0_BEAM_medial$combo) %>%
  mutate(Condition = case_when(color == "green" ~ "WT",
                            color == "red" ~ "OE"),
         Condition = fct_relevel(Condition, "OE", "WT"))

# plotting iues with medial electroporations
ggplot(axon_layer_p0_BEAM_medial, aes(x = surface_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = Condition, color = Condition, fill = Condition)) +
  facet_wrap(~ construct, scales = "free", nrow = 3) +
  theme_classic() +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Distance from surface (um)",
       y = "% Mean Intensity",
       color = "", fill = "")

# plotting iues with medial electroporations
sema3f_p0_BEAM_soma <- ggplot(axon_layer_p0_BEAM_medial %>%
         filter(construct == "Sema3fB"), 
         aes(x = surface_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = Condition, color = Condition, fill = Condition)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)) +
  labs(x = "um From Surface",
       y = "% Intensity",
       color = "", fill = ""); sema3f_p0_BEAM_soma
ggsave("supp_fig_6a_ii.pdf", sema3f_p0_BEAM_soma, units = "in",
       width = 1.44, height = 1.27)

nrp2_p0_BEAM_soma <- ggplot(axon_layer_p0_BEAM_medial %>%
         filter(construct == "Nrp2B"), 
         aes(x = surface_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = Condition, color = Condition, fill = Condition)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)) +
  labs(x = "um From Surface",
       y = "% Intensity",
       color = "", fill = ""); nrp2_p0_BEAM_soma
ggsave("supp_fig_6b_ii.pdf", nrp2_p0_BEAM_soma, units = "in",
       width = 1.44, height = 1.09)

plexa3_p0_BEAM_soma <- ggplot(axon_layer_p0_BEAM_medial %>%
         filter(construct == "PlexA3B"), 
         aes(x = surface_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = Condition, color = Condition, fill = Condition)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)) +
  labs(x = "um From Surface",
       y = "% Intensity",
       color = "", fill = ""); plexa3_p0_BEAM_soma
ggsave("supp_fig_6c_ii.pdf", plexa3_p0_BEAM_soma, units = "in",
       width = 1.44, height = 1.14)


## ----p3AxonExtendBEAM--------------------------------------------------------------------------------------------------------------------
# filtering for p3 BEAM brains
axon_p3_BEAM <- axon %>%
  filter(age == "p3") %>%
  filter(construct == "Sema3fB" |
           construct == "Nrp2B" |
           construct == "PlexA3B")

# calculating measurement sums for each ipsi section
axon_iue_p3_BEAM_sums <- axon_p3_BEAM %>%
  
  # restrict to less than 4000 microns
  filter(distance_microns < 4000) %>%
  filter(comp == "somata") %>%
  group_by(sample, color) %>%
  summarize(total_mean = sum(mean),
            total_totalIntensity = sum(total_intensity),
            total_threshFraction = sum(thresholdFraction),
            total_threshNumber = sum(thresholdNumber))

# calculating measurement percentages of each ipsi box per section
axon_iue_p3_BEAM <- axon_p3_BEAM %>%
  
  # restrict to less than 4000 microns
  filter(distance_microns < 4000) %>%
  filter(comp == "somata") %>%
  left_join(axon_iue_p3_BEAM_sums, by = c("sample", "color")) %>%
  mutate(percent_mean = mean / total_mean * 100,
         percent_totalIntensity = total_intensity / total_totalIntensity * 100,
         percent_fraction = thresholdFraction / total_threshFraction * 100,
         percent_number = thresholdNumber / total_threshNumber * 100)

# plotting all ipsi brains
ggplot(axon_iue_p3_BEAM, 
       aes(x = distance_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = combo, color = combo, fill = combo)) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  labs(x = "Distance from midline (microns)",
       y = "Mean Intensity (%)") +
  facet_wrap(~ construct + color, scales = "free")

# calculating measurement sums for each contra section
axon_axon_p3_BEAM_sums <- axon_p3_BEAM %>%
  
  # restrict to less than 4000 microns and greater than -2000 microns
  filter(full_microns < 4000 & full_microns > -2000) %>%
  filter(full_comp == "axon") %>%
  group_by(sample, color) %>%
  summarize(total_mean = sum(mean),
            total_totalIntensity = sum(total_intensity),
            total_threshFraction = sum(thresholdFraction),
            total_threshNumber = sum(thresholdNumber))

# calculating measurement percentages of each contra box per section
axon_axon_p3_BEAM <- axon_p3_BEAM %>%
  
  # restrict to less than 4000 microns and greater than -2000 microns
  filter(full_microns < 4000 & full_microns > -2000) %>%
  filter(full_comp == "axon") %>%
  left_join(axon_axon_p3_BEAM_sums, by = c("sample", "color")) %>%
  mutate(percent_mean = mean / total_mean * 100,
         percent_totalIntensity = total_intensity / total_totalIntensity * 100,
         percent_fraction = thresholdFraction / total_threshFraction * 100,
         percent_number = thresholdNumber / total_threshNumber * 100)

# plotting all contra brains
ggplot(axon_axon_p3_BEAM, 
       aes(x = full_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = combo, color = combo, fill = combo)) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  labs(x = "Distance from midline (microns)",
       y = "Mean Intensity (%)") +
  facet_wrap(~ construct + color, scales = "free")

# distinct axon samples
axon_p3_BEAM_constructs <- axon_iue_p3_BEAM %>%
  distinct(construct)

# plotting all constructs to identify bilateral IUEs
for (i in 1:nrow(axon_p3_BEAM_constructs)) {

  # select appropriate construct
  plot_construct <- axon_p3_BEAM_constructs %>%
    slice(i)

  # filter data for appropriate construct
  plot_axon_p3_BEAM <- axon_iue_p3_BEAM %>%
    filter(construct %in% plot_construct$construct)

  # plot cortical extension
  print(ggplot(plot_axon_p3_BEAM, aes(x = distance_microns, y = percent_mean)) +
    stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",

                # can play with span setting for different smoothness
                span = 0.1,
                aes(group = combo, color = combo, fill = combo)) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black")) +
    labs(x = "Distance from midline (microns)",
         y = "Mean Intensity (%)") +
      facet_wrap(~ color, scales = "free"))
}

# removing bilateral iues
axon_axon_p3_BEAM <- axon_axon_p3_BEAM
axon_iue_p3_BEAM <- axon_iue_p3_BEAM

# plotting all constructs with tdTomato average to identify medial IUEs
for (i in 1:nrow(axon_p3_BEAM_constructs)) {
  
  # select appropriate construct
  plot_construct <- axon_p3_BEAM_constructs %>%
    slice(i)
  
  # filter data for appropriate construct
  plot_axon_p3_BEAM <- axon_iue_p3_BEAM %>%
    filter(construct %in% plot_construct$construct)
  
  # plot electoroporation site
  print(ggplot(plot_axon_p3_BEAM, aes(x = distance_microns, y = percent_mean)) +
          stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
                      
                      # can play with span setting for different smoothness
                      span = 0.1,
                      aes(group = combo, color = combo, fill = combo)) +
          stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
                      
                      # plotting average for all constructs
                      data = axon_iue_p3_BEAM %>%  
                        transform(combo = NULL),
                      
                      # can play with span setting for different smoothness
                      span = 0.1,
                      aes(group = color, color = color, fill = color)) +
          facet_wrap(~ combo + color, scales = "free") +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black"),
                axis.text.y = element_text(color = "black")) +
          labs(x = "Distance from midline (microns)",
               y = "Mean Intensity (%)"))
}

# identifying brains without medial iues
axon_iue_p3_BEAM_medial <- axon_iue_p3_BEAM %>%
  filter(combo != "210_PlexA3B_p3_brain2") %>%
  filter(combo != "212_Sema3fB_p3_brain3") %>%
  mutate(Condition = case_when(color == "green" ~ "WT",
                            color == "red" ~ "OE"),
         Condition = fct_relevel(Condition, "OE", "WT"))
axon_axon_p3_BEAM_medial <- axon_axon_p3_BEAM %>%
  filter(combo %in% axon_iue_p3_BEAM_medial$combo) %>%
  mutate(Condition = case_when(color == "green" ~ "WT",
                            color == "red" ~ "OE"),
         Condition = fct_relevel(Condition, "OE", "WT"))

# plotting iues with medial electroporations
ggplot(axon_iue_p3_BEAM_medial, aes(x = distance_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = Condition, color = Condition, fill = Condition)) +
  facet_wrap(~ construct, scales = "free", nrow = 3) +
  theme_classic() +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Distance from midline (um)",
       y = "% Mean Intensity",
       color = "", fill = "")

sema3f_p3_BEAM_iue <- ggplot(axon_iue_p3_BEAM_medial %>%
                            filter(construct == "Sema3fB"), 
                          aes(x = distance_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = Condition, color = Condition, fill = Condition)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)) +
  labs(x = "um From Midline",
       y = "% Intensity",
       color = "", fill = ""); sema3f_p3_BEAM_iue
ggsave("supp_fig_6d_i.pdf", sema3f_p3_BEAM_iue, units = "in",
       width = 1.47, height = 1.22)
nrp2_p3_BEAM_iue <- ggplot(axon_iue_p3_BEAM_medial %>%
                            filter(construct == "Nrp2B"), 
                          aes(x = distance_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = Condition, color = Condition, fill = Condition)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)) +
  labs(x = "um From Midline",
       y = "% Intensity",
       color = "", fill = ""); nrp2_p3_BEAM_iue
ggsave("supp_fig_6e_i.pdf", nrp2_p3_BEAM_iue, units = "in",
       width = 1.47, height = 1.29)
plexa3_p3_BEAM_iue <- ggplot(axon_iue_p3_BEAM_medial %>%
                            filter(construct == "PlexA3B"), 
                          aes(x = distance_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = Condition, color = Condition, fill = Condition)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)) +
  labs(x = "um From Midline",
       y = "% Intensity",
       color = "", fill = ""); plexa3_p3_BEAM_iue
ggsave("supp_fig_6f_i.pdf", plexa3_p3_BEAM_iue, units = "in",
       width = 1.47, height = 0.93)

# plotting iues with medial electroporations
ggplot(axon_axon_p3_BEAM_medial, aes(x = full_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              # can play with span setting for different smoothness
              span = 0.05,
              aes(group = Condition, color = Condition, fill = Condition)) +
  facet_wrap(~ construct, scales = "free", nrow = 3) +
  theme_classic() +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Distance from midline (um)",
       y = "% Mean Intensity",
       color = "", fill = "")

sema3f_p3_BEAM_axon <- ggplot(axon_axon_p3_BEAM_medial %>%
                            filter(construct == "Sema3fB"), 
                          aes(x = full_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              
              # can play with span setting for different smoothness
              span = 0.05,
              aes(group = Condition, color = Condition, fill = Condition)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)) +
  labs(x = "um From Midline",
       y = "% Intensity",
       color = "", fill = ""); sema3f_p3_BEAM_axon
ggsave("supp_fig_5d.pdf", sema3f_p3_BEAM_axon, units = "in",
       width = 1.47, height = 1.50)

nrp2_p3_BEAM_axon <- ggplot(axon_axon_p3_BEAM_medial %>%
                            filter(construct == "Nrp2B"), 
                          aes(x = full_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              
              # can play with span setting for different smoothness
              span = 0.05,
              aes(group = Condition, color = Condition, fill = Condition)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)) +
  labs(x = "um From Midline",
       y = "% Intensity",
       color = "", fill = ""); nrp2_p3_BEAM_axon
ggsave("supp_fig_5e.pdf", nrp2_p3_BEAM_axon, units = "in",
       width = 1.47, height = 1.61)

plexa3_p3_BEAM_axon <- ggplot(axon_axon_p3_BEAM_medial %>%
                            filter(construct == "PlexA3B"), 
                          aes(x = full_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              
              # can play with span setting for different smoothness
              span = 0.05,
              aes(group = Condition, color = Condition, fill = Condition)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)) +
  labs(x = "um From Midline",
       y = "% Intensity",
       color = "", fill = ""); plexa3_p3_BEAM_axon
ggsave("supp_fig_5f.pdf", plexa3_p3_BEAM_axon, units = "in",
       width = 1.47, height = 1.68)

# finding representative sections
test_means <- axon_axon_p3_BEAM_medial %>%
  group_by(construct, color, full_box) %>%
  summarize(mean_percent = mean(percent_mean, na.rm = T))

test_rep <- axon_axon_p3_BEAM_medial %>%
  left_join(test_means, by = c("construct", "color", "full_box")) %>%
  mutate(mean_diff = abs(mean_percent - percent_mean)) %>%
  group_by(sample) %>%
  summarize(sum_diff = sum(mean_diff, na.rm = T))

# calculating measurement sums for each layer section
axon_layer_p3_BEAM_sums <- axon_p3_BEAM %>%
  
  filter(comp == "somataLayer") %>%
  group_by(sample, color) %>%
  summarize(total_mean = sum(mean),
            total_totalIntensity = sum(total_intensity),
            total_threshFraction = sum(thresholdFraction),
            total_threshNumber = sum(thresholdNumber),
            max_distance = max(distance_microns))

# calculating measurement percentages of each layer box per section
axon_layer_p3_BEAM <- axon_p3_BEAM %>%
  
  filter(comp == "somataLayer") %>%
  left_join(axon_layer_p3_BEAM_sums, by = c("sample", "color")) %>%
  mutate(percent_mean = mean / total_mean * 100,
         percent_totalIntensity = total_intensity / total_totalIntensity * 100,
         percent_fraction = thresholdFraction / total_threshFraction * 100,
         percent_number = thresholdNumber / total_threshNumber * 100,
         surface_box = -(box-499),
         surface_microns = -(distance_microns - max_distance))

# plotting all layer brains
ggplot(axon_layer_p3_BEAM, 
       aes(x = surface_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = combo, color = combo, fill = combo)) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        legend.position = "none") +
  labs(x = "Distance from surface (um)",
       y = "Mean Intensity (%)") +
  facet_wrap(~ construct + color, scales = "free")

axon_layer_p3_BEAM_medial <- axon_layer_p3_BEAM %>%
  filter(combo %in% axon_iue_p3_BEAM_medial$combo) %>%
  mutate(Condition = case_when(color == "green" ~ "WT",
                            color == "red" ~ "OE"),
         Condition = fct_relevel(Condition, "OE", "WT"))

# plotting iues with medial electroporations
ggplot(axon_layer_p3_BEAM_medial, aes(x = surface_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = Condition, color = Condition, fill = Condition)) +
  facet_wrap(~ construct, scales = "free", nrow = 3) +
  theme_classic() +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Distance from surface (um)",
       y = "% Mean Intensity",
       color = "", fill = "")

sema3f_p3_BEAM_soma <- ggplot(axon_layer_p3_BEAM_medial %>%
                            filter(construct == "Sema3fB"), 
                          aes(x = surface_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = Condition, color = Condition, fill = Condition)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)) +
  labs(x = "um From Surface",
       y = "% Intensity",
       color = "", fill = ""); sema3f_p3_BEAM_soma
ggsave("supp_fig_6d_ii.pdf", sema3f_p3_BEAM_soma, units = "in",
       width = 1.47, height = 1.22)
nrp2_p3_BEAM_soma <- ggplot(axon_layer_p3_BEAM_medial %>%
                            filter(construct == "Nrp2B"), 
                          aes(x = surface_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = Condition, color = Condition, fill = Condition)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)) +
  labs(x = "um From Surface",
       y = "% Intensity",
       color = "", fill = ""); nrp2_p3_BEAM_soma
ggsave("supp_fig_6e_ii.pdf", nrp2_p3_BEAM_soma, units = "in",
       width = 1.47, height = 1.29)
plexa3_p3_BEAM_soma <- ggplot(axon_layer_p3_BEAM_medial %>%
                            filter(construct == "PlexA3B"), 
                          aes(x = surface_microns, y = percent_mean)) +
  stat_smooth(se = T, level = 0.95, method = "loess", formula = "y ~ x",
              
              linewidth = 0.25,
              
              # can play with span setting for different smoothness
              span = 0.1,
              aes(group = Condition, color = Condition, fill = Condition)) +
  theme_classic(base_size = 7) +
  scale_fill_manual(values = c("#D55E00", "#009E73")) +
  scale_color_manual(values = c("#D55E00", "#009E73")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.9)) +
  labs(x = "um From Surface",
       y = "% Intensity",
       color = "", fill = ""); plexa3_p3_BEAM_soma
ggsave("supp_fig_6f_ii.pdf", plexa3_p3_BEAM_soma, units = "in",
       width = 1.47, height = 0.93)


## ----sessionInfo-------------------------------------------------------------------------------------------------------------------------
sessionInfo()

