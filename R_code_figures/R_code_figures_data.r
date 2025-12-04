#script used to generate figures for the publication Soil microbial diversity across tropical savanna ecosystems"
#responsible code author: Lea Nosalova
#this script is the first version, used to analyse and present the data


#to configurate the working space
# INIT ----
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #always set current dir as WD
source("0_configuration.R")


#to create figures and supplementary figures
###WHITTAKER BIOME PLOT
#https://rdrr.io/github/valentinitnelav/plotbiomes/f/html/Whittaker_biomes_examples.Rmd
library(plotbiomes)
library(ggplot2)
library(raster)
library(sp)

whittaker_base_plot() + theme_bw()

#low resolution of plotbiomes dataset
path <- system.file("extdata", "temp_pp.tif", package = "plotbiomes")
temp_pp <- raster::stack(path)
names(temp_pp) <- c("temperature", "precipitation")

#add locations
my_sites <- data.frame(site = c("Burkina Faso", "Colombia", "Ivory Coast", "Namibia","South Africa site 1 (Skukuza)", "South Africa site 2 (Satara)"),
                       lon = c(-4.340744, -71.337413, -5.0217472, 18.447108, 31.766753, 31.778191),
                       lat = c(11.087943, 4.572624, 6.2312639, -23.220712, -24.966832, -24.393179))

coordinates(my_sites) <- ~lon + lat
proj4string(my_sites) <- CRS("+proj=longlat +datum=WGS84")
path <- system.file("extdata", "temp_pp.tif", package = "plotbiomes")
temp_pp <- raster::stack(path)
names(temp_pp) <- c("temperature", "precipitation")

extractions <- raster::extract(temp_pp, my_sites, df = TRUE)

#temperature is stored x10, precipitation is mm convert to cm
extractions$temperature <- extractions$temperature / 10
extractions$precipitation <- extractions$precipitation / 10

extractions$site <- my_sites$site


plot(temp_pp[[1]] / 10, main = "Mean Annual Temperature (°C)")
points(my_sites, pch = 21, bg = "red", cex = 1.5)

plot(temp_pp[[2]] / 10, main = "Mean Annual Precipitation (cm)")
points(my_sites, pch = 21, bg = "blue", cex = 1.5)

#plot whittaker biome plot with extracted sites
whittaker_base_plot() +
  # add the temperature - precipitation data points
  geom_point(data = extractions,
             aes(x = temperature,
                 y = precipitation,
             size = 3,
             shape = 21) +
  geom_point(data = extractions,
             aes(x = temperature,
                 y = precipitation,
                 color = site),
             shape = 16,
             size  = 3) +
  theme_bw() + 
  scale_colour_manual(values = cbPalette) 



###ALPHA DIVERSITY
#https://www.spsanderson.com/steveondata/posts/2024-09-24/#introduction
#https://rpubs.com/RosaneRech/OneFactorBoxplot
library(phyloseq)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(car)
library(multcompView)


cbPalette <- c("#CC79A7","#009E73", "#999999",  "#E69F00", "#56B4E9", "#0072B2",  "#F0E442",  "#D55E00" )
fileAlpha <- file.path(alpha_diversity_folder, "alpha_diversity.csv") #path to output Alpha diversity folder

set.seed(123)
ps = readRDS(ps_rarefied_file) #important to load ps object which was already filtered and rarefied
measures <-  c("Observed","Shannon","InvSimpson", "Chao1")
alpha_diversity <- estimate_richness(ps, measures = measures)
alpha_diversity <- cbind(sample_data(ps), alpha_diversity)
write_csv(alpha_diversity, fileAlpha)

#to filter out outliers by tukey method
my_data <- alpha_diversity %>% dplyr::select(Shannon, Observed, InvSimpson, Chao1, sampling_point)
my_data <- my_data %>% filter(sampling_point == "satara")
my_data <- my_data %>% dplyr::select(Shannon, Observed, InvSimpson, Chao1)

my_data_cleaned <- my_data 

for(sampling_point in unique(my_data$sampling_point)) {
  
  #loop through the indexes (excluding the sampling_point)
  for(col in names(my_data)[-length(names(my_data))]) {
    
    #subset data for the current sampling point directly in the loop
    current_data <- my_data[my_data$sampling_point == sampling_point, col]
    #calculate quantiles and IQR for the current column (current_data contains only data for the sampling point)
    Q1 <- quantile(current_data, 0.25, na.rm = TRUE)
    Q3 <- quantile(current_data, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    #identify outliers for this column and sampling point
    outliers <- current_data < (Q1 - 1.5 * IQR) | current_data > (Q3 + 1.5 * IQR)
    
    #replace outliers with NA in the cleaned data
    my_data_cleaned[my_data$sampling_point == sampling_point & my_data[[col]] %in% current_data[outliers], col] <- NA
  }
}

write.csv(my_data_cleaned, "my_data_cleaned.csv", row.names = FALSE)

summary(my_data)
summary(my_data_cleaned)


#check for the normality and assumptions for anova
anova  <- aov(Shannon ~ sampling_point, data = my_data_cleaned)
summary(anova)

tukey <- TukeyHSD(anova)
print(tukey)

cld <- multcompLetters4(anova, tukey)
print(cld)


Tk <- my_data_cleaned %>%
  group_by(sampling_point) %>%
  summarise(mean = mean(Shannon, na.rm = TRUE), 
            quant = quantile(Shannon, probs = 0.75, na.rm = TRUE)) %>%
  arrange(desc(mean))

#extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$sampling_point)
Tk$cld <- cld$Letters
print(Tk)


#homogenity and normality of residuals check
par(mfrow = c(3,2))
plot(anova)
hist(resid(anova), las=1)
qqnorm(resid(anova), las=1)  
qqline(resid(anova), col = "red")
hist(resid(anova), las=1)  
plot(fitted(anova), resid(anova), las=1)  
abline(h = 0, col = "red", lwd = 2)  
plot(my_data_cleaned$sampling_point, resid(anova), las=1)  

#normality - Shapiro on residuals
shapiro.test(residuals(anova))

#homogeneity of variation assumption (Bartlette/Levane)
bartlett.test(InvSimpson~sampling_point, data=my_data_cleaned)
leveneTest(Shannon~sampling_point, data=my_data_cleaned)
residuals_data <- residuals(anova)
bartlett.test(residuals(anova) ~ my_data_cleaned$sampling_point)
bptest(anova)


#plot indexes
grid.arrange(p1, p2, p3, nrow = 1)
p1 <- ggplot(my_data_cleaned, aes(x = sampling_point, y = Observed, fill = sampling_point)) + 
  geom_boxplot(coef = 1.5) + 
  labs(x="", y="Richness (Observed ASV)") +
  theme_bw() +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_text(data = Tk, aes(x = sampling_point, y = quant  + 170 , label = cld), 
            size = 3, vjust = -1, hjust = 0.5) +
  scale_fill_manual(values =  cbPalette) +
  scale_x_discrete(labels = NULL)
p2 <- ggplot(my_data_cleaned, aes(x = sampling_point, y = Shannon, fill = sampling_point)) + 
  geom_boxplot(coef = 1.5) + 
  labs(x="", y="Shannon index") +
  theme_bw() +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_text(data = Tk, aes(x = sampling_point, y = quant  + 0.35 , label = cld), 
            size = 3, vjust = -1, hjust = 0.1) +
  scale_fill_manual(values =  cbPalette) +
  scale_x_discrete(labels = NULL)
p3 <- ggplot(my_data_cleaned, aes(x = sampling_point, y = InvSimpson, fill = sampling_point)) + 
  geom_boxplot(coef = 1.5) + 
  labs(x="", y="Inverse Simpson index") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_text(data = Tk, aes(x = sampling_point, y = quant + 110, label = cld), 
            size = 3, vjust = -1, hjust = 0.5) +
  scale_fill_manual(
    values = cbPalette,
    labels = c("Burkina Faso", "Colombia", "Ivory Coast", "Namibia", "South Africa 1 (site Satara)", "South Africa 2 (site Skukuza)" )
  ) +
  scale_x_discrete(labels = NULL) +
  theme(legend.position = "right" )
  
  
###BOXPLOTS bare vs plant
###significance letters per site for each aloha diversity index - boxplots
library(multcompView)
library(dplyr)
library(graphics)

summary(my_data_cleaned)
par(mfrow = c(2,2))
boxplot(Shannon ~ sampling_point, data=my_data, main="Boxplot for Shannon")
boxplot(Observed ~ sampling_point, data=my_data, main="Boxplot for Observed")
boxplot(Chao1 ~ sampling_point, data=my_data,  main="Boxplot for Chao1")
boxplot(InvSimpson ~ sampling_point, data=my_data,  main="Boxplot for InvSimpson")

my_data_cleaned$source_types <- alpha_diversity$source_types
my_alpha_data  
my_alpha_data <- my_data_cleaned %>% filter(sampling_point == "namibia") %>%   dplyr::select(InvSimpson, sampling_point, source_types)

#check for the normality and assumptions for anova
anova  <- aov(InvSimpson ~ source_types, data = my_alpha_data)
summary(anova)

tukey <- TukeyHSD(anova)
print(tukey)

cld <- multcompLetters4(anova, tukey)
print(cld)


Tk <- my_alpha_data %>%
  group_by(sampling_point, source_types) %>%
  summarise(mean = mean(InvSimpson, na.rm = TRUE), 
            quant = quantile(InvSimpson, probs = 0.75, na.rm = TRUE)) %>%
  arrange(desc(mean))
#extracting the compact letter display and adding to the Tk table
cld <- as.data.frame.list(cld$source_types)
Tk$cld <- cld$Letters
print(Tk)

write_csv(Tk, "Tk.csv")

###I calculared it per site and then added it to the dataframe
Tk <- Tk

palett_group <- c("springgreen4", "#cc6666")

level_order_subgroup <- c("plant", "bare")
p1 <- ggplot(my_data_cleaned, aes(x = sampling_point, y = InvSimpson, fill = factor(source_types, level = level_order_subgroup ))) + 
  geom_boxplot(coef = 1.5) + 
  labs(x="", y="InvSimpson index") +
  theme_bw() +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_text(data = Tk, aes( x = sampling_point, y = quant  + 50 ,   group = source_types, label = cld),
            size = 3, vjust = -1, hjust = 0.5) +
  scale_fill_manual(values =  palett_group) +
  #scale_x_discrete(labels = NULL)
  scale_x_discrete(labels = c("Burkina Faso", "Colombia", "Ivory Coast", "Namibia", "South Africa 1 (site Satara)", "South Africa 2 (site Skukuza)"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.position = "right", legend.title = element_blank()  )
p1  
    

###BETA DIVERSITY
#adapted from Lebre et al. 2023 
#https://zenodo.org/records/7997130
library(phyloseq)
library(vegan)
library(ggplot2)
library(pals)

cbPalette <- c("#CC79A7","#009E73", "#999999",  "#E69F00", "#56B4E9", "#0072B2",  "#F0E442",  "#D55E00" )
palette_specific <- c("burkina_faso" = "#CC79A7", "colombia" = "#009E73","ivory_coast" = "#999999","namibia" = "#E69F00",
                      "satara" = "#56B4E9", "skukuza" = "#0072B2")
set.seed(123)
ps = readRDS(ps_rarefied_file)
#transoform data first
physeq_transformed_Bact_Arch = transform_sample_counts(ps, function (x) log(x+1))
#distance matrix
bray_PCoA = phyloseq::distance(physeq_transformed_Bact_Arch, method = "bray")
#significance of dissimilarity
sampledf_bact = data.frame(sample_data(physeq_transformed_Bact_Arch))
grouping_factor <- sampledf_bact$sampling_point 
anosim_result <- anosim(dist_matrix, grouping_factor, permutations = 999)
print(anosim_result)

bray_PCoA_ordination <- ordinate(physeq = physeq_transformed_Bact_Arch, method = "PCoA", distance = bray_PCoA)
#extract PCoA values for plotting them on the axis labels
eig_values <- bray_PCoA_ordination$values$Eigenvalues
#calculate the percentage of variance explained by the first two axes
total_variance <- sum(eig_values)
percent_var1 <- round(eig_values[1] / total_variance * 100, 2)
percent_var2 <- round(eig_values[2] / total_variance * 100, 2)

#create new axis labels
x_label <- paste("PCoA 1 (", percent_var1, "%)", sep = "")
y_label <- paste("PCoA 2 (", percent_var2, "%)", sep = "")

print(x_label)
print(y_label)
names <- c( "burkina_faso" = "Burkina Faso",
            "colombia" = "Colombia",
            "ivory_coast" = "Ivory Coast",
            "satara" = "South Africa 1 (site Satara)",
            "skukuza" = "South Africa 2 (site Skukuza)",
            "namibia" = "Namibia")

#plot the results
PCoA_bray_sampling_point <- plot_ordination(
  physeq = physeq_transformed_Bact_Arch, 
  ordination = bray_PCoA_ordination, 
  color = "sampling_point", 
  axes = 1:2) + 
  theme_bw() +
  stat_ellipse(aes(color = sampling_point), alpha = 0.6) +  #add ellipses (0.95 by default) 
  labs(x = x_label, y = y_label) +   
  theme(legend.position = "right") +  
  geom_vline(aes(xintercept = 0), linetype="dotted", linewidth=0.7) +
  geom_hline(aes(yintercept = 0), linetype="dotted", linewidth=0.7) +
  scale_color_manual(values = cbPalette, labels = names) +
  scale_fill_manual(values = cbPalette, labels = names) 
  
PCoA_bray_sampling_point


  
###ANCOM-BC2
#https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
#adapted from Lin et al. 2023 https://s3.jcloud.sjtu.edu.cn/899a892efef34b1b944a19981040f55b-oss01/bioconductor/3.18/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
library(mia)
library(phyloseq)
library(ANCOMBC)
library(tidyverse)
library(phyloseq)

ps<- readRDS(ps_unrarefied_file)
#to check the structure before 
otu_table(ps)
sample_data(ps)
all(rownames(sample_data(ps)) == colnames(otu_table(ps)))
set.seed(123)
#change sampling_point to factor
sample_data(ps)$sampling_point = factor(sample_data(ps)$sampling_point, levels = c("ivory_coast", "burkina_faso", "colombia", "namibia", "satara", "skukuza"))

#run ancombc2
out <- ancombc2(
  data = ps,  
  assay.type = "counts",  
  fix_formula = "sampling_point", 
  rand_formula = NULL, 
  p_adj_method = "holm", 
  tax_level = "Family", prv_cut = 0.2, s0_perc = 0.05, struc_zero = TRUE, neg_lb = TRUE,
  group = "sampling_point", pairwise = TRUE, mdfdr_control = list(fwer_ctrl_method = "holm", B = 1000),
  alpha = 0.05, 
  n_cl = 12,
  verbose = TRUE)
saveRDS(out, file = "ancombc2_output.rds")
out <- ancombc2_output

#pairwise comparison 
res_pair = out$res_pair

df_fig_pair1 = res_pair %>%
  dplyr::filter(diff_sampling_pointburkina_faso  == 1 |
                  diff_sampling_pointcolombia   == 1 | 
                  diff_sampling_pointnamibia   == 1 |
                  diff_sampling_pointsatara    == 1 |
                  diff_sampling_pointskukuza  == 1 |
                  diff_sampling_pointcolombia_sampling_pointburkina_faso == 1 |
                  diff_sampling_pointnamibia_sampling_pointburkina_faso     == 1 |
                  diff_sampling_pointsatara_sampling_pointburkina_faso    == 1 | 
                  diff_sampling_pointskukuza_sampling_pointburkina_faso == 1 |
                  diff_sampling_pointnamibia_sampling_pointcolombia    == 1 |
                  diff_sampling_pointsatara_sampling_pointcolombia   == 1 | 
                  diff_sampling_pointskukuza_sampling_pointcolombia == 1 |
                  diff_sampling_pointsatara_sampling_pointnamibia    == 1 |
                  diff_sampling_pointskukuza_sampling_pointnamibia   == 1 | 
                  diff_sampling_pointskukuza_sampling_pointsatara == 1 ) %>%
  dplyr::mutate(lfc1 = ifelse(diff_sampling_pointburkina_faso == 1, round(lfc_sampling_pointburkina_faso, 2), NA),
                lfc2 = ifelse(diff_sampling_pointcolombia == 1, round(lfc_sampling_pointcolombia, 2), NA),
                lfc3 = ifelse(diff_sampling_pointnamibia == 1, round(lfc_sampling_pointnamibia, 2), NA),
                lfc4 = ifelse(diff_sampling_pointsatara == 1, round(lfc_sampling_pointsatara, 2), NA),
                lfc5 = ifelse(diff_sampling_pointskukuza == 1, round(lfc_sampling_pointskukuza, 2), NA),
                lfc6 = ifelse(diff_sampling_pointcolombia_sampling_pointburkina_faso == 1, round(lfc_sampling_pointcolombia_sampling_pointburkina_faso, 2), NA),
                lfc7 = ifelse(diff_sampling_pointnamibia_sampling_pointburkina_faso == 1, round(lfc_sampling_pointnamibia_sampling_pointburkina_faso, 2), NA),
                lfc8 = ifelse(diff_sampling_pointsatara_sampling_pointburkina_faso == 1, round(lfc_sampling_pointsatara_sampling_pointburkina_faso, 2), NA),
                lfc9 = ifelse(diff_sampling_pointskukuza_sampling_pointburkina_faso == 1, round(lfc_sampling_pointskukuza_sampling_pointburkina_faso, 2), NA),
                lfc10 = ifelse(diff_sampling_pointnamibia_sampling_pointcolombia == 1, round(lfc_sampling_pointnamibia_sampling_pointcolombia, 2), NA),
                lfc11 = ifelse(diff_sampling_pointsatara_sampling_pointcolombia == 1, round(lfc_sampling_pointsatara_sampling_pointcolombia, 2), NA),
                lfc12 = ifelse(diff_sampling_pointskukuza_sampling_pointcolombia == 1, round(lfc_sampling_pointskukuza_sampling_pointcolombia, 2), NA),
                lfc13 = ifelse(diff_sampling_pointsatara_sampling_pointnamibia == 1, round(lfc_sampling_pointsatara_sampling_pointnamibia, 2), NA),
                lfc14 = ifelse(diff_sampling_pointskukuza_sampling_pointnamibia == 1, round(lfc_sampling_pointskukuza_sampling_pointnamibia, 2), NA),
                lfc15 = ifelse(diff_sampling_pointskukuza_sampling_pointsatara == 1, round(lfc_sampling_pointskukuza_sampling_pointsatara, 2), NA)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc15, names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_pair1 <- df_fig_pair1 %>%
  dplyr::mutate(sig_label = case_when(
    group == "lfc1" & q_sampling_pointburkina_faso  < 0.001 ~ "***",
    group == "lfc1" & q_sampling_pointburkina_faso  < 0.01  ~ "**",
    group == "lfc1" & q_sampling_pointburkina_faso  < 0.05  ~ "*",
    group == "lfc2" & q_sampling_pointcolombia  < 0.001 ~ "***",
    group == "lfc2" & q_sampling_pointcolombia  < 0.01  ~ "**",
    group == "lfc2" & q_sampling_pointcolombia  < 0.05  ~ "*",
    group == "lfc3" & q_sampling_pointnamibia  < 0.001 ~ "***",
    group == "lfc3" & q_sampling_pointnamibia  < 0.01  ~ "**",
    group == "lfc3" & q_sampling_pointnamibia  < 0.05  ~ "*",
    group == "lfc4" & q_sampling_pointsatara  < 0.001 ~ "***",
    group == "lfc4" & q_sampling_pointsatara  < 0.01  ~ "**",
    group == "lfc4" & q_sampling_pointsatara  < 0.05  ~ "*",
    group == "lfc5" & q_sampling_pointskukuza  < 0.001 ~ "***",
    group == "lfc5" & q_sampling_pointskukuza  < 0.01  ~ "**",
    group == "lfc5" & q_sampling_pointskukuza  < 0.05  ~ "*",
    group == "lfc6" & q_sampling_pointcolombia_sampling_pointburkina_faso  < 0.001 ~ "***",
    group == "lfc6" & q_sampling_pointcolombia_sampling_pointburkina_faso  < 0.01  ~ "**",
    group == "lfc6" & q_sampling_pointcolombia_sampling_pointburkina_faso  < 0.05  ~ "*",
    group == "lfc7" & q_sampling_pointnamibia_sampling_pointburkina_faso  < 0.001 ~ "***",
    group == "lfc7" & q_sampling_pointnamibia_sampling_pointburkina_faso  < 0.01  ~ "**",
    group == "lfc7" & q_sampling_pointnamibia_sampling_pointburkina_faso  < 0.05  ~ "*",
    group == "lfc8" & q_sampling_pointsatara_sampling_pointburkina_faso  < 0.001 ~ "***",
    group == "lfc8" & q_sampling_pointsatara_sampling_pointburkina_faso  < 0.01  ~ "**",
    group == "lfc8" & q_sampling_pointsatara_sampling_pointburkina_faso  < 0.05  ~ "*",
    group == "lfc9" & q_sampling_pointskukuza_sampling_pointburkina_faso  < 0.001 ~ "***",
    group == "lfc9" & q_sampling_pointskukuza_sampling_pointburkina_faso  < 0.01  ~ "**",
    group == "lfc9" & q_sampling_pointskukuza_sampling_pointburkina_faso  < 0.05  ~ "*",
    group == "lfc10" & q_sampling_pointnamibia_sampling_pointcolombia  < 0.001 ~ "***",
    group == "lfc10" & q_sampling_pointnamibia_sampling_pointcolombia  < 0.01  ~ "**",
    group == "lfc10" & q_sampling_pointnamibia_sampling_pointcolombia  < 0.05  ~ "*",
    group == "lfc11" & q_sampling_pointsatara_sampling_pointcolombia  < 0.001 ~ "***",
    group == "lfc11" & q_sampling_pointsatara_sampling_pointcolombia  < 0.01  ~ "**",
    group == "lfc11" & q_sampling_pointsatara_sampling_pointcolombia  < 0.05  ~ "*",
    group == "lfc12" & q_sampling_pointskukuza_sampling_pointcolombia  < 0.001 ~ "***",
    group == "lfc12" & q_sampling_pointskukuza_sampling_pointcolombia  < 0.01  ~ "**",
    group == "lfc12" & q_sampling_pointskukuza_sampling_pointcolombia  < 0.05  ~ "*",
    group == "lfc13" & q_sampling_pointsatara_sampling_pointnamibia  < 0.001 ~ "***",
    group == "lfc13" & q_sampling_pointsatara_sampling_pointnamibia  < 0.01  ~ "**",
    group == "lfc13" & q_sampling_pointsatara_sampling_pointnamibia  < 0.05  ~ "*",
    group == "lfc14" & q_sampling_pointskukuza_sampling_pointnamibia  < 0.001 ~ "***",
    group == "lfc14" & q_sampling_pointskukuza_sampling_pointnamibia  < 0.01  ~ "**",
    group == "lfc14" & q_sampling_pointskukuza_sampling_pointnamibia  < 0.05  ~ "*",
    group == "lfc15" & q_sampling_pointskukuza_sampling_pointsatara  < 0.001 ~ "***",
    group == "lfc15" & q_sampling_pointskukuza_sampling_pointsatara  < 0.01  ~ "**",
    group == "lfc15" & q_sampling_pointskukuza_sampling_pointsatara  < 0.05  ~ "*" , 
    TRUE ~ ""  ),
  label = ifelse(!is.na(value), paste0(value, sig_label), "")  )
    

df_fig_pair1$group = recode(df_fig_pair1$group, 
                           `lfc1` = "BF - IC",
                           `lfc2` = "COL - IC",
                           `lfc3` = "NAM - IC",
                           `lfc4` = "SAT - IC",
                           `lfc5` = "SKU - IC",
                           `lfc6` = "COL - BF",
                           `lfc7` = "NAM - BF",
                           `lfc8` = "SAT - BF",
                           `lfc9` = "SKU - BF",
                           `lfc10` = "NAM - COL",
                           `lfc11` = "SAT - COL",
                           `lfc12` = "SKU - COL",
                           `lfc13` = "SAT - NAM",
                           `lfc14` = "SKU - NAM",
                           `lfc15` = "SKU - SAT")
df_fig_pair1$group = factor(df_fig_pair1$group, 
                           levels = c("BF - IC",
                                      "COL - IC", 
                                      "NAM - IC",
                                      "SAT - IC",
                                      "SKU - IC", 
                                      "COL - BF",
                                      "NAM - BF",
                                      "SAT - BF", 
                                      "SKU - BF",
                                      "NAM - COL",
                                      "SAT - COL", 
                                      "SKU - COL",
                                      "SAT - NAM", 
                                      "SKU - NAM",
                                      "SKU - SAT"))



lo  <- floor(min(df_fig_pair1$value, na.rm = TRUE))
up  <- ceiling(max(df_fig_pair1$value, na.rm = TRUE))
mid <- (lo + up) / 2

taxa_order <- sort(unique(df_fig_pair1$taxon), decreasing = TRUE)
df_fig_pair1$taxon <- factor(df_fig_pair1$taxon, levels = taxa_order)


fig_pair <- df_fig_pair1 %>%
  ggplot(aes(x = group, y = taxon, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(label = label), size = 2.5) +
  labs(x = NULL, y = NULL, title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 6.5))+
  theme(plot.title = element_text(hjust = 0.5))

fig_pair

###FAPROTAX plot
#adapted from Zhang et al. 2022 https://github.com/PlantNutrition/Liyu
library(readr)
library(reshape2)
library(doBy)
library(pheatmap)

metadata <- readRDS(metadata_file)
func_table_rel <- read_delim("RESULTS/function/func_table.csv", delim = "\t")
func_table_rel <- as.data.frame(func_table_rel) 
rownames(func_table_rel) <- func_table_rel$group
original_order <- rownames(func_table_rel)    
func_table_rel <- func_table_rel[, -1]
func_table_rel$func_table_rel <- rownames(func_table_rel)
func_melt <- melt(func_table_rel,id.vars = c("func_table_rel"),variable.name = "Sample_ID",value.name = "RA")
group <- metadata
func_merge <-merge(func_melt,group,by="Sample_ID")
func_mean <- summaryBy(RA~func_table_rel+sampling_point, func_merge, FUN = c(mean))
names(func_mean) <- c('func_table_rel', 'sampling_point', 'RA')

func_dcast <- dcast(func_mean,func_table_rel~sampling_point,value.var="RA")
rownames(func_dcast)<-func_dcast[,1]
func_dcast<-func_dcast[,-1]

#scale_test <- apply(func_dcast, 2, function(x) log2(x+1))
log_fap <- log2(func_dcast + 1e-6)
z_scaled_func <- t(scale(t(func_dcast)))
z_scaled_func <- z_scaled_func[original_order[original_order %in% rownames(z_scaled_func)], ]

names <- c( "burkina_faso" = "Burkina Faso",
            "colombia" = "Colombia",
            "ivory_coast" = "Ivory Coast",
            "satara" = "South Africa 1 (site Satara)",
            "skukuza" = "South Africa 2 (site Skukuza)",
            "namibia" = "Namibia")
colnames(z_scaled_func) <- names[colnames(z_scaled_func)]

p <-  pheatmap(z_scaled_func, cluster_row = FALSE, cluster_col = TRUE, cutree_col = 2, show_colnames = TRUE, show_rownames = TRUE)

p

###ENV VARIABLES PCA
#https://www.davidzeleny.net/anadat-r/doku.php/en:pca_examples
library(vegan)
library(ggplot2)

cbPalette <- c("#CC79A7","#009E73", "#999999",  "#E69F00", "#56B4E9", "#0072B2",  "#F0E442",  "#D55E00" )
#import the table with environmental variables
#import the table with environmental variables
table_environ <- readRDS(metadata_file)

#select variable that should be tested (only continuous variables)
env.bio <- table_environ[, c(17, 19:22)]
sampling_point <- table_environ$sampling_point


#runPCA
prc <- rda (env.bio, scale = TRUE)
head (summary (prc))
#to get PCA coordinates of samples
pca_scores <- as.data.frame(prc$CA$u) 
colnames(pca_scores) <- c("PC1", "PC2")
pca_scores$sampling_point <- sampling_point

labels <- c("burkina_faso" = "Burkina Faso",
            "colombia" = "Colombia",
            "ivory_coast" = "Ivory Coast",
            "namibia" = "Namibia",
            "satara" = "South Africa 1 (site Satara)",
            "skukuza" = "South Africa 2 (site Skukuza)")

# Create the PCA plot with arrows
p <- ggplot(pca_scores, aes(x = PC1, y = PC2, color=sampling_point)) +
  geom_point(size = 1) +  
  theme_bw() +
  labs(title = "", x = x_label, y = y_label) +
  theme(legend.position = "right") +
  #annotate('text', label = expression(italic("R"^2) ~ "= 0.839, p-value < 0.001"), x = 0.05, y = 0.25, size =3) +
  scale_color_manual(values = cbPalette, labels = labels) +
  scale_fill_manual(values = cbPalette, labels = labels) +
  geom_vline(aes(xintercept = 0), linetype = "dotted", linewidth = 0.7) +
  geom_hline(aes(yintercept = 0), linetype = "dotted", linewidth = 0.7) +
  stat_ellipse(aes(group= sampling_point), level = 0.95,  size = 0.5) 
  
p

#extract PCA values to plot them
eig_values <- prc$CA$eig
total_variance <- sum(eig_values)

#calculate variance explained by PC1 and PC2
percent_var1 <- round(eig_values[1] / total_variance * 100, 2)
percent_var2 <- round(eig_values[2] / total_variance * 100, 2)

x_label <- paste("PC 1 (", percent_var1, "%)", sep = "")
y_label <- paste("PC 2 (", percent_var2, "%)", sep = "")

print(x_label)
print(y_label)

#calculate the significane and R-sured of sample distribution 
dtp <- data.frame('Site' = sampling_point, prc$CA$u[, 1:2])  #PCA scores from vegan
PCAcoords <- dtp[, c("PC1", "PC2")]
PCAdist <- dist(PCAcoords)
results_env_PCA <- adonis2(PCAdist  ~ Site, data = dtp, permutations = 999)
print(results_env_PCA)


###CORELATIONS
#information obtained from: https://rpubs.com/MajstorMaestro/240657
#https://stackoverflow.com/questions/59018156/position-of-two-regression-equations-ggplot-r
#https://rcompanion.org/rcompanion/e_01.html

cbPalette <- c("#CC79A7","#009E73", "#999999",  "#E69F00", "#56B4E9", "#0072B2",  "#F0E442",  "#D55E00" )

library(tidyverse)
library(phyloseq)
library(ggpubr)
library(vegan)
library(ggplot2)
library(dplyr)
library(corrplot)
library(Hmisc)
library(ggrepel)
library(ggpmisc)
library(car)
library(rstatix)
library(PerformanceAnalytics)

#multiple linear regression
#first remove <LD from datase
table_environ <- readRDS(metadata_file)
table_environ$available_P[table_environ$available_P < 3.71] <- NA
env.bio <- table_environ[, c(17:22)]
Sample_ID <- table_environ$Sample_ID
sampling_point <- table_environ$sampling_point
env <- cbind(Sample_ID, sampling_point, env.bio)

#outliers as NA (outlier were identified by Tueky method)
table_alpha_diversity <- alpha_diversity
alpha.div <- table_alpha_diversity[, c(23,24, 26:28)]
Sample_ID <- table_alpha_diversity$Sample_ID
sampling_point <- table_alpha_diversity$sampling_point
transf_alpha_div <- cbind(Sample_ID, sampling_point, alpha.div_norm)

merged_data_env_div <- merge(transf_alpha_div, env, by = "Sample_ID")

merged_data_env_div <- read_csv("RESULTS/merged_data_env_div.csv")
my_data <- select(merged_data_env_div, Observed, Shannon, InvSimpson, pH, organic_matter, organic_C, available_P, total_N, organic_C_N_ratio)
cor_1 <- cor(my_data, method = "spearman")
cor_2 <- rcorr(as.matrix(my_data, method = "spearman"))
cor_p <- cor_2$P
cor_p_df <- data.frame(cor_p)
write.table(cor_p_df, file = "spearman_cor_p_out2.txt", row.names = FALSE, col.names = TRUE, sep = "\t")
cor_R <- cor_2$r
cor_R_df <- data.frame(cor_R)
write.table(cor_R_df, file = "spearman_cor_Rout2.txt", row.names = FALSE, col.names = TRUE, sep = "\t")
cor_n <- cor_2$n
cor_n_df <- data.frame(cor_n)
write.table(cor_n_df, file = "spearman_cor_n_out2.txt", row.names = FALSE, col.names = TRUE, sep = "\t")

#for each sampling_point separately
my_data_selected <- dplyr::select(my_data, Observed, Shannon, InvSimpson, pH, organic_matter, 
                           organic_C, available_P, total_N, organic_C_N_ratio, sampling_point)
unique_sampling_points <- unique(my_data_selected$sampling_point)
cor_p_list <- list()
cor_r_list <- list()
cor_n_list <- list()

for (point in unique_sampling_points) {
  subset_data <- my_data_selected %>%
    filter(sampling_point == point) %>%
    select(-sampling_point) 
  
  cor_result <- rcorr(as.matrix(subset_data), type = "spearman")
  
  cor_p_list[[point]] <- data.frame(cor_result$P)
  cor_r_list[[point]] <- data.frame(cor_result$r)
  cor_n_list[[point]] <- data.frame(cor_result$n)
  
  write.table(cor_p_list[[point]], file = paste0("spearman_cor_p_", point, ".txt"), 
              row.names = TRUE, col.names = TRUE, sep = "\t")
  write.table(cor_r_list[[point]], file = paste0("spearman_cor_r_", point, ".txt"), 
              row.names = TRUE, col.names = TRUE, sep = "\t")
  write.table(cor_n_list[[point]], file = paste0("spearman_cor_n_", point, ".txt"), 
              row.names = TRUE, col.names = TRUE, sep = "\t")
}


#to show only significant 
k <- ggplot(filtered_data, aes(pH, Shannon, color = sampling_point)) +
  geom_point() +
  labs(x = "Organic C [%]", y = "Inverse Simpson index") +
  stat_smooth(aes(fill = sampling_point, color = sampling_point), method = "lm", formula = model) +
  stat_poly_eq(
    aes(label = paste( ..p.value.label.., ..rr.label.., sep = "~~~~")), 
    rsquared.conf.level = 0.95,
    formula = model, 
    size = 3,
    label.x = "left",  
    label.y = "bottom",
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = palette_specific) +  
  scale_fill_manual(values = palette_specific) +
  theme_bw()

k

selected_sampling_points <- c("namibia", "ivory_coast", "colombia", "burkina_faso", "satara", "skukuza")
filtered_data <- my_data[my_data$sampling_point %in% selected_sampling_points, ]

#check the model correlation
model <- lm(Shannon ~ pH  , data = my_data)
my$data
summary(model)
par(mfrow = c(2, ))
plot(model) 
shapiro.test(residuals(model))

hist(residuals(model))    
#homogenity and normality of residuals check
par(mfrow = c(2,2))
qqnorm(resid(model), las=1)  
qqline(resid(model), col = "red")
hist(resid(model), las=1)  
plot(fitted(model), residuals(model))
     abline(h = 0, col = "red", lwd = 2) 

shapiro.test(residuals(model))
residuals_data <- residuals(model)
bartlett.test(residuals(model) ~ alpha_diversity$sampling_point)


leveneTest(residuals(model) ~ alpha_diversity$sampling_point)


###MODELS
###model using stepAIC package
#https://www.sthda.com/english/articles/37-model-selection-essentials-in-r/154-stepwise-regression-essentials-in-r/
#fit the full model 
#use data from previous step
full.model <- lm((InvSimpson) ~., data = my_data)
step.model <- stepAIC(full.model, direction = "both",  trace = FALSE)

summary(full.model)
summary(step.model)

#model using traincontrol 
train.control <- trainControl(method = "cv", number = 10)
step.model <- train((InvSimpson) ~., data = my_data, na.rm = TRUE,
                    method = "leapSeq", 
                    tuneGrid = data.frame(nvmax = 1:5),
                    trControl = train.control)

#model accuracy
step.model$results
step.model$bestTune
#final model coefficients
step.model$finalModel
#summary of the model
summary(step.model$finalModel)

#regression coeficients
coef(step.model$finalModel, 4)

final_model_fit <- lm((InvSimpson) ~  pH +  organic_C + total_N +  organic_C_N_ratio +  sampling_point + pH:sampling_point +
                        organic_C:sampling_point + total_N:sampling_point + organic_C_N_ratio:sampling_point,  data = my_data)

summary(final_model_fit)

stepwise_model <- stepAIC(final_model_fit, direction = "both")

#choose the best option
final_model_fit <- lm((InvSimpson) ~ total_N + sampling_point + total_N:sampling_point,  data = my_data)

par(mfrow=c(2,3))
plot(final_model_fit)

Anova(final_model_fit,Type="II")
hist(residuals(final_model_fit))
shapiro.test(residuals(final_model_fit))
bptest(final_model_fit)
bartlett.test(residuals(final_model_fit) ~ as.factor(my_data$sampling_point))
leveneTest(residuals(final_model_fit) ~ as.factor(my_data$sampling_point))


# test
#leverage and cooks distance
leverage <- hatvalues(final_model_fit)
cooks_distance <- cooks.distance(final_model_fit)
plot_data <- data.frame(
  Sample = rownames(my_data),
  Leverage = leverage,
  CooksDistance = cooks_distance
)

ggplot(plot_data, aes(x = Leverage, y = CooksDistance)) +
  geom_point(color = "blue") + 
  geom_text(aes(label = Sample), hjust = 1.5, vjust = 1.5, size = 3, color = "red") +  
  labs(title = "Cook's Distance vs Leverage",
       x = "Leverage",
       y = "Cook's Distance") +
  theme_minimal() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  #cooks distance threshold
  geom_vline(xintercept = 2 * mean(leverage), linetype = "dashed", color = "red")  #leverage threshold
  
###BARPLOTS AT PHYLUM AND CLASS LEVEL
#https://www.geeksforgeeks.org/r-language/how-to-save-time-with-data-visualization-using-stack-in-r-with-ggplot2/
library(ggplot2)
library(dplyr)
library(phyloseq)
library(microViz)
library(forcats)
library(pals)
library(stringr)
#collapse the phyloseq object at the Phylum level + sampling_sites
ps = readRDS(ps_rarefied_file)
#CLEAN PS OBJECT ----
ps <- ps %>% 
  tax_mutate(
    Species = case_when(
      !is.na(Genus) & !is.na(Species) ~ paste(sub('.*__', '', Genus), sub('.*__', '', Species), sep = " "),
      !is.na(Genus) & is.na(Species) ~ paste(sub('.*__', '', Genus), "sp.", sep = " "),
      TRUE ~ "Unknown"),
    Species = sub(".*__","",Species),
    Genus = sub(".*__","",Genus),
    Family = sub(".*__","",Family),
    Order = sub(".*__","",Order),
    Class = sub(".*__","",Class),
    Phylum = sub(".*__","",Phylum),
    Kingdom = sub(".*__","",Kingdom)) %>% 
  tax_fix() %>%
  phyloseq_validate(verbose=TRUE)

#agglomerate taxa
ps_phylum <- tax_glom(ps, taxrank = "Phylum", NArm = FALSE)
ps_df <- psmelt(ps_phylum)

#group by sampling_point and Phylum, sum the counts per group
ps_grouped <- ps_df %>%
  group_by(sampling_point, Phylum, ) %>%
  summarise(Abundance = sum(Abundance))

#calculate total abundance for each sampling point
ps_grouped_total <- ps_grouped %>%
  group_by(sampling_point) %>%
  summarise(Total_abundance = sum(Abundance))

#merge total abundance with grouped data to calculate relative abundance
ps_grouped <- ps_grouped %>%
  left_join(ps_grouped_total, by = "sampling_point") %>%
  mutate(Relative_abundance = Abundance / Total_abundance)

ps_total_abundance <- ps_grouped %>%
  group_by(Phylum) %>%
  summarise(OverallRelativeAbundance = sum(Abundance) / sum(Total_abundance), .groups = 'drop')


#identify phyla with overall relative abundance < 1% and group them as "< 0.01"
low_abund_label <- "< 0.01"
ps_grouped <- ps_grouped %>%
  left_join(ps_total_abundance, by = "Phylum") %>%
  mutate(Phylum = ifelse(OverallRelativeAbundance < 0.01, low_abund_label, Phylum))

ps_grouped$Phylum[is.na(ps_grouped$Phylum)] <- low_abund_label
#reorder based on overall rel. abundance
ps_total_abundance_updated <- ps_grouped %>%
  group_by(Phylum) %>%
  summarise(OverallRelativeAbundance = sum(Abundance) / sum(Total_abundance), .groups = 'drop')

ps_grouped_reorder <- ps_total_abundance_updated %>%
  arrange(desc(OverallRelativeAbundance)) %>%
  pull(Phylum)
ps_grouped_reorder <- c(setdiff(ps_grouped_reorder, low_abund_label), low_abund_label)
ps_grouped$Phylum <- factor(ps_grouped$Phylum, levels = ps_grouped_reorder)

#to reorder sampling points
ps_grouped <- ps_grouped %>%
  group_by(sampling_point) %>%
  mutate(sampling_point = fct_reorder(sampling_point, Relative_abundance, .fun = sum, .desc = FALSE))



#plot relative abundance at phylum, grouped by sampling point
p <-ggplot(ps_grouped, aes(x = sampling_point, y = Relative_abundance, fill = (Phylum ))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_discrete(labels = c("Burkina Faso", "Colombia", "Ivory Coast", "Namibia", "South Africa 1 (site Satara)", "South Africa 2 (site Skukuza)")) +
  scale_fill_manual(values = c( "#99cc99","#cc6666","orchid4",  "#E69F00","gold1","darkcyan", "lightpink", "#ff9966",  "#56B4E9",
                                "springgreen4", "firebrick","#666666","lightblue", "mediumvioletred", 
                                  "#99cc00",   "darkblue")) +
  theme_minimal() +
  theme(legend.title=element_blank()) +
  labs(y = "Relative abundances", x = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10, vjust=1))+
  theme(panel.grid.major = element_blank())
p

#collapse the phyloseq object at the Class level + sampling points
#begin with the ps object created as in the first case

#agglomerate taxa
ps_class <- tax_glom(ps, taxrank = "Class", NArm = FALSE)
ps_df <- psmelt(ps_class)

#group by sampling point and Class, sum the counts per group
ps_grouped <- ps_df %>%
  group_by(sampling_point, Class) %>%
  summarise(Abundance = sum(Abundance))

#calculate the total abundance for each sampling point
ps_grouped_total <- ps_grouped %>%
  group_by(sampling_point) %>%
  summarise(Total_abundance = sum(Abundance))

#merge total abundance with grouped data to calculate relative abundance
ps_grouped <- ps_grouped %>%
  left_join(ps_grouped_total, by = "sampling_point") %>%
  mutate(Relative_abundance = Abundance / Total_abundance)

ps_total_abundance <- ps_grouped %>%
  group_by(Class) %>%
  summarise(OverallRelativeAbundance = sum(Abundance) / sum(Total_abundance))

#identify taxa with overall relative abundance < 1% and group them as "< 0.01"
low_abund_label <- "< 0.01"
ps_grouped <- ps_grouped %>%
  left_join(ps_total_abundance, by = "Class") %>%
  mutate(Class = ifelse(OverallRelativeAbundance < 0.01, low_abund_label, Class))

ps_grouped$Class[is.na(ps_grouped$Class)] <- low_abund_label
#reorder based on overall rel. abundance
ps_total_abundance_updated <- ps_grouped %>%
  group_by(Class) %>%
  summarise(OverallRelativeAbundance = sum(Abundance) / sum(Total_abundance), .groups = 'drop')

ps_grouped_reorder <- ps_total_abundance_updated %>%
  arrange(desc(OverallRelativeAbundance)) %>%
  pull(Class)
ps_grouped_reorder <- c(setdiff(ps_grouped_reorder, low_abund_label), low_abund_label)
ps_grouped$Class <- factor(ps_grouped$Class, levels = ps_grouped_reorder)
#to reorder sampling points
ps_grouped <- ps_grouped %>%
  group_by(sampling_point) %>%
  mutate(sampling_point = fct_reorder(sampling_point, Relative_abundance, .fun = sum, .desc = TRUE))

#plot relative abundance at Class, grouped by sampling point
p <- ggplot(ps_grouped, aes(x = sampling_point, y = Relative_abundance, fill = Class)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_discrete(labels = c( "Burkina Faso","Colombia","Ivory Coast", "Namibia","South Africa 1 (site Satara)","South Africa 2 (site Skukuza)")) +
  scale_fill_manual(values = c("#cc6666", "springgreen4", "#ff9966", "orchid4", "mediumvioletred","lightpink", "#99cc00", "darkblue", "gold1","darkcyan",
                              "lightblue" ,"firebrick", "#E69F00","#99cc99", "#56B4E9","#666666")) +
  theme_minimal() +
  theme(legend.title=element_blank()) +
  labs(y = "Relative abundances", x = "") +
  theme(axis.text.x = element_text(angle = 90,hjust = 1, size =10, vjust=1))+
  theme(axis.text.y = element_text( hjust = 1, size =10, vjust=0.1)) +
  theme(panel.grid.major = element_blank())

p 

###db-RDA analysis
#https://github.com/gvMicroarctic/AntarcticBiogeographyPaper/blob/main/dbRDA_variation_partitioning.R
# which is at the based on this https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html
#https://github.com/KatjaKo/PacBio_AMF/blob/v1.2/scripts/script_multivariate_stats.Rmd

library(vegan)
library(phyloseq)
library(ggplot2)
library(dplyr)
set.seed(123)
ps <- ps_rarefied
#taxa data prepare
ps_taxa <- tax_glom(ps, taxrank = "Genus")
#rel abudnaces
ps_taxa_rel <- transform_sample_counts(ps_taxa, function(x) x / sum(x))

rel_abund <- otu_table(ps_taxa_rel)
if (!taxa_are_rows(rel_abund)) {
  rel_abund <- t(rel_abund)
}
#hellinger transformed genera dataset
taxa_abund_hellinger <- decostand(rel_abund, method = "hellinger")
taxa_abund_hellinger <- as.data.frame(t(taxa_abund_hellinger))
dim(taxa_abund_hellinger)

#standardize env variables
env_vars <- c("pH", "total_N", "organic_C", "organic_matter", "organic_C_N_ratio")
env <- metadata[, env_vars]
env$Sample_ID <- metadata$Sample_ID
rownames(env) <- env$Sample_ID
env$Sample_ID <- NULL
env_stand <- decostand(env, method = "standardize")
round(apply(env_stand, 2, mean), 1)
apply(env_stand, 2, sd)
env_stand <- data.frame(env_stand)

#sanity check, all must be in the same order
env_stand <- env_stand[rownames(taxa_abund_hellinger), ]
all(rownames(env_stand) == rownames(taxa_abund_hellinger))
stopifnot(all(rownames(env_stand) == rownames(taxa_abund_hellinger)))


#run RDA
dbrda_upr <- capscale(taxa_abund_hellinger ~ ., data = data.frame(env_stand), dist="bray") # upper model limit (the "full" model; it takes in account all variables)
dbrda_lwr <- capscale(taxa_abund_hellinger ~ 1, data = data.frame(env_stand), dist="bray") # lower model limit

model <- ordiR2step(dbrda_lwr, #lower model limit
                    scope = formula(dbrda_upr), #upper model limit
                    direction = "both",
                    R2scope = TRUE, #can't surpass the "full" model's R2
                    pstep = 999,
                    trace = TRUE) #change to TRUE to see the selection process

model$call

dbrda_final <- capscale(formula = taxa_abund_hellinger ~ pH +  organic_matter + organic_C_N_ratio, data = data.frame(env_stand), distance = "bray")
summary(dbrda_final)

vif_values <- vif.cca(dbrda_final)
print(vif_values)


#get statistical significance of the model
#same results as anova()
anova.cca(model, permutations = 999) 
#get statistical significance of the model for each variable
anova.cca(model, step = 999, by = "term")
anova.cca(model, step = 999, by = "axis")

#check constrained and unconstrained variance
RsquareAdj(dbrda_final)
#check constrained and unconstrained variance again
dbrda_final$CCA$tot.chi/dbrda_final$tot.chi
constrained_eig <- dbrda_final$CCA$eig/dbrda_final$tot.chi*100
unconstrained_eig <- dbrda_final$CA$eig/dbrda_final$tot.chi*100
expl_var <- c(constrained_eig, unconstrained_eig)
barplot (expl_var[1:20], col = c(rep ('red', length (constrained_eig)), rep ('black', length (unconstrained_eig))),
         las = 2, ylab = '% variation')


#plot rda
perc <- round(100*(summary(dbrda_final)$cont$importance[2, 1:2]), 2)

#extract scores - these are coordinates in the RDA space
sc_si <- scores(dbrda_final, display="sites", choices=c(1,2), scaling=1)
sc_sp <- scores(dbrda_final, display="species", choices=c(1,2), scaling=1)
sc_bp <- scores(dbrda_final, display="bp", choices=c(1, 2), scaling=1)

plot(dbrda_final,
     scaling = 1, 
     type = "none", 
     frame = FALSE,
     # set axis limits
     xlim = c(-1.2, 1.2), 
     ylim = c(-1.2, 1.2),
     xlab = paste0("RDA1 (", perc[1], "%)"), 
     ylab = paste0("RDA2 (", perc[2], "%)"))

#to give colors based sites
sites <- metadata$sampling_point[match(rownames(sc_si), metadata$Sample_ID)]
sites <- factor(sites)
stopifnot(length(sites) == nrow(sc_si))
table(sites)
palette_specific <- c("burkina_faso" = "#CC79A7", "colombia" = "#009E73","ivory_coast" = "#999999","namibia" = "#E69F00",
                      "satara" = "#56B4E9", "skukuza" = "#0072B2")

ordihull(dbrda_final, sites, draw = "lines", scaling = 2, col = palette_specific, alpha = 40, lty = 0)
#centres and areas of the hulls
summary(pl)
#add points for samples
points(sc_si, pch = 21,  bg = palette_specific[sites], cex = 0.8)

set.seed(123)
#extract env variables
envfit_res <- envfit(dbrda_final ~ pH + organic_C_N_ratio + organic_matter, data = env_stand, permutations = 999)
envfit_res
envfit_df <- as.data.frame(scores(envfit_res, display = "vectors"))
envfit_df$Variable <- rownames(envfit_df)
envfit_df$pval <- envfit_res$vectors$pvals
head(envfit_df)
envfit_sig <- envfit_df %>% filter(pval <= 0.05)

#add env varibales
arrows(0, 0, envfit_sig$CAP1, envfit_sig$CAP2, col = "black", lwd = 1, length = 0.1)
text(envfit_sig$CAP1 * 1.1, envfit_sig$CAP2 * 1.1, labels = envfit_sig$Variable, col = "black", cex = 0.8)

#text(sc_sp[,1], sc_sp[,2], labels = rownames(sc_sp),  col = "blue", cex = 0.7, font = 3)  
#text(sc_sp[,1] + 0.03, sc_sp[,2] + 0.03, labels = rownames(sc_sp), col = "blue", cex = 0.7, font = 3)

#rename columns to something clearer -this is only visualization. it is not correlation between env variable and phyla. I only extracted values assesed for genera,
#then classified them at the phylum level, calculated centroid and ploted them in relation to the main rda axes
colnames(sc_sp) <- c("CAP1", "CAP2")
tax_table_df <- as.data.frame(tax_table(ps_taxa_rel))

#match and combine
taxa_phylum <- tax_table_df$Phylum[match(rownames(sc_sp), rownames(tax_table_df))]
phylum_scores <- data.frame(sc_sp, Phylum = taxa_phylum)

phylum_avg <- phylum_scores %>%
  group_by(Phylum) %>%
  summarise(
    CAP1 = mean(CAP1, na.rm = TRUE),
    CAP2 = mean(CAP2, na.rm = TRUE)) %>%
  mutate(arrow_length = sqrt(CAP1^2 + CAP2^2)) %>%
  filter(arrow_length > 0.3) 

#plot phyla arrows
arrows(0, 0, phylum_avg$CAP1, phylum_avg$CAP2, col = "red", lwd = 1, length = 0.1)
text(phylum_avg$CAP1 * 1.1, phylum_avg$CAP2 * 1.1, labels = phylum_avg$Phylum, col = "red", cex = 0.8, font = 3)



