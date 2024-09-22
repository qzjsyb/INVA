# load library  -------------------------------------------------------------
library(tidyverse)
library(mixOmics)
library(gridExtra)
library(lme4)
library(lmerTest)
library(performance)
library(ggrepel)
library(vegan)
library(stats)
setwd('C:/Users/mekay/Desktop/Thesis/R code/KEGG')
metadata <- read.csv("metadata_maldi.csv", header = TRUE, row.names = 1)
metadata <- metadata[,c(1,2)]
compounds <- read.csv("Compound_Link.csv", header = TRUE, row.names = 1)

#Load Functions ------------------------------------
#Pareto Scaling:
PS_helper <- function(x) {
  (x - mean(x)) / sqrt(sd(x, na.rm = T))
}	

pareto_scale <- function(x){
  mtb_scaled <- data.frame(apply(x, 2, PS_helper))
  return(mtb_scaled)
}

### can be sensitive to large range. The difference between Z-score is it squar root the sd again.

#Auto Scaling:
AS_helper <- function(x) {
  (x - mean(x)) / sd(x, na.rm = T)
} 

auto_scale <- function(x){
  mtb_scaled <- apply(x, 2, AS_helper) 
  return(mtb_scaled)
}

#Log Transformation Functions:
log_helper <- function(x, min.val) {
  log2((x + sqrt(x ^ 2 + min.val ^ 2)) / 2)
}


#Log Scaling:
log_transform <- function(x){
  x_nz <- x[ ,which(apply(x, 2, sum) != 0)] # remove 0 features
  min.val <- min(abs(x_nz[x_nz!=0]))/10
  x_log_trans <- data.frame(apply(x_nz, 2, log_helper, min.val))
  return(x_log_trans)
}

set.seed(123)



## APEBA---------------------
APEBA <- read_csv("APEBA.csv") %>%
  #Make the compound IDs the row names
  column_to_rownames('ID') %>%
  #Transpose the matrix to make it so samples are rows and compounds are columns
  t()  %>%
  as.data.frame()

#boxplot(APEBA)

APEBA_pareto <- APEBA %>%
  log_transform() %>%
  pareto_scale()


#Combine metadata with metabolite data
APEBA_final  <- cbind(APEBA_pareto,metadata)

#NMDS plot
set.seed(123)

#plot NMDS
# Compute distance matrix
X <- APEBA_final  %>%
  dplyr::select(starts_with('APEBA'))
dist_matrix <- vegdist(X, method = "euclidean")

# Perform NMDS
nmds_result <- metaMDS(dist_matrix)

# Add NMDS coordinates to dataframe
APEBA_final$NMDS1 <- nmds_result$points[, 1]
APEBA_final$NMDS2 <- nmds_result$points[, 2]

# Plot NMDS
plot_nmds <- ggplot(APEBA_final, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(shape = CanopyStatus, color = Plant, fill = Plant), size = 3) +
  theme_bw(base_size = 15) +
  scale_color_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_fill_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_shape_manual(values = c("Canopy" = 16, "Open" = 1)) +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  annotate("text", x = Inf, y = -Inf, vjust = 0, hjust = 1, label = paste("Stress =", round(nmds_result$stress, 3)))

plot_nmds

# Perform PERMANOVA
permanova_result <- adonis2(dist_matrix ~ Plant * CanopyStatus, data = APEBA_final)

# Round numeric values in PERMANOVA result
permanova_table <- as.data.frame(permanova_result)
permanova_table[, sapply(permanova_table, is.numeric)] <- round(permanova_table[, sapply(permanova_table, is.numeric)], 3)

# Create table plot
table_plot <- tableGrob(permanova_table)

# Save NMDS plot and table
ggsave("NMDS_APEBA.png", 
       grid.arrange(plot_nmds + ggtitle("APEBA"), table_plot, ncol = 1, heights = c(3, 1)), 
       width = 5, height = 7)

# differential analysis ------------
p.plant <- apply(APEBA_final[, 1:151], 2, function(x){kruskal.test(x ~ Plant, data = APEBA_final)$p.value})
#p.adj <- p.adjust(p.plant, method = "fdr")

#calculate p value

APEBA_LC2  <- cbind(APEBA, metadata)

DICA <- APEBA_LC2[which(APEBA_LC2$Plant == "DICA"),1:ncol(APEBA)]
ELRE <- APEBA_LC2[which(APEBA_LC2$Plant == "ERLE"),1:ncol(APEBA)]
LC2_APEBA <- log2(colMeans(DICA)/colMeans(ELRE))
volcano_APEBA.kw <- data.frame(pvalue = p.plant,
                               Log2FC = LC2_APEBA) 


Plot_volcano <- function(dataframe) {
  # Create a copy of the input dataframe
  volcano_data <- dataframe
  filename <- deparse(substitute(dataframe))
  
  # Perform data transformations and calculations
  volcano_data$diffexpressed <- "NO"
  volcano_data$diffexpressed[volcano_data$Log2FC > 2 & volcano_data$pvalue < 0.05] <- "DICA"
  volcano_data$diffexpressed[volcano_data$Log2FC < -2 & volcano_data$pvalue < 0.05] <- "ERLE"
  volcano_data$delabel <- NA
  volcano_data$ID <- row.names(volcano_data)
  compounds$ID <- row.names(compounds)
  volcano_data <- left_join(volcano_data,compounds, by = "ID")
  
  volcano_data$delabel <- volcano_data$moleculeNames
  volcano_data$delabel[volcano_data$diffexpressed == "NO"] <- NA
  
  # Save the processed dataframe with a name corresponding to the input dataframe
  assign(filename, volcano_data, envir = .GlobalEnv)
  
  # Create volcano plot
  plot_volcano <- ggplot(data = volcano_data, aes(x = Log2FC, y = -log10(pvalue), col = diffexpressed)) + 
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", col = "red", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "red", alpha = 0.5) +
    geom_point() + 
    theme_minimal() +
    scale_color_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A", "NO" = "grey"))
  
  return(plot_volcano)
} 

Plot_volcano(volcano_APEBA.kw)

open <- APEBA_final[APEBA_final$CanopyStatus == "Open",]
canopy <- APEBA_final[APEBA_final$CanopyStatus == "Canopy",]

p.plant.o <- apply(open[, 1:151], 2, function(x){kruskal.test(x ~ Plant, data = open)$p.value})
p.adj.o <- p.adjust(p.plant.o, method = "fdr")
p.plant.c <- apply(canopy[, 1:151], 2, function(x){kruskal.test(x ~ Plant, data = open)$p.value})
p.adj.c <- p.adjust(p.plant.c, method = "fdr")

e.open <- APEBA_LC2[APEBA_LC2$Plant == "ERLE" & APEBA_LC2$CanopyStatus == "Open",]
e.canopy <- APEBA_LC2[APEBA_LC2$Plant == "ERLE" & APEBA_LC2$CanopyStatus == "Canopy",]
d.open <- APEBA_LC2[APEBA_LC2$Plant == "DICA" & APEBA_LC2$CanopyStatus == "Open",]
d.canopy <- APEBA_LC2[APEBA_LC2$Plant == "DICA" & APEBA_LC2$CanopyStatus == "Canopy",]

LC2_open <- log2(colMeans(d.open[,c(1:151)])/colMeans(e.open[,c(1:151)]))
LC2_canopy <- log2(colMeans(d.canopy[,c(1:151)])/colMeans(e.canopy[,c(1:151)]))

volcano_APEBA.open<-data.frame(pvalue = p.plant.o,
                               Log2FC = LC2_open) 
volcano_APEBA.canopy<-data.frame(pvalue = p.plant.c,
                                 Log2FC = LC2_canopy) 

Plot_volcano(volcano_APEBA.canopy)
Plot_volcano(volcano_APEBA.open)


APEBA_LC2 %>%
  ggplot(aes(x = Plant, y = APEBA045)) +
  geom_jitter(aes(color = Plant, shape = CanopyStatus), position = position_jitterdodge(dodge.width = 0.75)) +
  geom_boxplot(aes(color = Plant, fill = Plant, shape = CanopyStatus, alpha = CanopyStatus), outlier.shape = NA) +
  scale_fill_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_color_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_shape_manual(values = c("Canopy" = 16, "Open" = 1)) +
  scale_alpha_manual(values = c("Canopy" = 0.7, "Open" = 0)) +
  labs(title = "Root_asparagine_MALDI_APEBA") +
  xlab(NULL) +
  ylab("MALDI_intensity_average") +
  theme_bw(base_size = 15) +
  theme(legend.position = "none")

APEBA_LC2 %>%
  ggplot(aes(x = Plant, y = APEBA058)) +
  geom_jitter(aes(color = Plant, shape = CanopyStatus), position = position_jitterdodge(dodge.width = 0.75)) +
  geom_boxplot(aes(color = Plant, fill = Plant, shape = CanopyStatus, alpha = CanopyStatus), outlier.shape = NA) +
  scale_fill_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_color_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_shape_manual(values = c("Canopy" = 16, "Open" = 1)) +
  scale_alpha_manual(values = c("Canopy" = 0.7, "Open" = 0)) +
  labs(title = "Root_glutamine_MALDI_APEBA") +
  xlab(NULL) +
  ylab("MALDI_intensity_average") +
  theme_bw(base_size = 15) +
  theme(legend.position = "none")

## NEDC---------------------
NEDC <- read_csv("NEDC.csv") %>%
  #Make the compound IDs the row names
  column_to_rownames('ID') %>%
  #Transpose the matrix to make it so samples are rows and compounds are columns
  t()  %>%
  as.data.frame()

NEDC_pareto <- NEDC %>%
  log_transform() %>%
  pareto_scale()

#Combine metadata with metabolite data
NEDC_final  <- cbind(NEDC_pareto,metadata)

#NMDS plot
set.seed(123)

#plot NMDS
# Compute distance matrix
X <- NEDC_final  %>%
  dplyr::select(starts_with('NEDC'))
dist_matrix <- vegdist(X, method = "euclidean")

# Perform NMDS
nmds_result <- metaMDS(dist_matrix)

# Add NMDS coordinates to dataframe
NEDC_final$NMDS1 <- nmds_result$points[, 1]
NEDC_final$NMDS2 <- nmds_result$points[, 2]

# Plot NMDS
plot_nmds <- ggplot(NEDC_final, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(shape = CanopyStatus, color = Plant, fill = Plant), size = 3) +
  theme_bw(base_size = 15) +
  scale_color_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_fill_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_shape_manual(values = c("Canopy" = 16, "Open" = 1)) +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  annotate("text", x = Inf, y = -Inf, vjust = 0, hjust = 1, label = paste("Stress =", round(nmds_result$stress, 3)))

plot_nmds

# Perform PERMANOVA
permanova_result <- adonis2(dist_matrix ~ Plant * CanopyStatus, data = NEDC_final)

# Round numeric values in PERMANOVA result
permanova_table <- as.data.frame(permanova_result)
permanova_table[, sapply(permanova_table, is.numeric)] <- round(permanova_table[, sapply(permanova_table, is.numeric)], 3)

# Create table plot
table_plot <- tableGrob(permanova_table)

# Save NMDS plot and table
ggsave("NMDS_NEDC.png", 
       grid.arrange(plot_nmds + ggtitle("NEDC"), table_plot, ncol = 1, heights = c(3, 1)), 
       width = 5, height = 7)


# differential analysis ------------
p.plant <- apply(NEDC_final[, 1:225], 2, function(x){kruskal.test(x ~ Plant, data = NEDC_final)$p.value})
#p.adj <- p.adjust(p.plant, method = "fdr")
#calculate p value


NEDC_LC2  <- cbind(NEDC, metadata)

DICA <- NEDC_LC2[which(NEDC_LC2$Plant == "DICA"),1:ncol(NEDC)]
ELRE <- NEDC_LC2[which(NEDC_LC2$Plant == "ERLE"),1:ncol(NEDC)]
LC2_NEDC <- log2(colMeans(DICA)/colMeans(ELRE))
#volcano_NEDC.lm <- data.frame(pvalue = NEDC_results$p_values_plant,
#                               Log2FC = LC2_NEDC)  
volcano_NEDC.kw <- data.frame(pvalue = p.plant,
                               Log2FC = LC2_NEDC) 


Plot_volcano(volcano_NEDC.kw)

open <- NEDC_final[NEDC_final$CanopyStatus == "Open",]
canopy <- NEDC_final[NEDC_final$CanopyStatus == "Canopy",]

p.plant.o <- apply(open[, 1:225], 2, function(x){kruskal.test(x ~ Plant, data = open)$p.value})
p.plant.c <- apply(canopy[, 1:225], 2, function(x){kruskal.test(x ~ Plant, data = open)$p.value})

e.open <- NEDC_LC2[NEDC_LC2$Plant == "ERLE" & NEDC_LC2$CanopyStatus == "Open",]
e.canopy <- NEDC_LC2[NEDC_LC2$Plant == "ERLE" & NEDC_LC2$CanopyStatus == "Canopy",]
d.open <- NEDC_LC2[NEDC_LC2$Plant == "DICA" & NEDC_LC2$CanopyStatus == "Open",]
d.canopy <- NEDC_LC2[NEDC_LC2$Plant == "DICA" & NEDC_LC2$CanopyStatus == "Canopy",]

LC2_open <- log2(colMeans(d.open[,c(1:225)])/colMeans(e.open[,c(1:225)]))
LC2_canopy <- log2(colMeans(d.canopy[,c(1:225)])/colMeans(e.canopy[,c(1:225)]))

volcano_NEDC.open<-data.frame(pvalue = p.plant.o,
                               Log2FC = LC2_open) 
volcano_NEDC.canopy<-data.frame(pvalue = p.plant.c,
                                 Log2FC = LC2_canopy) 

Plot_volcano(volcano_NEDC.canopy)
Plot_volcano(volcano_NEDC.open)

NEDC_LC2 %>%
  ggplot(aes(x = Plant, y = NEDC006)) +
  geom_jitter(aes(color = Plant, shape = CanopyStatus), position = position_jitterdodge(dodge.width = 0.75)) +
  geom_boxplot(aes(color = Plant, fill = Plant, shape = CanopyStatus, alpha = CanopyStatus), outlier.shape = NA) +
  scale_fill_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_color_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_shape_manual(values = c("Canopy" = 16, "Open" = 1)) +
  scale_alpha_manual(values = c("Canopy" = 0.7, "Open" = 0)) +
  labs(title = "Root_asparagine_MALDI_NEDC") +
  xlab(NULL) +
  ylab("MALDI_intensity_average") +
  theme_bw(base_size = 15) +
  theme(legend.position = "none")

NEDC_LC2 %>%
  ggplot(aes(x = Plant, y = NEDC013)) +
  geom_jitter(aes(color = Plant, shape = CanopyStatus), position = position_jitterdodge(dodge.width = 0.75)) +
  geom_boxplot(aes(color = Plant, fill = Plant, shape = CanopyStatus, alpha = CanopyStatus), outlier.shape = NA) +
  scale_fill_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_color_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_shape_manual(values = c("Canopy" = 16, "Open" = 1)) +
  scale_alpha_manual(values = c("Canopy" = 0.7, "Open" = 0)) +
  labs(title = "Root_glutamine_MALDI_NEDC") +
  xlab(NULL) +
  ylab("MALDI_intensity_average") +
  theme_bw(base_size = 15) +
  theme(legend.position = "none")

## DHB---------------------
DHB <- read_csv("DHB.csv") %>%
  #Make the compound IDs the row names
  column_to_rownames('ID') %>%
  #Transpose the matrix to make it so samples are rows and compounds are columns
  t()  %>%
  as.data.frame()

#boxplot(APEBA)

DHB_pareto <- DHB %>%
  log_transform() %>%
  pareto_scale()

#boxplot(APEBA_pareto)

#Combine metadata with metabolite data
DHB_final  <- cbind(DHB_pareto,metadata)

#NMDS plot
set.seed(123)

#plot NMDS
# Compute distance matrix
X <- DHB_final  %>%
  dplyr::select(starts_with('DHB'))
dist_matrix <- vegdist(X, method = "euclidean")

# Perform NMDS
nmds_result <- metaMDS(dist_matrix)

# Add NMDS coordinates to dataframe
DHB_final$NMDS1 <- nmds_result$points[, 1]
DHB_final$NMDS2 <- nmds_result$points[, 2]

# Plot NMDS
plot_nmds <- ggplot(DHB_final, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(shape = CanopyStatus, color = Plant, fill = Plant), size = 3) +
  theme_bw(base_size = 15) +
  scale_color_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_fill_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_shape_manual(values = c("Canopy" = 16, "Open" = 1)) +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  annotate("text", x = Inf, y = -Inf, vjust = 0, hjust = 1, label = paste("Stress =", round(nmds_result$stress, 3)))

plot_nmds

# Perform PERMANOVA
permanova_result <- adonis2(dist_matrix ~ Plant * CanopyStatus, data = DHB_final)

# Round numeric values in PERMANOVA result
permanova_table <- as.data.frame(permanova_result)
permanova_table[, sapply(permanova_table, is.numeric)] <- round(permanova_table[, sapply(permanova_table, is.numeric)], 3)

# Create table plot
table_plot <- tableGrob(permanova_table)

# Save NMDS plot and table
ggsave("NMDS_DHB.png", 
       grid.arrange(plot_nmds + ggtitle("DHB"), table_plot, ncol = 1, heights = c(3, 1)), 
       width = 5, height = 7)


# differential analysis ------------
p.plant <- apply(DHB_final[, 1:122], 2, function(x){kruskal.test(x ~ Plant, data = DHB_final)$p.value})

#calculate p value

DHB_LC2  <- cbind(DHB, metadata)

DICA <- DHB_LC2[which(DHB_LC2$Plant == "DICA"),1:ncol(DHB)]
ELRE <- DHB_LC2[which(DHB_LC2$Plant == "ERLE"),1:ncol(DHB)]
LC2_DHB <- log2(colMeans(DICA)/colMeans(ELRE))
#volcano_DHB.lm <- data.frame(pvalue = DHB_results$p_values_plant,
#                               Log2FC = LC2_NEDC)  
volcano_DHB.kw <- data.frame(pvalue = p.plant,
                              Log2FC = LC2_DHB) 


Plot_volcano(volcano_DHB.kw)

open <- DHB_final[DHB_final$CanopyStatus == "Open",]
canopy <- DHB_final[DHB_final$CanopyStatus == "Canopy",]

p.plant.o <- apply(open[, 1:122], 2, function(x){kruskal.test(x ~ Plant, data = open)$p.value})
p.plant.c <- apply(canopy[, 1:122], 2, function(x){kruskal.test(x ~ Plant, data = open)$p.value})

e.open <- DHB_LC2[DHB_LC2$Plant == "ERLE" & DHB_LC2$CanopyStatus == "Open",]
e.canopy <- DHB_LC2[DHB_LC2$Plant == "ERLE" & DHB_LC2$CanopyStatus == "Canopy",]
d.open <- DHB_LC2[DHB_LC2$Plant == "DICA" & DHB_LC2$CanopyStatus == "Open",]
d.canopy <- DHB_LC2[DHB_LC2$Plant == "DICA" & DHB_LC2$CanopyStatus == "Canopy",]

LC2_open <- log2(colMeans(d.open[,c(1:122)])/colMeans(e.open[,c(1:122)]))
LC2_canopy <- log2(colMeans(d.canopy[,c(1:122)])/colMeans(e.canopy[,c(1:122)]))

volcano_DHB.open<-data.frame(pvalue = p.plant.o,
                              Log2FC = LC2_open) 
volcano_DHB.canopy<-data.frame(pvalue = p.plant.c,
                                Log2FC = LC2_canopy) 

Plot_volcano(volcano_DHB.canopy)
Plot_volcano(volcano_DHB.open)

DHB_LC2 %>%
  ggplot(aes(x = Plant, y = DHB015)) +
  geom_jitter(aes(color = Plant, shape = CanopyStatus), position = position_jitterdodge(dodge.width = 0.75)) +
  geom_boxplot(aes(color = Plant, fill = Plant, shape = CanopyStatus, alpha = CanopyStatus), outlier.shape = NA) +
  scale_fill_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_color_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_shape_manual(values = c("Canopy" = 16, "Open" = 1)) +
  scale_alpha_manual(values = c("Canopy" = 0.7, "Open" = 0)) +
  labs(title = "Root_asparagine_MALDI_DHB") +
  xlab(NULL) +
  ylab("MALDI_intensity_average") +
  theme_bw(base_size = 15) +
  theme(legend.position = "none")

DHB_LC2 %>%
  ggplot(aes(x = Plant, y = DHB015)) +
  geom_jitter(aes(color = Plant, shape = CanopyStatus), position = position_jitterdodge(dodge.width = 0.75)) +
  geom_boxplot(aes(color = Plant, fill = Plant, shape = CanopyStatus, alpha = CanopyStatus), outlier.shape = NA) +
  scale_fill_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_color_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_shape_manual(values = c("Canopy" = 16, "Open" = 1)) +
  scale_alpha_manual(values = c("Canopy" = 0.7, "Open" = 0)) +
  labs(title = "Root_asparagine_MALDI_DHB") +
  xlab(NULL) +
  ylab("MALDI_intensity_average") +
  theme_bw(base_size = 15) +
  theme(legend.position = "none")

#View differentially expressed metabolites
diff_expressed <- volcano_APEBA.kw[volcano_APEBA.kw$diffexpressed != "NO", ]

print(diff_expressed)

