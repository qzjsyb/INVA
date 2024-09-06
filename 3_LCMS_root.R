# Environment -------------------------------------------------------------
library(tidyverse)
library(mixOmics)
library(gridExtra)
library(ggrepel)
library(vegan)
library(forcats)
library(ggbreak)

setwd('~/INVA/INVA')
metadata <- read.table("metadata.txt", header = TRUE, sep = "\t", row.names = 1)
metadata$block <- as.factor(metadata$block)
Root_compounds <- read.csv("Root_compound_link.csv", header = TRUE, row.names = 1)
Root_compounds <- rownames_to_column(Root_compounds, var = "compoundID")

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
#plot NMDS
Plot_nmds <- function(data) {
  # Extract dataframe name as filename
  filename <- deparse(substitute(data))
  filename <- sub("_final","", filename)
  
  # Compute distance matrix
  X <- data  %>%
    dplyr::select(starts_with('FT_'))
  
  dist_matrix <- vegdist(X, method = "euclidean")
  
  # Save distance matrix for future distance
  saveRDS(dist_matrix, file = paste0(filename, "_dis.rds"))
  
  # Perform NMDS
  nmds_result <- metaMDS(dist_matrix)
  
  # Add NMDS coordinates to dataframe
  data$NMDS1 <- nmds_result$points[, 1]
  data$NMDS2 <- nmds_result$points[, 2]
  
  # Plot NMDS
  plot_nmds <- ggplot(data, aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(shape = canopy, color = plant, fill = plant), size = 3) +
    theme_bw(base_size = 15) +
    scale_color_manual(values = c("ERLE" = "#50124C", "DICA" = "#4DAF4A")) +
    scale_fill_manual(values = c("ERLE" = "#50124C", "DICA" = "#4DAF4A")) +
    scale_shape_manual(values = c("close" = 16, "open" = 1)) +
    theme(legend.position = "none") +
    labs(x = NULL, y = NULL) +
    annotate("text", x = Inf, y = -Inf, vjust = 0, hjust = 1, label = paste("Stress =", round(nmds_result$stress, 3)))
  
  # Perform PERMANOVA
  permanova_result <- adonis2(dist_matrix ~ plant * canopy, data = data)
  
  # Round numeric values in PERMANOVA result
  permanova_table <- as.data.frame(permanova_result)
  permanova_table[, sapply(permanova_table, is.numeric)] <- round(permanova_table[, sapply(permanova_table, is.numeric)], 3)
  
  # Create table plot
  table_plot <- tableGrob(permanova_table)
  
  # Save NMDS plot and table
  ggsave(paste0("NMDS_", filename, ".pdf"), plot_nmds, width = 4, height = 4)
}
  
#    ggsave(paste0("NMDS_", filename, ".pdf"), 
 #        grid.arrange(plot_nmds + ggtitle(filename), table_plot, ncol = 1, heights = c(3, 1)), 
  #       width = 5, height = 7)
#}

Plot_nmds_oc <- function(data) {
  filename <- deparse(substitute(data))
  filename <- sub("_final","", filename) 
  # Compute distance matrix
  X <- data  %>%
    dplyr::select(starts_with('FT_'))
  
  dist_matrix <- vegdist(X, method = "euclidean")
  
  assign(paste0("NMDS_", filename, "_dist"), dist_matrix, envir = .GlobalEnv)
  
  # Perform NMDS
  nmds_result <- metaMDS(dist_matrix)
  
  # Add NMDS coordinates to dataframe
  data$NMDS1 <- nmds_result$points[, 1]
  data$NMDS2 <- nmds_result$points[, 2]

plot_nmds_oc <- ggplot(data, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes( color = canopy, fill = canopy), size = 3) +
  theme_bw(base_size = 15) +
  scale_color_manual(values = c("open" = "#FFAF00", "close" = "#00B6BD")) +
  scale_fill_manual(values = c("open" = "#FFAF00", "close" = "#00B6BD")) +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  annotate("text", x = Inf, y = -Inf, vjust = 0, hjust = 1, label = paste("Stress =", round(nmds_result$stress, 3)))

# Perform PERMANOVA
permanova_result <- adonis2(dist_matrix ~ canopy, data = data)

# Round numeric values in PERMANOVA result
permanova_table <- as.data.frame(permanova_result)
permanova_table[, sapply(permanova_table, is.numeric)] <- round(permanova_table[, sapply(permanova_table, is.numeric)], 3)

# Create table plot
table_plot <- tableGrob(permanova_table)

# Save NMDS plot and table
ggsave(paste0("NMDS_", filename, ".pdf"), 
       grid.arrange(plot_nmds_oc + ggtitle(filename), table_plot, ncol = 1, heights = c(3, 1)), 
       width = 5, height = 7)
}



#plot volcano between DICA and ERLE
Plot_volcano <- function(dataframe) {
  # Create a copy of the input dataframe
  volcano_data <- dataframe
  filename <- deparse(substitute(dataframe))
  
  # Perform data transformations and calculations
  volcano_data$diffexpressed <- "NO"
  volcano_data$diffexpressed[volcano_data$Log2FC > 1.5 & volcano_data$pvalue < 0.05] <- "DICA"
  volcano_data$diffexpressed[volcano_data$Log2FC < -1.5 & volcano_data$pvalue < 0.05] <- "ERLE"
  volcano_data$delabel <- NA
  volcano_data$compoundID <- row.names(volcano_data)
  volcano_data <- left_join(volcano_data, Root_compounds, by = "compoundID")
  
  # Label the data
  volcano_data[grepl("Confirmed ID", volcano_data$Tags),]$delabel <- volcano_data[grepl("Confirmed ID", volcano_data$Tags),]$Name
  volcano_data[grepl("Potential ID", volcano_data$Tags),]$delabel <- volcano_data[grepl("Potential ID", volcano_data$Tags),]$Formula
  volcano_data[grepl("only", volcano_data$Tags),]$delabel <- paste("?", volcano_data[grepl("only", volcano_data$Tags),]$Name, sep = "")
  volcano_data$delabel[volcano_data$diffexpressed == "NO"] <- NA
  
  # Save the processed dataframe with a name corresponding to the input dataframe
  assign(filename, volcano_data, envir = .GlobalEnv)
  
  # Create volcano plot
  plot_volcano <- ggplot(data = volcano_data, aes(x = Log2FC, y = -log10(pvalue), col = diffexpressed, label = delabel)) + 
    geom_vline(xintercept=c(-1.5, 1.5), col="dimgrey") +
    geom_hline(yintercept=-log10(0.05), col="dimgrey") +
    geom_point() + 
    theme_minimal() +
    geom_text_repel(max.overlaps = Inf) +
    scale_color_manual(values = c("ERLE" = "#50124C", "DICA" = "#4AAF49", "NO" = "grey"))
  
  return(plot_volcano)
}
#plot volcano between open and close
Plot_volcano_oc <- function(dataframe) {
  # Create a copy of the input dataframe
  volcano_data <- dataframe
  filename <- deparse(substitute(dataframe))
  
  # Perform data transformations and calculations
  volcano_data$diffexpressed <- "NO"
  volcano_data$diffexpressed[volcano_data$Log2FC > 1.5 & volcano_data$pvalue < 0.05] <- "CLOSE"
  volcano_data$diffexpressed[volcano_data$Log2FC < -1.5 & volcano_data$pvalue < 0.05] <- "OPEN"
  volcano_data$delabel <- NA
  volcano_data$compoundID <- row.names(volcano_data)
  volcano_data <- left_join(volcano_data, Root_compounds, by = "compoundID")
  
  # Label the data
  volcano_data[grepl("Confirmed ID", volcano_data$Tags),]$delabel <- volcano_data[grepl("Confirmed ID", volcano_data$Tags),]$Name
  volcano_data[grepl("Potential ID", volcano_data$Tags),]$delabel <- volcano_data[grepl("Potential ID", volcano_data$Tags),]$Formula
  volcano_data[grepl("only", volcano_data$Tags),]$delabel <- paste("?", volcano_data[grepl("only", volcano_data$Tags),]$Name, sep = "")
  volcano_data$delabel[volcano_data$diffexpressed == "NO"] <- NA
  
  # Save the processed dataframe with a name corresponding to the input dataframe
  assign(filename, volcano_data, envir = .GlobalEnv)
  
  # Create volcano plot
  plot_volcano <- ggplot(data = volcano_data, aes(x = Log2FC, y = -log10(pvalue), col = diffexpressed, label = delabel)) + 
    geom_vline(xintercept=c(-1.5, 1.5), col="dimgrey") +
    geom_hline(yintercept=-log10(0.05), col="dimgrey") +
    geom_point() + 
    theme_minimal() +
    geom_text_repel(max.overlaps = Inf) +
    scale_color_manual(values=c("OPEN" = "#FBB413", "CLOSE" = "#06B4BC", "NO" = "grey"))
  
  return(plot_volcano)
}

## Rp positive---------------------
rppos_raw <- read_csv("Root_LCMS_RPPOS.csv") %>%
  #Make the compound IDs the row names
  column_to_rownames('Compound ID') %>%
  #Transpose the matrix to make it so samples are rows and compounds are columns
  t()  %>%
  as.data.frame()

#R27 has low read on everything that makes the analysis crazy.
rppos_raw <- rppos_raw[-which(rownames(rppos_raw) == "R27_RP_Pos"), ]

rppos_pareto <- rppos_raw %>%
  log_transform() %>%
  pareto_scale()

#Combine metadata with metabolite data
Root_rppos_final  <- rppos_pareto %>%
  rownames_to_column('SampleID_root') %>%
  mutate(SampleID_root = sub("_RP_Pos$", "", SampleID_root)) %>% 
  left_join(metadata)%>%
  column_to_rownames("SampleID_root")



#Plot NMDS
Plot_nmds(Root_rppos_final)

Root_rppos_final_ERLE <- Root_rppos_final[Root_rppos_final$plant== "ERLE",]
Root_rppos_final_DICA <- Root_rppos_final[Root_rppos_final$plant== "DICA",]
Plot_nmds_oc(Root_rppos_final_ERLE)
Plot_nmds_oc(Root_rppos_final_DICA)



## Rp negative (select and replace rppos to rpneg; RP_Pos to RP_Neg)---------------------
rpneg_raw <- read_csv("Root_LCMS_RPNEG.csv") %>%
    #Make the compound IDs the row names
    column_to_rownames('Compound ID') %>%
    #Transnege the matrix to make it so samples are rows and compounds are columns
    t()  %>%
    as.data.frame()
  
### here is for testing different transforming methods.
boxplot(rpneg_raw) # big variaation.



rpneg_pareto <- rpneg_raw %>%
  log_transform() %>%
  pareto_scale()

#Combine metadata with metabolite data
Root_rpneg_final  <- rpneg_pareto %>%
  rownames_to_column('SampleID_root') %>%
  mutate(SampleID_root = sub("_RP_Neg$", "", SampleID_root)) %>% 
  right_join(metadata)%>%
  column_to_rownames("SampleID_root")

Plot_nmds(Root_rpneg_final)

Root_rpneg_final_ERLE <- Root_rpneg_final[Root_rpneg_final$plant== "ERLE",]
Root_rpneg_final_DICA <- Root_rpneg_final[Root_rpneg_final$plant== "DICA",]
Plot_nmds_oc(Root_rpneg_final_ERLE)
Plot_nmds_oc(Root_rpneg_final_DICA)

## HILIC positive-----------------
  hpos_raw <- read_csv("Root_LCMS_HPOS.csv") %>%
    #Make the compound IDs the row names
    column_to_rownames('Compound ID') %>%
    #Transnege the matrix to make it so samples are rows and compounds are columns
    t()  %>%
    as.data.frame()
  
### here is for testing different transforming methods.
boxplot(hpos_raw) # big variaation.

hpos_pareto <- hpos_raw %>%
  log_transform() %>%
  pareto_scale()

#Combine metadata with metabolite data
Root_hpos_final  <- hpos_pareto %>%
  rownames_to_column('SampleID_root') %>%
  mutate(SampleID_root = sub("_HILIC_Pos$", "", SampleID_root)) %>% 
  right_join(metadata)%>%
  column_to_rownames("SampleID_root")

Plot_nmds(Root_hpos_final)


  
Root_hpos_final_ERLE <- Root_hpos_final[Root_hpos_final$plant== "ERLE",]
Root_hpos_final_DICA <- Root_hpos_final[Root_hpos_final$plant== "DICA",]
Plot_nmds_oc(Root_hpos_final_ERLE)
Plot_nmds_oc(Root_hpos_final_DICA)
# HILIC negative-----------------
  hneg_raw <- read_csv("Root_LCMS_HNEG.csv") %>%
    #Make the compound IDs the row names
    column_to_rownames('Compound ID') %>%
    #Transnege the matrix to make it so samples are rows and compounds are columns
    t()  %>%
    as.data.frame()
  
  ### here is for testing different transforming methods.
  boxplot(hneg_raw) # big variaation.
  
  hneg_pareto <- hneg_raw %>%
    log_transform() %>%
    pareto_scale()
  
  #Combine metadata with metabolite data
  Root_hneg_final  <- hneg_pareto %>%
    rownames_to_column('SampleID_root') %>%
    mutate(SampleID_root = sub("_HILIC_Neg$", "", SampleID_root)) %>% 
    left_join(metadata)%>%
    column_to_rownames("SampleID_root")
 
Plot_nmds(Root_hneg_final)


Root_hneg_final_ERLE <- Root_hneg_final[Root_hneg_final$plant== "ERLE",]
Root_hneg_final_DICA <- Root_hneg_final[Root_hneg_final$plant== "DICA",]
Plot_nmds_oc(Root_hneg_final_ERLE)
Plot_nmds_oc(Root_hneg_final_DICA)

# differential analysis ------------
##rppos########
p.plant <- apply(Root_rppos_final[,1:ncol(rppos_raw)], 2, function(x){kruskal.test(x~plant, data = Root_rppos_final)$p.value})  

  Root_rppos_LC2  <- rppos_raw %>%
    rownames_to_column('SampleID_root') %>%
    mutate(SampleID_root = sub("_RP_Pos$", "", SampleID_root)) %>% 
    left_join(metadata) %>%
    column_to_rownames('SampleID_root')
  
  DICA <- Root_rppos_LC2[which(Root_rppos_LC2$plant == "DICA"),1:ncol(rppos_raw)]
  ELRE <- Root_rppos_LC2[which(Root_rppos_LC2$plant == "ERLE"),1:ncol(rppos_raw)]
  LC2_rppos <- log2(colMeans(DICA)/colMeans(ELRE))
  volcano_rppos <- data.frame(pvalue = p.plant,
                              Log2FC = LC2_rppos)  
  
  Plot_Root_volcano_rppos<-Plot_volcano(volcano_rppos)
  
##rpneg-----------
  p.plant <- apply(Root_rpneg_final[,1:ncol(rpneg_raw)], 2, function(x){kruskal.test(x~plant, data = Root_rpneg_final)$p.value})  
  
  Root_rpneg_LC2  <- rpneg_raw %>%
    rownames_to_column('SampleID_root') %>%
    mutate(SampleID_root = sub("_RP_Neg$", "", SampleID_root)) %>% 
    right_join(metadata)%>%
    column_to_rownames('SampleID_root')
  
  DICA <- Root_rpneg_LC2[which(Root_rpneg_LC2$plant == "DICA"),1:ncol(rpneg_raw)]
  ELRE <- Root_rpneg_LC2[which(Root_rpneg_LC2$plant == "ERLE"),1:ncol(rpneg_raw)]
  LC2_rpneg <- log2(colMeans(DICA)/colMeans(ELRE))
  volcano_rpneg <- data.frame(pvalue = p.plant,
                              Log2FC = LC2_rpneg)  

  Plot_Root_volcano_rpneg <- Plot_volcano(volcano_rpneg)
  
  
  
##hpos-----------
  p.plant <- apply(Root_hpos_final[,1:ncol(hpos_raw)], 2, function(x){kruskal.test(x~plant, data = Root_hpos_final)$p.value})  
  
  Root_hpos_LC2  <- hpos_raw %>%
    rownames_to_column('SampleID_root') %>%
    mutate(SampleID_root = sub("_HILIC_Pos$", "", SampleID_root)) %>% 
    right_join(metadata)%>%
    column_to_rownames('SampleID_root')
  
  DICA <- Root_hpos_LC2[which(Root_hpos_LC2$plant == "DICA"),1:ncol(hpos_raw)]
  ELRE <- Root_hpos_LC2[which(Root_hpos_LC2$plant == "ERLE"),1:ncol(hpos_raw)]
  LC2_hpos <- log2(colMeans(DICA)/colMeans(ELRE))
  volcano_hpos <- data.frame(pvalue = p.plant,
                             Log2FC = LC2_hpos)  
  
  Plot_Root_volcano_hpos <- Plot_volcano(volcano_hpos)
  
  
  
  ##hneg-----------
  p.plant <- apply(Root_hneg_final[,1:ncol(hneg_raw)], 2, function(x){kruskal.test(x~plant, data = Root_hneg_final)$p.value})  
  
  Root_hneg_LC2  <- hneg_raw %>%
    rownames_to_column('SampleID_root') %>%
    mutate(SampleID_root = sub("_HILIC_Neg$", "", SampleID_root)) %>% 
    right_join(metadata)%>%
    column_to_rownames('SampleID_root')
  
  DICA <- Root_hneg_LC2[which(Root_hneg_LC2$plant == "DICA"),1:ncol(hneg_raw)]
  ELRE <- Root_hneg_LC2[which(Root_hneg_LC2$plant == "ERLE"),1:ncol(hneg_raw)]
  LC2_hneg <- log2(colMeans(DICA)/colMeans(ELRE))
  volcano_hneg <- data.frame(pvalue = p.plant,
                             Log2FC = LC2_hneg)  
  
  
  Plot_Root_volcano_hneg <- Plot_volcano(volcano_hneg)
  
  
  ########### output ----------
  ggsave("Root_volcano_hneg.pdf", Plot_Root_volcano_hneg + ggtitle("hneg plant") + theme(legend.position = "none"), width = 6, height = 4)
  ggsave("Root_volcano_hpos.pdf", Plot_Root_volcano_hpos + ggtitle("hpos plant") + theme(legend.position = "none"), width = 6, height = 4)
  ggsave("Root_volcano_rppos.pdf", Plot_Root_volcano_rppos + ggtitle("rppos plant") + theme(legend.position = "none"), width = 6, height = 4)
  ggsave("Root_volcano_rpneg.pdf", Plot_Root_volcano_rpneg + ggtitle("rpneg plant") + theme(legend.position = "none"), width = 6, height = 4)
  
  
  Diff_Root_plant <- rbind(volcano_rppos, volcano_hpos, volcano_rpneg, volcano_hneg)
  # Remove duplicates in "Name" column, prioritizing rows with non-NA values in "delabel" column
 
  write.table(Diff_Root_plant, file = "Diff_Root_plant.csv", sep = ",", row.names = FALSE)
  

  ### now compare ERLE open vs close and DICA open vs close -----
  #rppos-------
  e <- Root_rppos_final[Root_rppos_final$plant == "ERLE", ]
  d <- Root_rppos_final[Root_rppos_final$plant == "DICA", ]
  
  p_eoc <- apply(e[,c(1:ncol(rppos_raw))],2, function(x){kruskal.test(x~canopy, data = e)$p.value})
  p_doc <- apply(d[,c(1:ncol(rppos_raw))],2, function(x){kruskal.test(x~canopy, data = d)$p.value})
  
  LC2.eo <- Root_rppos_LC2[Root_rppos_LC2$plant == "ERLE" & Root_rppos_LC2$canopy == "open", 1:ncol(rppos_raw)]
  LC2.ec <- Root_rppos_LC2[Root_rppos_LC2$plant == "ERLE" & Root_rppos_LC2$canopy == "close", 1:ncol(rppos_raw)]
  LC2.do <- Root_rppos_LC2[Root_rppos_LC2$plant == "DICA" & Root_rppos_LC2$canopy == "open", 1:ncol(rppos_raw)]
  LC2.dc <- Root_rppos_LC2[Root_rppos_LC2$plant == "DICA" & Root_rppos_LC2$canopy == "close", 1:ncol(rppos_raw)]
  
  LC2.eoc <- log2(colMeans(LC2.ec)/colMeans(LC2.eo))
  LC2.doc <- log2(colMeans(LC2.dc)/colMeans(LC2.do))
  
  volcano_rppos_eoc <- data.frame(pvalue = p_eoc,
                                  Log2FC = LC2.eoc) 
  volcano_rppos_doc <- data.frame(pvalue = p_doc,
                                  Log2FC = LC2.doc) 
  
  Plot_Root_volcano_rppos_eoc <- Plot_volcano_oc(volcano_rppos_eoc)
  Plot_Root_volcano_rppos_doc <- Plot_volcano_oc(volcano_rppos_doc)
  ggsave("Plot_Root_volcano_rppos_eoc.pdf", Plot_Root_volcano_rppos_eoc + ggtitle("rppos ERLE by canopy") + theme(legend.position = "none"), width = 6, height = 4)
  ggsave("Plot_Root_volcano_rppos_doc.pdf", Plot_Root_volcano_rppos_doc + ggtitle("rppos DICA by canopy") + theme(legend.position = "none"), width = 6, height = 4)
  
  
  #rpneg-------
  e <- Root_rpneg_final[Root_rpneg_final$plant == "ERLE", ]
  d <- Root_rpneg_final[Root_rpneg_final$plant == "DICA", ]
  
  p_eoc <- apply(e[,c(1:ncol(rpneg_raw))],2, function(x){kruskal.test(x~canopy, data = e)$p.value})
  p_doc <- apply(d[,c(1:ncol(rpneg_raw))],2, function(x){kruskal.test(x~canopy, data = d)$p.value})
  
  
  LC2.eo <- Root_rpneg_LC2[Root_rpneg_LC2$plant == "ERLE" & Root_rpneg_LC2$canopy == "open", 1:ncol(rpneg_raw)]
  LC2.ec <- Root_rpneg_LC2[Root_rpneg_LC2$plant == "ERLE" & Root_rpneg_LC2$canopy == "close", 1:ncol(rpneg_raw)]
  LC2.do <- Root_rpneg_LC2[Root_rpneg_LC2$plant == "DICA" & Root_rpneg_LC2$canopy == "open", 1:ncol(rpneg_raw)]
  LC2.dc <- Root_rpneg_LC2[Root_rpneg_LC2$plant == "DICA" & Root_rpneg_LC2$canopy == "close", 1:ncol(rpneg_raw)]
  
  LC2.eoc <- log2(colMeans(LC2.ec)/colMeans(LC2.eo))
  LC2.doc <- log2(colMeans(LC2.dc)/colMeans(LC2.do))
  
  volcano_rpneg_eoc <- data.frame(pvalue = p_eoc,
                                  Log2FC = LC2.eoc) 
  volcano_rpneg_doc <- data.frame(pvalue = p_doc,
                                  Log2FC = LC2.doc) 
  
  Plot_Root_volcano_rpneg_eoc <- Plot_volcano_oc(volcano_rpneg_eoc)
  Plot_Root_volcano_rpneg_doc <- Plot_volcano_oc(volcano_rpneg_doc)
  ggsave("Plot_Root_volcano_rpneg_eoc.pdf", Plot_Root_volcano_rpneg_eoc + ggtitle("rpneg ERLE by canopy") + theme(legend.position = "none"), width = 6, height = 4)
  ggsave("Plot_Root_volcano_rpneg_doc.pdf", Plot_Root_volcano_rpneg_doc + ggtitle("rpneg DICA by canopy") + theme(legend.position = "none"), width = 6, height = 4)
  
  
  #hpos-------
  e <- Root_hpos_final[Root_hpos_final$plant == "ERLE", ]
  d <- Root_hpos_final[Root_hpos_final$plant == "DICA", ]
  
  p_eoc <- apply(e[,c(1:ncol(hpos_raw))],2, function(x){kruskal.test(x~canopy, data = e)$p.value})
  p_doc <- apply(d[,c(1:ncol(hpos_raw))],2, function(x){kruskal.test(x~canopy, data = d)$p.value})
  
  LC2.eo <- Root_hpos_LC2[Root_hpos_LC2$plant == "ERLE" & Root_hpos_LC2$canopy == "open", 1:ncol(hpos_raw)]
  LC2.ec <- Root_hpos_LC2[Root_hpos_LC2$plant == "ERLE" & Root_hpos_LC2$canopy == "close", 1:ncol(hpos_raw)]
  LC2.do <- Root_hpos_LC2[Root_hpos_LC2$plant == "DICA" & Root_hpos_LC2$canopy == "open", 1:ncol(hpos_raw)]
  LC2.dc <- Root_hpos_LC2[Root_hpos_LC2$plant == "DICA" & Root_hpos_LC2$canopy == "close", 1:ncol(hpos_raw)]
  
  LC2.eoc <- log2(colMeans(LC2.ec)/colMeans(LC2.eo))
  LC2.doc <- log2(colMeans(LC2.dc)/colMeans(LC2.do))
  
  volcano_hpos_eoc <- data.frame(pvalue = p_eoc,
                                 Log2FC = LC2.eoc) 
  volcano_hpos_doc <- data.frame(pvalue = p_doc,
                                 Log2FC = LC2.doc) 
  
  Plot_Root_volcano_hpos_eoc <- Plot_volcano_oc(volcano_hpos_eoc)
  Plot_Root_volcano_hpos_doc <- Plot_volcano_oc(volcano_hpos_doc)
  ggsave("Plot_Root_volcano_hpos_eoc.pdf", Plot_Root_volcano_hpos_eoc + ggtitle("hpos ERLE by canopy") + theme(legend.position = "none"), width = 6, height = 4)
  ggsave("Plot_Root_volcano_hpos_doc.pdf", Plot_Root_volcano_hpos_doc + ggtitle("hpos DICA by canopy") + theme(legend.position = "none"), width = 6, height = 4)
  
  
  
  #hneg-------
  e <- Root_hneg_final[Root_hneg_final$plant == "ERLE", ]
  d <- Root_hneg_final[Root_hneg_final$plant == "DICA", ]
  
  p_eoc <- apply(e[,c(1:ncol(hneg_raw))],2, function(x){kruskal.test(x~canopy, data = e)$p.value})
  p_doc <- apply(d[,c(1:ncol(hneg_raw))],2, function(x){kruskal.test(x~canopy, data = d)$p.value})
  
  
  LC2.eo <- Root_hneg_LC2[Root_hneg_LC2$plant == "ERLE" & Root_hneg_LC2$canopy == "open", 1:ncol(hneg_raw)]
  LC2.ec <- Root_hneg_LC2[Root_hneg_LC2$plant == "ERLE" & Root_hneg_LC2$canopy == "close", 1:ncol(hneg_raw)]
  LC2.do <- Root_hneg_LC2[Root_hneg_LC2$plant == "DICA" & Root_hneg_LC2$canopy == "open", 1:ncol(hneg_raw)]
  LC2.dc <- Root_hneg_LC2[Root_hneg_LC2$plant == "DICA" & Root_hneg_LC2$canopy == "close", 1:ncol(hneg_raw)]
  
  LC2.eoc <- log2(colMeans(LC2.ec)/colMeans(LC2.eo))
  LC2.doc <- log2(colMeans(LC2.dc)/colMeans(LC2.do))
  
  volcano_hneg_eoc <- data.frame(pvalue = p_eoc,
                                 Log2FC = LC2.eoc) 
  volcano_hneg_doc <- data.frame(pvalue = p_doc,
                                 Log2FC = LC2.doc) 
  
  Plot_Root_volcano_hneg_eoc <- Plot_volcano_oc(volcano_hneg_eoc)
  Plot_Root_volcano_hneg_doc <- Plot_volcano_oc(volcano_hneg_doc)
  ggsave("Plot_Root_volcano_hneg_eoc.pdf", Plot_Root_volcano_hneg_eoc + ggtitle("hneg ERLE by canopy") + theme(legend.position = "none"), width = 6, height = 4)
  ggsave("Plot_Root_volcano_hneg_doc.pdf", Plot_Root_volcano_hneg_doc + ggtitle("hneg DICA by canopy") + theme(legend.position = "none"), width = 6, height = 4)
  
  
  # make table -------
  
  Diff_Root_ERLE_canopy <- rbind(volcano_rppos_eoc, volcano_hpos_eoc, volcano_rpneg_eoc, volcano_hneg_eoc)
  

  # Remove duplicates in "Name" column, prioritizing rows with non-NA values in "delabel" column
  
  Diff_Root_DICA_canopy <- rbind(volcano_rppos_doc, volcano_hpos_doc, volcano_rpneg_doc, volcano_hneg_doc)
  
  write.table(Diff_Root_ERLE_canopy, file = "Diff_Root_ERLE_canopy.csv", sep = ",", row.names = FALSE)
  write.table(Diff_Root_DICA_canopy, file = "Diff_Root_DICA_canopy.csv", sep = ",", row.names = FALSE)
  

  
  ### now compare open ERLE vs DICA and close ERLE vs DICA -----
  #rppos-------
  o <- Root_rppos_final[Root_rppos_final$canopy == "open", ]
  c <- Root_rppos_final[Root_rppos_final$canopy == "close", ]
  
  p_oed <- apply(o[,c(1:ncol(rppos_raw))],2, function(x){kruskal.test(x~plant, data = o)$p.value})
  p_ced <- apply(c[,c(1:ncol(rppos_raw))],2, function(x){kruskal.test(x~plant, data = c)$p.value})
  
  LC2.eo <- Root_rppos_LC2[Root_rppos_LC2$plant == "ERLE" & Root_rppos_LC2$canopy == "open", 1:ncol(rppos_raw)]
  LC2.ec <- Root_rppos_LC2[Root_rppos_LC2$plant == "ERLE" & Root_rppos_LC2$canopy == "close", 1:ncol(rppos_raw)]
  LC2.do <- Root_rppos_LC2[Root_rppos_LC2$plant == "DICA" & Root_rppos_LC2$canopy == "open", 1:ncol(rppos_raw)]
  LC2.dc <- Root_rppos_LC2[Root_rppos_LC2$plant == "DICA" & Root_rppos_LC2$canopy == "close", 1:ncol(rppos_raw)]
  
  LC2.oed <- log2(colMeans(LC2.do)/colMeans(LC2.eo))
  LC2.ced <- log2(colMeans(LC2.dc)/colMeans(LC2.ec))
  
  volcano_rppos_oed <- data.frame(pvalue = p_oed,
                                  Log2FC = LC2.oed) 
  volcano_rppos_ced <- data.frame(pvalue = p_ced,
                                  Log2FC = LC2.ced) 
  
  Plot_Root_volcano_rppos_oed <- Plot_volcano(volcano_rppos_oed)
  Plot_Root_volcano_rppos_ced <- Plot_volcano(volcano_rppos_ced)
  ggsave("Plot_Root_volcano_rppos_oed.pdf", Plot_Root_volcano_rppos_oed + ggtitle("rppos open ERLE vs DICA") + theme(legend.position = "none"), width = 6, height = 4)
  ggsave("Plot_Root_volcano_rppos_ced.pdf", Plot_Root_volcano_rppos_ced + ggtitle("rppos close ERLE vs DICA") + theme(legend.position = "none"), width = 6, height = 4)
  
  
  #rpneg-------
  o <- Root_rpneg_final[Root_rpneg_final$canopy == "open", ]
  c <- Root_rpneg_final[Root_rpneg_final$canopy == "close", ]
  
  p_oed <- apply(o[,c(1:ncol(rpneg_raw))],2, function(x){kruskal.test(x~plant, data = o)$p.value})
  p_ced <- apply(c[,c(1:ncol(rpneg_raw))],2, function(x){kruskal.test(x~plant, data = c)$p.value})
  
  LC2.eo <- Root_rpneg_LC2[Root_rpneg_LC2$plant == "ERLE" & Root_rpneg_LC2$canopy == "open", 1:ncol(rpneg_raw)]
  LC2.ec <- Root_rpneg_LC2[Root_rpneg_LC2$plant == "ERLE" & Root_rpneg_LC2$canopy == "close", 1:ncol(rpneg_raw)]
  LC2.do <- Root_rpneg_LC2[Root_rpneg_LC2$plant == "DICA" & Root_rpneg_LC2$canopy == "open", 1:ncol(rpneg_raw)]
  LC2.dc <- Root_rpneg_LC2[Root_rpneg_LC2$plant == "DICA" & Root_rpneg_LC2$canopy == "close", 1:ncol(rpneg_raw)]
  
  LC2.oed <- log2(colMeans(LC2.do)/colMeans(LC2.eo))
  LC2.ced <- log2(colMeans(LC2.dc)/colMeans(LC2.ec))
  
  volcano_rpneg_oed <- data.frame(pvalue = p_oed,
                                  Log2FC = LC2.oed) 
  volcano_rpneg_ced <- data.frame(pvalue = p_ced,
                                  Log2FC = LC2.ced) 
  
  Plot_Root_volcano_rpneg_oed <- Plot_volcano(volcano_rpneg_oed)
  Plot_Root_volcano_rpneg_ced <- Plot_volcano(volcano_rpneg_ced)
  ggsave("Plot_Root_volcano_rpneg_oed.pdf", Plot_Root_volcano_rpneg_oed + ggtitle("rpneg open ERLE vs DICA") + theme(legend.position = "none"), width = 6, height = 4)
  ggsave("Plot_Root_volcano_rpneg_ced.pdf", Plot_Root_volcano_rpneg_ced + ggtitle("rpneg close ERLE vs DICA") + theme(legend.position = "none"), width = 6, height = 4)
  
  #hpos-------
  o <- Root_hpos_final[Root_hpos_final$canopy == "open", ]
  c <- Root_hpos_final[Root_hpos_final$canopy == "close", ]
  
  p_oed <- apply(o[,c(1:ncol(hpos_raw))],2, function(x){kruskal.test(x~plant, data = o)$p.value})
  p_ced <- apply(c[,c(1:ncol(hpos_raw))],2, function(x){kruskal.test(x~plant, data = c)$p.value})
  
  LC2.eo <- Root_hpos_LC2[Root_hpos_LC2$plant == "ERLE" & Root_hpos_LC2$canopy == "open", 1:ncol(hpos_raw)]
  LC2.ec <- Root_hpos_LC2[Root_hpos_LC2$plant == "ERLE" & Root_hpos_LC2$canopy == "close", 1:ncol(hpos_raw)]
  LC2.do <- Root_hpos_LC2[Root_hpos_LC2$plant == "DICA" & Root_hpos_LC2$canopy == "open", 1:ncol(hpos_raw)]
  LC2.dc <- Root_hpos_LC2[Root_hpos_LC2$plant == "DICA" & Root_hpos_LC2$canopy == "close", 1:ncol(hpos_raw)]
  
  LC2.oed <- log2(colMeans(LC2.do)/colMeans(LC2.eo))
  LC2.ced <- log2(colMeans(LC2.dc)/colMeans(LC2.ec))
  
  volcano_hpos_oed <- data.frame(pvalue = p_oed,
                                  Log2FC = LC2.oed) 
  volcano_hpos_ced <- data.frame(pvalue = p_ced,
                                  Log2FC = LC2.ced) 
  
  Plot_Root_volcano_hpos_oed <- Plot_volcano(volcano_hpos_oed)
  Plot_Root_volcano_hpos_ced <- Plot_volcano(volcano_hpos_ced)
  ggsave("Plot_Root_volcano_hpos_oed.pdf", Plot_Root_volcano_hpos_oed + ggtitle("hpos open ERLE vs DICA") + theme(legend.position = "none"), width = 6, height = 4)
  ggsave("Plot_Root_volcano_hpos_ced.pdf", Plot_Root_volcano_hpos_ced + ggtitle("hpos close ERLE vs DICA") + theme(legend.position = "none"), width = 6, height = 4)
  
  
  #hneg-------
  o <- Root_hneg_final[Root_hneg_final$canopy == "open", ]
  c <- Root_hneg_final[Root_hneg_final$canopy == "close", ]
  
  p_oed <- apply(o[,c(1:ncol(hneg_raw))],2, function(x){kruskal.test(x~plant, data = o)$p.value})
  p_ced <- apply(c[,c(1:ncol(hneg_raw))],2, function(x){kruskal.test(x~plant, data = c)$p.value})
  
  LC2.eo <- Root_hneg_LC2[Root_hneg_LC2$plant == "ERLE" & Root_hneg_LC2$canopy == "open", 1:ncol(hneg_raw)]
  LC2.ec <- Root_hneg_LC2[Root_hneg_LC2$plant == "ERLE" & Root_hneg_LC2$canopy == "close", 1:ncol(hneg_raw)]
  LC2.do <- Root_hneg_LC2[Root_hneg_LC2$plant == "DICA" & Root_hneg_LC2$canopy == "open", 1:ncol(hneg_raw)]
  LC2.dc <- Root_hneg_LC2[Root_hneg_LC2$plant == "DICA" & Root_hneg_LC2$canopy == "close", 1:ncol(hneg_raw)]
  
  LC2.oed <- log2(colMeans(LC2.do)/colMeans(LC2.eo))
  LC2.ced <- log2(colMeans(LC2.dc)/colMeans(LC2.ec))
  
  volcano_hneg_oed <- data.frame(pvalue = p_oed,
                                 Log2FC = LC2.oed) 
  volcano_hneg_ced <- data.frame(pvalue = p_ced,
                                 Log2FC = LC2.ced) 
  
  Plot_Root_volcano_hneg_oed <- Plot_volcano(volcano_hneg_oed)
  Plot_Root_volcano_hneg_ced <- Plot_volcano(volcano_hneg_ced)
  ggsave("Plot_Root_volcano_hneg_oed.pdf", Plot_Root_volcano_hneg_oed + ggtitle("hneg open ERLE vs DICA") + theme(legend.position = "none"), width = 6, height = 4)
  ggsave("Plot_Root_volcano_hneg_ced.pdf", Plot_Root_volcano_hneg_ced + ggtitle("hneg close ERLE vs DICA") + theme(legend.position = "none"), width = 6, height = 4)
  
  
  
  # make table -------
  
  Diff_Root_open_plant <- rbind(volcano_rppos_oed, volcano_hpos_oed, volcano_rpneg_oed, volcano_hneg_oed)
  # Remove duplicates in "Name" column, prioritizing rows with non-NA values in "delabel" column

  
  Diff_Root_close_plant <- rbind(volcano_rppos_ced, volcano_hpos_ced, volcano_rpneg_ced, volcano_hneg_ced)

  
  write.table(Diff_Root_open_plant, file = "Diff_Root_open_plant.csv", sep = ",", row.names = FALSE)
  write.table(Diff_Root_close_plant, file = "Diff_Root_close_plant.csv", sep = ",", row.names = FALSE)
  
  
  
  