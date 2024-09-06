library(dplyr)
library(ggplot2)

setwd("~/INVA/INVA")

kegg_link_rhizo <- read.csv("Rhizo_kegg_link.csv")
kegg_link_rhizo <- kegg_link_rhizo[,c(1,6)]
kegg_link_root <- read.csv("root_kegg.csv")
kegg_link_root <- kegg_link_root[,c(1,6)]

#functions ----

#plot volcano between DICA and ERLE
Plot_volcano <- function(dataframe) {
  
  color_palette <- get_palette(palette = 'Dark2', 6)
  maingroups<-c("Amino acids, peptides, and analogues","Carboxylic acids and derivatives","Fatty acids and conjugates","Lineolic acids and derivatives","Medium-chain hydroxy acids and derivatives","Purines and purine derivatives")
  colors <- setNames(color_palette,maingroups)
  
  acid <- unique(rhizo_plant_all[which(rhizo_plant_all$Class == "Carboxylic acids and derivatives"),]$Subclass)
  acid <- acid[-c(1,6)]
  
  
  # Create a copy of the input dataframe
  volcano_data <- dataframe
  
  volcano_data$group <- volcano_data$Subclass
  
  volcano_data <- volcano_data %>%
    mutate(group = ifelse(group %in% acid,
                          "Carboxylic acids and derivatives",
                          group))
  
  volcano_data <- volcano_data %>%
    mutate(group = case_when(
      group %in% maingroups ~ group,
      TRUE ~ "other"))
  
  volcano_data[which(volcano_data$diffexpressed == "NO"),]$group <- "NO"
      
  filename <- deparse(substitute(dataframe))
  
  # Create volcano plot
  plot_volcano <- ggplot(data = volcano_data, aes(x = Log2FC, y = -log10(pvalue), col = group)) + 
    geom_vline(xintercept=c(-1.5, 1.5), col="dimgrey") +
    geom_hline(yintercept=-log10(0.05), col="dimgrey") +
    geom_point(size = 1) + 
    theme_minimal() +
    theme(legend.position = "none")+
    scale_color_manual(values = c(colors, "other" = "bisque2", "NO" = "grey"))
  
  return(plot_volcano)
}

#root ----
#root: total plant 
root_NMR_plant<-read.csv("Diff_root_plant_NMR.csv", header = T)
root_NMR_plant<-root_NMR_plant[,c(1:5,7,9:11)]
colnames <-  c("pvalue","Log2FC","diffexpressed","delabel","Name","Formula","Superclass","Class","Subclass")
colnames(root_NMR_plant) <- colnames
root_NMR_plant$Tags <- "Confirmed ID (HIgh Confidence)"
root_NMR_plant$compoundID <- paste0("NMR_",root_NMR_plant$Name)

root_LCMS_plant <- read.csv("Diff_root_plant.csv", header = T)
root_LCMS_plant$delabel[which(root_LCMS_plant$Tags != "Confirmed ID (HIgh Confidence)")] <- NA
root_plant_all <- rbind(root_LCMS_plant[,intersect(names(root_LCMS_plant), names(root_NMR_plant))],root_NMR_plant)

#root: open DvE
root_NMR_ode<-read.csv("Diff_root_OPEN_DvE_NMR.csv", header = T)
root_NMR_ode<-root_NMR_ode[,c(1:5,7,9:11)]
colnames <-  c("pvalue","Log2FC","diffexpressed","delabel","Name","Formula","Superclass","Class","Subclass")
colnames(root_NMR_ode) <- colnames
root_NMR_ode$Tags <- "Confirmed ID (HIgh Confidence)"
root_NMR_ode$compoundID <- paste0("NMR_",root_NMR_ode$Name)

root_LCMS_ode <- read.csv("Diff_Root_open_Plant.csv", header = T)
root_LCMS_ode$delabel[which(root_LCMS_ode$Tags != "Confirmed ID (HIgh Confidence)")] <- NA
root_ode_all <- rbind(root_LCMS_ode[,intersect(names(root_LCMS_ode), names(root_NMR_ode))],root_NMR_ode)

#root: close DvE cde
root_NMR_cde<-read.csv("Diff_root_CLOSE_DvE_NMR.csv", header = T)
root_NMR_cde<-root_NMR_cde[,c(1:5,7,9:11)]
colnames <-  c("pvalue","Log2FC","diffexpressed","delabel","Name","Formula","Superclass","Class","Subclass")
colnames(root_NMR_cde) <- colnames
root_NMR_cde$Tags <- "Confirmed ID (HIgh Confidence)"
root_NMR_cde$compoundID <- paste0("NMR_",root_NMR_cde$Name)

root_LCMS_cde <- read.csv("Diff_Root_close_Plant.csv", header = T)
root_LCMS_cde$delabel[which(root_LCMS_cde$Tags != "Confirmed ID (HIgh Confidence)")] <- NA
root_cde_all <- rbind(root_LCMS_cde[,intersect(names(root_LCMS_cde), names(root_NMR_cde))],root_NMR_cde)
#root: DICA OvC
root_NMR_doc<-read.csv("Diff_root_DICA_CvO_NMR.csv", header = T)
root_NMR_doc<-root_NMR_doc[,c(1:5,7,9:11)]
colnames <-  c("pvalue","Log2FC","diffexpressed","delabel","Name","Formula","Superclass","Class","Subclass")
colnames(root_NMR_doc) <- colnames
root_NMR_doc$Tags <- "Confirmed ID (HIgh Confidence)"

root_NMR_doc$diffexpressed <- "NO"
# Set diffexpressed to "OPEN" for rows where pvalue < 0.05 and Log2FC < -1.5
root_NMR_doc$diffexpressed[root_NMR_doc$pvalue < 0.05 & root_NMR_doc$Log2FC < -1.5] <- "OPEN"
# Set diffexpressed to "CLOSE" for rows where pvalue < 0.05 and Log2FC > 1.5
root_NMR_doc$diffexpressed[root_NMR_doc$pvalue < 0.05 & root_NMR_doc$Log2FC > 1.5] <- "CLOSE"
root_NMR_doc$compoundID <- paste0("NMR_",root_NMR_doc$Name)

root_LCMS_doc <- read.csv("Diff_Root_DICA_canopy.csv", header = T)
root_LCMS_doc$delabel[which(root_LCMS_doc$Tags != "Confirmed ID (HIgh Confidence)")] <- NA

root_doc_all<-  rbind(root_LCMS_doc[,intersect(names(root_LCMS_doc), names(root_NMR_doc))],root_NMR_doc)

#FT_291.12151_0.946_HPOS
#FT_293.13715_0.965_HPOS
#FT_287.12668_0.965_HPOS
#FT_243.10064_0.942_HPOS
#FT_275.12659_0.942_HPOS
#FT_354.18937_0.965_HPOS

#root: ERLE OvC
root_NMR_eoc<-read.csv("Diff_root_ERLE_CvO_NMR.csv", header = T)
root_NMR_eoc<-root_NMR_eoc[,c(1:5,7,9:11)]
colnames <-  c("pvalue","Log2FC","diffexpressed","delabel","Name","Formula","Superclass","Class","Subclass")
colnames(root_NMR_eoc) <- colnames
root_NMR_eoc$Tags <- "Confirmed ID (HIgh Confidence)"

root_NMR_eoc$diffexpressed <- "NO"
# Set diffexpressed to "OPEN" for rows where pvalue < 0.05 and Log2FC < -1.5
root_NMR_eoc$diffexpressed[root_NMR_eoc$pvalue < 0.05 & root_NMR_eoc$Log2FC < -1.5] <- "OPEN"
# Set diffexpressed to "CLOSE" for rows where pvalue < 0.05 and Log2FC > 1.5
root_NMR_eoc$diffexpressed[root_NMR_eoc$pvalue < 0.05 & root_NMR_eoc$Log2FC > 1.5] <- "CLOSE"
root_NMR_eoc$compoundID <- paste0("NMR_",root_NMR_eoc$Name)

root_LCMS_eoc <- read.csv("Diff_Root_ERLE_canopy.csv", header = T)
root_LCMS_eoc$delabel[which(root_LCMS_eoc$Tags != "Confirmed ID (HIgh Confidence)")] <- NA

root_eoc_all<-  rbind(root_LCMS_eoc[,intersect(names(root_LCMS_eoc), names(root_NMR_eoc))],root_NMR_eoc)


#open 116, 052, 073,
#FT_236.14788_6.82_HPOS,FT_193.03492_4.305_HNEG, FT_401.181_4.454_HNEG

#close 072, 064, 071, 072
#FT_308.16026_4.783_RPPOS
#FT_210.1114_4.007_HPOS
#FT_306.14536_5.058_RPNEG
#FT_518.16768_5.153_RPNEG

#rhizo------
#rhizo: total plant
rhizo_NMR_plant<-read.csv("Diff_rhizo_plant_NMR.csv", header = T)
rhizo_NMR_plant<-rhizo_NMR_plant[,c(1:7,9:11)]
colnames <-  c("pvalue","Log2FC","diffexpressed","delabel","Name", "KEGG", "Formula","Superclass","Class","Subclass")
colnames(rhizo_NMR_plant) <- colnames
rhizo_NMR_plant$Tags <- "Confirmed ID (HIgh Confidence)"
rhizo_NMR_plant$compoundID <- paste0("NMR_",rhizo_NMR_plant$Name)

rhizo_LCMS_plant <- read.csv("Diff_Rhizo_plant.csv", header = T)
rhizo_LCMS_plant <- left_join(rhizo_LCMS_plant, kegg_link_rhizo)
rhizo_LCMS_plant <- rhizo_LCMS_plant[,c(1:5, 8, 10:15)]
rhizo_LCMS_plant$delabel[which(rhizo_LCMS_plant$Tags != "Confirmed ID (HIgh Confidence)")] <- NA

rhizo_plant_all<- rhizo_LCMS_plant[,c(1,2,3,4,6,12,8,9,10,11,7,5)]
rhizo_plant_all <- rbind(rhizo_plant_all,rhizo_NMR_plant)


#rhizo: open DvE
rhizo_NMR_ode<-read.csv("Diff_rhizo_OPEN_DvE_NMR.csv", header = T)
rhizo_NMR_ode<-rhizo_NMR_ode[,c(1:7,9:11)]
colnames <-  c("pvalue","Log2FC","diffexpressed","delabel","Name", "KEGG", "Formula","Superclass","Class","Subclass")
colnames(rhizo_NMR_ode) <- colnames
rhizo_NMR_ode$Tags <- "Confirmed ID (HIgh Confidence)"
rhizo_NMR_ode$compoundID <- paste0("NMR_",rhizo_NMR_ode$Name)

rhizo_LCMS_ode <- read.csv("Diff_rhizo_open_Plant.csv", header = T)
rhizo_LCMS_ode <- left_join(rhizo_LCMS_ode, kegg_link_rhizo)
rhizo_LCMS_ode <- rhizo_LCMS_ode[,c(1:5, 8, 10:15)]
rhizo_LCMS_ode$delabel[which(rhizo_LCMS_ode$Tags != "Confirmed ID (HIgh Confidence)")] <- NA

rhizo_ode_all<- rhizo_LCMS_ode[,c(1,2,3,4,6,12,8,9,10,11,7,5)]
rhizo_ode_all <- rbind(rhizo_ode_all,rhizo_NMR_ode)
#rhizo: close DvE cde
rhizo_NMR_cde<-read.csv("Diff_rhizo_CLOSE_DvE_NMR.csv", header = T)
rhizo_NMR_cde<-rhizo_NMR_cde[,c(1:7,9:11)]
colnames <-  c("pvalue","Log2FC","diffexpressed","delabel","Name", "KEGG", "Formula","Superclass","Class","Subclass")
colnames(rhizo_NMR_cde) <- colnames
rhizo_NMR_cde$Tags <- "Confirmed ID (HIgh Confidence)"
rhizo_NMR_cde$compoundID <- paste0("NMR_",rhizo_NMR_cde$Name)

rhizo_LCMS_cde <- read.csv("Diff_rhizo_close_Plant.csv", header = T)
rhizo_LCMS_cde <- left_join(rhizo_LCMS_cde, kegg_link_rhizo)
rhizo_LCMS_cde <- rhizo_LCMS_cde[,c(1:5, 8, 10:15)]
rhizo_LCMS_cde$delabel[which(rhizo_LCMS_cde$Tags != "Confirmed ID (HIgh Confidence)")] <- NA

rhizo_cde_all<- rhizo_LCMS_cde[,c(1,2,3,4,6,12,8,9,10,11,7,5)]
rhizo_cde_all <- rbind(rhizo_cde_all,rhizo_NMR_cde)

#rhizo: DICA OvC
rhizo_NMR_doc<-read.csv("Diff_rhizo_DICA_CvO_NMR.csv", header = T)
rhizo_NMR_doc<-rhizo_NMR_doc[,c(1:7,9:11)]
colnames <-  c("pvalue","Log2FC","diffexpressed","delabel","Name", "KEGG", "Formula","Superclass","Class","Subclass")
colnames(rhizo_NMR_doc) <- colnames
rhizo_NMR_doc$Tags <- "Confirmed ID (HIgh Confidence)"
rhizo_NMR_doc$compoundID <- paste0("NMR_",rhizo_NMR_doc$Name)

rhizo_LCMS_doc <- read.csv("Diff_rhizo_DICA_canopy.csv", header = T)
rhizo_LCMS_doc <- left_join(rhizo_LCMS_doc, kegg_link_rhizo)
rhizo_LCMS_doc <- rhizo_LCMS_doc[,c(1:5, 8, 10:15)]
rhizo_LCMS_doc$delabel[which(rhizo_LCMS_doc$Tags != "Confirmed ID (HIgh Confidence)")] <- NA

rhizo_doc_all<- rhizo_LCMS_doc[,c(1,2,3,4,6,12,8,9,10,11,7,5)]
rhizo_doc_all <- rbind(rhizo_doc_all,rhizo_NMR_doc)
#rhizo: ERLE OvC
rhizo_NMR_eoc<-read.csv("Diff_rhizo_ERLE_CvO_NMR.csv", header = T)
rhizo_NMR_eoc<-rhizo_NMR_eoc[,c(1:7,9:11)]
colnames <-  c("pvalue","Log2FC","diffexpressed","delabel","Name", "KEGG", "Formula","Superclass","Class","Subclass")
colnames(rhizo_NMR_eoc) <- colnames
rhizo_NMR_eoc$Tags <- "Confirmed ID (HIgh Confidence)"
rhizo_NMR_eoc$compoundID <- paste0("NMR_",rhizo_NMR_eoc$Name)

rhizo_LCMS_eoc <- read.csv("Diff_rhizo_ERLE_canopy.csv", header = T)
rhizo_LCMS_eoc <- left_join(rhizo_LCMS_eoc, kegg_link_rhizo)
rhizo_LCMS_eoc <- rhizo_LCMS_eoc[,c(1:5, 8, 10:15)]
rhizo_LCMS_eoc$delabel[which(rhizo_LCMS_eoc$Tags != "Confirmed ID (HIgh Confidence)")] <- NA

rhizo_eoc_all<- rhizo_LCMS_eoc[,c(1,2,3,4,6,12,8,9,10,11,7,5)]
rhizo_eoc_all <- rbind(rhizo_eoc_all,rhizo_NMR_eoc)




#stack plot-----
#root close ----
defense_chemicals <- c("Cinnamic acids and derivatives", "Flavonoids",
                       "Prenol lipids",  "Benzene and substituted derivatives", "Azoles", 	
                       "Fatty Acyls", "Indoles and derivatives","Benzopyran", "Phenols" )

counts <- root_cde_all %>%
  filter(diffexpressed != "NO" & !is.na(Class)) %>%
  group_by(Class, diffexpressed) %>%
  summarise(count = n())

# Calculate total counts
sum_counts <- aggregate(count ~ Class, data = counts, FUN = sum)

# Calculate total counts for each Class
total <- as.data.frame(table(root_cde_all$Class)) %>%
  filter(Var1 != "")
colnames(total) <- c("Class", "count")
total$diffexpressed <- "total"

total <- total[total$Class %in% defense_chemicals,]
sum_counts <- sum_counts[sum_counts$Class %in% defense_chemicals,]

# Merge total counts with sum_counts
total <- merge(total, as.data.frame(sum_counts), by = "Class", suffixes = c(".total", ".count"), all.x = TRUE)
total <- replace_na(total, list(count.count = 0))
total$count <- total$count.total - total$count.count

# Determine order of superclasses
order <- total[order(total$count.total), ]$Class

# Combine counts and total counts
combined <- rbind(counts[counts$Class %in% defense_chemicals,], total[, c(1, 3, 5)])

# Generate all combinations of Super.class and diffexpressed
all_combinations <- expand.grid(Class = unique(total$Class),
                                diffexpressed = c("ERLE", "DICA", "total"),
                                stringsAsFactors = FALSE)

# Merge all_combinations with combined
combined_df <- merge(all_combinations, combined, by = c("Class", "diffexpressed"), all.x = TRUE)
combined_df$count[is.na(combined_df$count)] <- 0

combined_df$Class <- factor(combined_df$Class, levels = order) 
combined_df$diffexpressed <- factor(combined_df$diffexpressed, levels = c("total", "DICA", "ERLE", "close", "CLOSE")) 

# Define color palette
color_palette <- c("ERLE" = "#50124C", "DICA" = "#4DAF4A", "close" = "#FBB413", "CLOSE" = "#06B4BC", "total" = "grey")

# Create the plot
Root_stack_close <- ggplot(combined_df, aes(fill = diffexpressed, x = count, y = Class)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = color_palette) +
  theme_minimal()+
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank())+
  theme(axis.title = element_blank(), legend.position = "none")

#root open ----
defense_chemicals <- c("Cinnamic acids and derivatives", "Flavonoids",
                       "Prenol lipids",  "Benzene and substituted derivatives", "Azoles", 	
                       "Fatty Acyls", "Indoles and derivatives","Benzopyran", "Phenols" )

counts <- root_ode_all %>%
  filter(diffexpressed != "NO" & !is.na(Class)) %>%
  group_by(Class, diffexpressed) %>%
  summarise(count = n())

# Calculate total counts
sum_counts <- aggregate(count ~ Class, data = counts, FUN = sum)

# Calculate total counts for each Class
total <- as.data.frame(table(root_cde_all$Class)) %>%
  filter(Var1 != "")
colnames(total) <- c("Class", "count")
total$diffexpressed <- "total"

total <- total[total$Class %in% defense_chemicals,]
sum_counts <- sum_counts[sum_counts$Class %in% defense_chemicals,]

# Merge total counts with sum_counts
total <- merge(total, as.data.frame(sum_counts), by = "Class", suffixes = c(".total", ".count"), all.x = TRUE)
total <- replace_na(total, list(count.count = 0))
total$count <- total$count.total - total$count.count

# Determine order of superclasses
order <- total[order(total$count.total), ]$Class

# Combine counts and total counts
combined <- rbind(counts[counts$Class %in% defense_chemicals,], total[, c(1, 3, 5)])

# Generate all combinations of Super.class and diffexpressed
all_combinations <- expand.grid(Class = unique(total$Class),
                                diffexpressed = c("ERLE", "DICA", "total"),
                                stringsAsFactors = FALSE)

# Merge all_combinations with combined
combined_df <- merge(all_combinations, combined, by = c("Class", "diffexpressed"), all.x = TRUE)
combined_df$count[is.na(combined_df$count)] <- 0

combined_df$Class <- factor(combined_df$Class, levels = order) 
combined_df$diffexpressed <- factor(combined_df$diffexpressed, levels = c("total", "DICA", "ERLE", "close", "CLOSE")) 

# Define color palette
color_palette <- c("ERLE" = "#50124C", "DICA" = "#4DAF4A", "close" = "#FBB413", "CLOSE" = "#06B4BC", "total" = "grey")

# Create the plot
Root_stack_open <- ggplot(combined_df, aes(fill = diffexpressed, x = count, y = Class)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = color_palette) +
  theme_minimal()+
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank())+
  theme(axis.title = element_blank(), legend.position = "none") 


ggsave("Root_open_stack.pdf", Root_stack_open, height = 4, width = 6)
ggsave("Root_close_stack.pdf", Root_stack_close , height = 4, width = 6)

#rhizo close ----
defense_chemicals <- c("Cinnamic acids and derivatives", "Flavonoids",
                       "Prenol lipids",  "Benzene and substituted derivatives", "Azoles", 	
                       "Fatty Acyls", "Indoles and derivatives","Benzopyran", "Phenols" )

counts <- rhizo_cde_all %>%
  filter(diffexpressed != "NO" & !is.na(Class)) %>%
  group_by(Class, diffexpressed) %>%
  summarise(count = n())

# Calculate total counts
sum_counts <- aggregate(count ~ Class, data = counts, FUN = sum)

# Calculate total counts for each Class
total <- as.data.frame(table(rhizo_cde_all$Class)) %>%
  filter(Var1 != "")
colnames(total) <- c("Class", "count")
total$diffexpressed <- "total"

total <- total[total$Class %in% defense_chemicals,]
sum_counts <- sum_counts[sum_counts$Class %in% defense_chemicals,]

# Merge total counts with sum_counts
total <- merge(total, as.data.frame(sum_counts), by = "Class", suffixes = c(".total", ".count"), all.x = TRUE)
total <- replace_na(total, list(count.count = 0))
total$count <- total$count.total - total$count.count

# Determine order of superclasses
order <- total[order(total$count.total), ]$Class

# Combine counts and total counts
combined <- rbind(counts[counts$Class %in% defense_chemicals,], total[, c(1, 3, 5)])

# Generate all combinations of Super.class and diffexpressed
all_combinations <- expand.grid(Class = unique(total$Class),
                                diffexpressed = c("ERLE", "DICA", "total"),
                                stringsAsFactors = FALSE)

# Merge all_combinations with combined
combined_df <- merge(all_combinations, combined, by = c("Class", "diffexpressed"), all.x = TRUE)
combined_df$count[is.na(combined_df$count)] <- 0

combined_df$Class <- factor(combined_df$Class, levels = order) 
combined_df$diffexpressed <- factor(combined_df$diffexpressed, levels = c("total", "DICA", "ERLE", "close", "CLOSE")) 

# Define color palette
color_palette <- c("ERLE" = "#50124C", "DICA" = "#4DAF4A", "close" = "#FBB413", "CLOSE" = "#06B4BC", "total" = "grey")

# Create the plot
rhizo_stack_close <- ggplot(combined_df, aes(fill = diffexpressed, x = count, y = Class)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = color_palette) +
  theme_minimal()+
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank())+
  theme(axis.title = element_blank(), legend.position = "none") 

#rhizo open ----
defense_chemicals <- c("Cinnamic acids and derivatives", "Flavonoids",
                       "Prenol lipids",  "Benzene and substituted derivatives", "Azoles", 	
                       "Fatty Acyls", "Indoles and derivatives","Benzopyran", "Phenols" )

counts <- rhizo_ode_all %>%
  filter(diffexpressed != "NO" & !is.na(Class)) %>%
  group_by(Class, diffexpressed) %>%
  summarise(count = n())

# Calculate total counts
sum_counts <- aggregate(count ~ Class, data = counts, FUN = sum)

# Calculate total counts for each Class
total <- as.data.frame(table(rhizo_cde_all$Class)) %>%
  filter(Var1 != "")
colnames(total) <- c("Class", "count")
total$diffexpressed <- "total"

total <- total[total$Class %in% defense_chemicals,]
sum_counts <- sum_counts[sum_counts$Class %in% defense_chemicals,]

# Merge total counts with sum_counts
total <- merge(total, as.data.frame(sum_counts), by = "Class", suffixes = c(".total", ".count"), all.x = TRUE)
total <- replace_na(total, list(count.count = 0))
total$count <- total$count.total - total$count.count

# Determine order of superclasses
order <- total[order(total$count.total), ]$Class

# Combine counts and total counts
combined <- rbind(counts[counts$Class %in% defense_chemicals,], total[, c(1, 3, 5)])

# Generate all combinations of Super.class and diffexpressed
all_combinations <- expand.grid(Class = unique(total$Class),
                                diffexpressed = c("ERLE", "DICA", "total"),
                                stringsAsFactors = FALSE)

# Merge all_combinations with combined
combined_df <- merge(all_combinations, combined, by = c("Class", "diffexpressed"), all.x = TRUE)
combined_df$count[is.na(combined_df$count)] <- 0

combined_df$Class <- factor(combined_df$Class, levels = order) 
combined_df$diffexpressed <- factor(combined_df$diffexpressed, levels = c("total", "DICA", "ERLE", "close", "CLOSE")) 

# Define color palette
color_palette <- c("ERLE" = "#50124C", "DICA" = "#4DAF4A", "close" = "#FBB413", "CLOSE" = "#06B4BC", "total" = "grey")

# Create the plot
rhizo_stack_open <- ggplot(combined_df, aes(fill = diffexpressed, x = count, y = Class)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = color_palette) +
  theme_minimal()+
  theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank())+
  theme(axis.title = element_blank(), legend.position = "none")



ggsave("rhizo_open_stack.pdf", rhizo_stack_open, height = 3.5, width = 6)
ggsave("rhizo_close_stack.pdf", rhizo_stack_close , height = 3.5, width = 6)








#volcano plots ----------

ggsave("root_plant.pdf", Plot_volcano(root_plant_all), width =5, height = 4)
ggsave("root_eoc.pdf", Plot_volcano(root_eoc_all), width =5, height = 4)
ggsave("root_doc.pdf", Plot_volcano(root_doc_all), width =5, height = 4)

ggsave("rhizo_plant.pdf", Plot_volcano(rhizo_plant_all), width =5, height = 4)
ggsave("rhizo_eoc.pdf", Plot_volcano(rhizo_eoc_all), width =5, height = 4)
ggsave("rhizo_doc.pdf", Plot_volcano(rhizo_doc_all), width =5, height = 4)

# List of center parts of your dataframe names
center_parts <- c("plant", "cde", "ode", "doc", "eoc")


save_df <- function(df, prefix, center_part) {
  filename <- paste0(prefix, "_", center_part, "_all.csv") 
  write.csv(df, filename, row.names = FALSE)  
}

# Iterate over names and save
for (part in center_parts) {
  # Get dataframes assuming they exist in your workspace
  root_df <- get(paste0("root_", part, "_all")) 
  rhizo_df <- get(paste0("rhizo_", part, "_all")) 
  
  # Save both "root_" and "rhizo_" variations
  save_df(root_df, prefix = "root", center_part = part) 
  save_df(rhizo_df, prefix = "rhizo", center_part = part) 
}
