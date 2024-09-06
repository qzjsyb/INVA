library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(FSA)

#pathway viwer and enrichment analysis.
#read all root data.
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

#Load data ------------------------------------
setwd('~/INVA/INVA')
metadata <- read.table("metadata.txt", header = TRUE, sep = "\t", row.names = 1)
metadata$block <- as.factor(metadata$block)
Root_compounds <- read.csv("Root_compound_link.csv", header = TRUE, row.names = 1)
Root_compounds <- rownames_to_column(Root_compounds, var = "compoundID")
LCMS_kegg <- read.csv("root_kegg.csv", header = TRUE, row.names = 1)
LCMS_kegg$compoundID <- rownames(LCMS_kegg)

rppos_raw <- read_csv("Root_LCMS_RPPOS.csv") %>%
  #Make the compound IDs the row names
  column_to_rownames('Compound ID') %>%
  #Transpose the matrix to make it so samples are rows and compounds are columns
  t()  %>%
  as.data.frame()
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

hpos_raw <- read_csv("Root_LCMS_HPOS.csv") %>%
  #Make the compound IDs the row names
  column_to_rownames('Compound ID') %>%
  #Transnege the matrix to make it so samples are rows and compounds are columns
  t()  %>%
  as.data.frame()

### here is for testing different transforming methods.

hpos_pareto <- hpos_raw %>%
  log_transform() %>%
  pareto_scale()

#Combine metadata with metabolite data
Root_hpos_final  <- hpos_pareto %>%
  rownames_to_column('SampleID_root') %>%
  mutate(SampleID_root = sub("_HILIC_Pos$", "", SampleID_root)) %>% 
  right_join(metadata)%>%
  column_to_rownames("SampleID_root")

hneg_raw <- read_csv("Root_LCMS_HNEG.csv") %>%
  #Make the compound IDs the row names
  column_to_rownames('Compound ID') %>%
  #Transnege the matrix to make it so samples are rows and compounds are columns
  t()  %>%
  as.data.frame()


hneg_pareto <- hneg_raw %>%
  log_transform() %>%
  pareto_scale()

#Combine metadata with metabolite data
Root_hneg_final  <- hneg_pareto %>%
  rownames_to_column('SampleID_root') %>%
  mutate(SampleID_root = sub("_HILIC_Neg$", "", SampleID_root)) %>% 
  left_join(metadata)%>%
  column_to_rownames("SampleID_root")

NMR_root <- as.data.frame(t(read.table("NMR_root_clean.txt", header = TRUE, sep = "\t", row.names =1)))
NMR_chem <- read.csv("NMR_chem.csv", header = TRUE)
NMR_chem <- NMR_chem %>% 
  arrange(Class)
NMR_chem$metabolite <- gsub("^([0-9]+)", "X\\1", NMR_chem$metabolite)  # Add X in the beginning
NMR_chem$metabolite <- gsub("-", ".", NMR_chem$metabolite)              # Change - to .
NMR_chem$metabolite <- gsub(" ", ".", NMR_chem$metabolite)              # Replace spaces with .


NMR_pareto <- NMR_root %>%
  log_transform() %>%
  pareto_scale()
rownames(NMR_pareto) <- gsub("\\.","-",rownames(NMR_pareto))

NMR_final <- cbind(NMR_pareto[rownames(metadata),], metadata)
Root_rppos_final <-Root_rppos_final[metadata[-which(metadata$SampleID_root == "R27"),]$SampleID_root,]
Root_rpneg_final <-Root_rpneg_final[metadata$SampleID_root,]
Root_hpos_final <-Root_hpos_final[metadata$SampleID_root,]
Root_hneg_final <-Root_hneg_final[metadata$SampleID_root,]


#ggsave("adenine.png", plot_point("C00049"), height = 1.5, width = 1.5, dpi = 300)
#ggsave("guanine.png", plot_point("C00242"), height = 1.5, width = 1.5, dpi = 300)
#ggsave("Hypoxanthine.png", plot_point("C00262"), height = 1.5, width = 1.5, dpi = 300)
#ggsave("allantoin.png", plot_point("C02350"), height = 1.5, width = 1.5, dpi = 300)

#ggsave("arginine.png", plot_point("C00062"), height = 1.5, width = 1.5, dpi = 300)
#ggsave("citrulline.png", plot_point("C00327"), height = 1.5, width = 1.5, dpi = 300)
#ggsave("glutamine.png", plot_point("C00064"), height = 1.5, width = 1.5, dpi = 300)
#ggsave("glutamate.png", plot_point("C00025"), height = 1.5, width = 1.5, dpi = 300)

adenine_ids <- LCMS_kegg[LCMS_kegg$KEGG2 == "C00049",]$compoundID
guanine_ids <- LCMS_kegg[LCMS_kegg$KEGG1 == "C00042",]$compoundID
Hypoxanthine_ids <- LCMS_kegg[LCMS_kegg$KEGG1 == "C00262",]$compoundID
allantoin_ids <- LCMS_kegg[LCMS_kegg$KEGG1 == "C02350",]$compoundID
arginine_ids <- LCMS_kegg[LCMS_kegg$KEGG2 == "C00062",]$compoundID
citrulline_ids <- LCMS_kegg[LCMS_kegg$KEGG1 == "C00327",]$compoundID
glutamine_ids <- LCMS_kegg[LCMS_kegg$KEGG2 == "C00064",]$compoundID
glutamate_ids <- LCMS_kegg[LCMS_kegg$KEGG2 == "C00025",]$compoundID

extract_compound_data <- function(compound_ids, dataframes) {
  compound_data <- list()
  max_rows <- max(sapply(dataframes, nrow))
  
  for (df_name in names(dataframes)) {
    df <- dataframes[[df_name]]
    matching_columns <- df[, colnames(df) %in% compound_ids, drop = FALSE]
    
    if (ncol(matching_columns) > 0) {
      # Pad with NAs if the number of rows is less than the maximum number of rows
      if (nrow(matching_columns) < max_rows) {
        padding <- matrix(NA, nrow = max_rows - nrow(matching_columns), ncol = ncol(matching_columns))
        colnames(padding) <- colnames(matching_columns)
        matching_columns <- rbind(matching_columns, padding)
      }
      compound_data[[df_name]] <- matching_columns
    }
  }
  
  if (length(compound_data) > 0) {
    combined_data <- do.call(cbind, compound_data)
  } else {
    combined_data <- data.frame()
  }
  return(combined_data)
}

# Dataframes to search
dataframes_to_search <- list(
  Root_rpneg_final = Root_rpneg_final,
  Root_hpos_final = Root_hpos_final,
  Root_hneg_final = Root_hneg_final,
  Root_rppos_final = Root_rppos_final
)

adenine <- extract_compound_data(adenine_ids, dataframes_to_search)
guanine <- extract_compound_data(guanine_ids, dataframes_to_search)
Hypoxanthine <- extract_compound_data(Hypoxanthine_ids, dataframes_to_search)
allantoin <- extract_compound_data(allantoin_ids, dataframes_to_search)
arginine <- extract_compound_data(arginine_ids, dataframes_to_search)
citrulline <- extract_compound_data(citrulline_ids, dataframes_to_search)
glutamine <- extract_compound_data(glutamine_ids, dataframes_to_search)
glutamate <- extract_compound_data(glutamate_ids, dataframes_to_search)

# List of compound names to search for
compounds <- c("Adenine", "Guanine", "Hypoxanthine", "Allantoin", "Arginine", "Citrulline", "Glutamine", "Glutamate")

# Check if the compounds exist in the column names of NMR_final
matching_columns <- colnames(NMR_final)[colnames(NMR_final) %in% compounds]

# Extract the matching columns
extracted_columns <- NMR_final[, matching_columns, drop = FALSE]
glutamine$NMR <-  extracted_columns$Glutamine
glutamate$NMR <- extracted_columns$Glutamate

#all data are there, now run kruskal test

df <- adenine %>%
  rownames_to_column("SampleID_root")%>%
  left_join(., metadata[,c(2:4)])
df <- melt(df, id.vars = c("SampleID_root", "plant", "canopy"))
df$group <- paste0(df$plant,df$canopy)

dunnTest(value ~ group, data = df) 

df <- guanine %>%
  rownames_to_column("SampleID_root")%>%
  left_join(., metadata[,c(2:4)])
df <- melt(df, id.vars = c("SampleID_root", "plant", "canopy"))
df$group <- paste0(df$plant,df$canopy)

dunnTest(value ~ group, data = df) 



metacol<- c("block","plant","canopy","SampleID_root","SampleID_bulk","SampleID_rhizo","paring")

calculate_average_intensity <- function(data) {
  # Define unique combinations of plant and canopy
  treatments <- unique(data[, c("plant", "canopy")])
  # Initialize a list to store results
  results <- list()
  # Loop through each treatment combination
  for (i in 1:nrow(treatments)) {
    plant <- treatments$plant[i]
    canopy <- treatments$canopy[i]
    # Subset data for the current treatment combination
    subset_data <- data[data$plant == plant & data$canopy == canopy, ]
    cols_to_remove <- which(names(subset_data) %in% metacol)
    # Calculate average intensity for each metabolite
    average_intensity <- colMeans(subset_data[, -cols_to_remove], na.rm = TRUE)
    # Store results in a named vector
    results[[paste(plant, canopy, sep = "_")]] <- average_intensity
  }
  # Return the list of results
  return(results)
}

NMR.ave <- calculate_average_intensity(NMR_final)
rppos.ave <- calculate_average_intensity(Root_rppos_final)
rpneg.ave <- calculate_average_intensity(Root_rpneg_final)
hpos.ave <- calculate_average_intensity(Root_hpos_final)
hneg.ave <- calculate_average_intensity(Root_hneg_final)

for (result_name in names(rppos.ave)) {
  # Transpose the dataframe
  df <- as.data.frame(rppos.ave[[result_name]])
  colnames(df) <- "value"  # Set column name as "value"
  df$compoundID <- rownames(df)  # Create a new column "compoundID" from rownames
  
  # Left join with LCMS_kegg by compoundID
  joined_df <- merge(df, LCMS_kegg, by = "compoundID", all.x = TRUE)
  
  #keep KEGG and value
  kegg_df <- joined_df[, c("value", "KEGG1", "KEGG2", "KEGG3", "KEGG4")]
  melt_df <- melt(kegg_df, id.vars = "value",value.name = "kegg_id")
  final_df <- melt_df[!is.na(melt_df$kegg_id) & melt_df$kegg_id != "", ]
  
   # Save the joined dataframe
  assign(paste0(result_name, "_kegg_rppos"), final_df)
}

for (result_name in names(rppos.ave)) {
  # Transpose the dataframe
  df <- as.data.frame(rppos.ave[[result_name]])
  colnames(df) <- "value"  # Set column name as "value"
  df$compoundID <- rownames(df)  # Create a new column "compoundID" from rownames
  
  # Left join with LCMS_kegg by compoundID
  joined_df <- merge(df, LCMS_kegg, by = "compoundID", all.x = TRUE)
  
  #keep KEGG and value
  kegg_df <- joined_df[, c("value", "KEGG1", "KEGG2", "KEGG3", "KEGG4")]
  melt_df <- melt(kegg_df, id.vars = "value",value.name = "kegg_id")
  final_df <- melt_df[!is.na(melt_df$kegg_id) & melt_df$kegg_id != "", ]
  final_df <- final_df[,c("value","kegg_id")]
  
   # Save the joined dataframe
  assign(paste0(result_name, "_kegg_rppos"), final_df)
}
for (result_name in names(rpneg.ave)) {
  # Transpose the dataframe
  df <- as.data.frame(rpneg.ave[[result_name]])
  colnames(df) <- "value"  # Set column name as "value"
  df$compoundID <- rownames(df)  # Create a new column "compoundID" from rownames
  
  # Left join with LCMS_kegg by compoundID
  joined_df <- merge(df, LCMS_kegg, by = "compoundID", all.x = TRUE)
  
  #keep KEGG and value
  kegg_df <- joined_df[, c("value", "KEGG1", "KEGG2", "KEGG3", "KEGG4")]
  melt_df <- melt(kegg_df, id.vars = "value",value.name = "kegg_id")
  final_df <- melt_df[!is.na(melt_df$kegg_id) & melt_df$kegg_id != "", ]
  final_df <- final_df[,c("value","kegg_id")]
  
  # Save the joined dataframe
  assign(paste0(result_name, "_kegg_rpneg"), final_df)
}
for (result_name in names(hpos.ave)) {
  # Transpose the dataframe
  df <- as.data.frame(hpos.ave[[result_name]])
  colnames(df) <- "value"  # Set column name as "value"
  df$compoundID <- rownames(df)  # Create a new column "compoundID" from rownames
  
  # Left join with LCMS_kegg by compoundID
  joined_df <- merge(df, LCMS_kegg, by = "compoundID", all.x = TRUE)
  
  #keep KEGG and value
  kegg_df <- joined_df[, c("value", "KEGG1", "KEGG2", "KEGG3", "KEGG4")]
  melt_df <- melt(kegg_df, id.vars = "value",value.name = "kegg_id")
  final_df <- melt_df[!is.na(melt_df$kegg_id) & melt_df$kegg_id != "", ]
  final_df <- final_df[,c("value","kegg_id")]
  
  # Save the joined dataframe
  assign(paste0(result_name, "_kegg_hpos"), final_df)
}
for (result_name in names(hneg.ave)) {
  # Transpose the dataframe
  df <- as.data.frame(hneg.ave[[result_name]])
  colnames(df) <- "value"  # Set column name as "value"
  df$compoundID <- rownames(df)  # Create a new column "compoundID" from rownames
  
  # Left join with LCMS_kegg by compoundID
  joined_df <- merge(df, LCMS_kegg, by = "compoundID", all.x = TRUE)
  
  #keep KEGG and value
  kegg_df <- joined_df[, c("value", "KEGG1", "KEGG2", "KEGG3", "KEGG4")]
  melt_df <- melt(kegg_df, id.vars = "value",value.name = "kegg_id")
  final_df <- melt_df[!is.na(melt_df$kegg_id) & melt_df$kegg_id != "", ]
  final_df <- final_df[,c("value","kegg_id")]
  
  # Save the joined dataframe
  assign(paste0(result_name, "_kegg_hneg"), final_df)
}
for (result_name in names(NMR.ave)) {
  # Transpose the dataframe
  df <- as.data.frame(NMR.ave[[result_name]])
  colnames(df) <- "value"  # Set column name as "value"
  df$metabolite <- rownames(df)  # Create a new column "compoundID" from rownames
  
  
  
  # Left join with LCMS_kegg by compoundID
  joined_df <- merge(df, NMR_chem, by = "metabolite", all.x = TRUE)
  
  #keep KEGG and value
  kegg_df <- joined_df[, c("value", "KEGG.Compound.ID")]
  colnames(kegg_df) <- c("value", "kegg_id")
  final_df <- kegg_df[!is.na(kegg_df$kegg_id) & kegg_df$kegg_id != "", ]
  
  # Save the joined dataframe
  assign(paste0(result_name, "_kegg_NMR"), final_df)
}

dc <- rbind(DICA_close_kegg_hneg, DICA_close_kegg_hpos, DICA_close_kegg_rppos, DICA_close_kegg_rpneg, DICA_close_kegg_NMR)
do <- rbind(DICA_open_kegg_hneg, DICA_open_kegg_hpos, DICA_open_kegg_rppos, DICA_open_kegg_rpneg, DICA_open_kegg_NMR)
ec <- rbind(ERLE_close_kegg_hneg, ERLE_close_kegg_hpos, ERLE_close_kegg_rppos, ERLE_close_kegg_rpneg, ERLE_close_kegg_NMR)
eo <- rbind(ERLE_open_kegg_hneg, ERLE_open_kegg_hpos, ERLE_open_kegg_rppos, ERLE_open_kegg_rpneg, ERLE_open_kegg_NMR)

dc_agg <- aggregate(value ~ kegg_id, data = dc, FUN = mean)
do_agg <- aggregate(value ~ kegg_id, data = do, FUN = mean)
ec_agg <- aggregate(value ~ kegg_id, data = ec, FUN = mean)
eo_agg <- aggregate(value ~ kegg_id, data = eo, FUN = mean)

combined <- merge(dc_agg, do_agg, by = "kegg_id", all = TRUE)
names(combined)[names(combined) == "value.x"] <- "dc"
names(combined)[names(combined) == "value.y"] <- "do"
combined <- merge(combined, ec_agg, by = "kegg_id", all = TRUE)
names(combined)[names(combined) == "value"] <- "ec"
combined <- merge(combined, eo_agg, by = "kegg_id", all = TRUE)
names(combined)[names(combined) == "value"] <- "eo"
combined <- combined[,c("kegg_id", "do","eo","dc", "ec" )]

write.table(combined, "root_intensity_kegg.csv", row.names = FALSE, sep = ",")

plot_point <- function(kegg_id) {
  # Subset data for the specified KEGG ID
  subset_data <- combined[combined$kegg_id == kegg_id, ]
  
  subset_data <- as.data.frame(t(subset_data))
  subset_data <- subset_data[-1, , drop = FALSE]
  colnames(subset_data) <- "value"
  subset_data$canopy <- c("open","open","close","close")
  subset_data$plant <- c("DICA","ERLE","DICA","ERLE")
  
# calculate CI 
  list_NMR <- NMR_chem[which(NMR_chem$KEGG.Compound.ID == kegg_id),]$metabolite
  selected_NMR <- NMR_final[, names(NMR_final) %in% list_NMR, drop = FALSE]
  
  list_LCMS <- c(LCMS_kegg[which(LCMS_kegg$KEGG1 == kegg_id),]$compoundID, 
            LCMS_kegg[which(LCMS_kegg$KEGG2 == kegg_id),]$compoundID,  
            LCMS_kegg[which(LCMS_kegg$KEGG3 == kegg_id),]$compoundID,
            LCMS_kegg[which(LCMS_kegg$KEGG4 == kegg_id),]$compoundID)
  selected_rppos <- Root_rppos_final[, names(Root_rppos_final) %in% list_LCMS, drop = FALSE]
  selected_rpneg <- Root_rpneg_final[, names(Root_rpneg_final) %in% list_LCMS, drop = FALSE]
  selected_hpos <- Root_hpos_final[, names(Root_hpos_final) %in% list_LCMS, drop = FALSE]
  selected_hneg <- Root_hneg_final[, names(Root_hneg_final) %in% list_LCMS, drop = FALSE]
  
dataframes <- list(selected_NMR, selected_rppos, selected_rpneg, selected_hpos, selected_hneg)
  
  # Define combinations of plant and canopy
  combinations <- expand.grid(plant = c("DICA", "ERLE"), canopy = c("open", "close"))
  
  # Initialize an empty list to store results
  result_list <- list()
  
  # Loop over each dataframe and combination of plant and canopy
  for (i in 1:nrow(combinations)) {
    plant <- combinations$plant[i]
    canopy <- combinations$canopy[i]
    combined_result <- NULL
    # Loop over each dataframe
    for (df in dataframes) {
      # Subset the dataframe based on the current combination of plant and canopy
      result <- df[which(metadata$plant == plant & metadata$canopy == canopy), ]
      # If combined_result is NULL, assign the result directly
      if (is.null(combined_result)) {
        combined_result <- result
      } else {
        # Otherwise, combine the result with the existing combined_result using rbind
        combined_result <- rbind(combined_result, result)
      }
    }
    # Assign the combined result to the key (combination of plant and canopy)
    result_list[[paste0(plant, ".", canopy)]] <- combined_result
  }
  
  ci.do <- as.vector(unlist(result_list[["DICA.open"]]))
  ci.dc <- as.vector(unlist(result_list[["DICA.close"]]))
  ci.eo <- as.vector(unlist(result_list[["ERLE.open"]]))
  ci.ec <- as.vector(unlist(result_list[["ERLE.close"]]))
  
  subset_data$ave <- c(mean(ci.do,na.rm = T),mean(ci.eo, na.rm = T),mean(ci.dc,na.rm = T),mean(ci.ec,na.rm = T))
  subset_data$std <- c(sd(ci.do,na.rm = T),sd(ci.eo, na.rm = T),sd(ci.dc,na.rm = T),sd(ci.ec,na.rm = T))
  subset_data$upper <- subset_data$ave + 1.96 * subset_data$std
  subset_data$lower <- subset_data$ave - 1.96 * subset_data$std 
  
  subset_data$canopy <- factor(subset_data$canopy, levels = c("open", "close"))
  
  # Create the ggplot
  p <- ggplot(subset_data, aes(x = canopy, y = ave, color = plant)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, alpha = 0.6,position = position_dodge(width = 0.2)) +
    geom_point(size = 3, position = position_dodge(width = 0.2)) +
    scale_color_manual(values = c("DICA" = "#4AAF49", "ERLE" = "#50124C")) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 15) +
    theme(plot.background = element_blank())+
    theme(legend.position = "none", panel.grid = element_blank())
  
  # Return the ggplot object
  return(p)
}




ggsave("adenine.png", plot_point("C00147"), height = 1.5, width = 1.5, dpi = 300)
ggsave("guanine.png", plot_point("C00242"), height = 1.5, width = 1.5, dpi = 300)
ggsave("Hypoxanthine.png", plot_point("C00262"), height = 1.5, width = 1.5, dpi = 300)
ggsave("allantoin.png", plot_point("C02350"), height = 1.5, width = 1.5, dpi = 300)

ggsave("arginine.png", plot_point("C00062"), height = 1.5, width = 1.5, dpi = 300)
ggsave("citrulline.png", plot_point("C00327"), height = 1.5, width = 1.5, dpi = 300)
ggsave("glutamine.png", plot_point("C00064"), height = 1.5, width = 1.5, dpi = 300)
ggsave("glutamate.png", plot_point("C00025"), height = 1.5, width = 1.5, dpi = 300)



ggsave("adenine.pdf", plot_point("C00147"), height = 1.5, width = 1.5)
ggsave("guanine.pdf", plot_point("C00242"), height = 1.5, width = 1.5)
ggsave("Hypoxanthine.pdf", plot_point("C00262"), height = 1.5, width = 1.5)
ggsave("allantoin.pdf", plot_point("C02350"), height = 1.5, width = 1.5)

ggsave("arginine.pdf", plot_point("C00062"), height = 1.5, width = 1.5)
ggsave("citrulline.pdf", plot_point("C00327"), height = 1.5, width = 1.5)
ggsave("glutamine.pdf", plot_point("C00064"), height = 1.5, width = 1.5)
ggsave("glutamate.pdf", plot_point("C00025"), height = 1.5, width = 1.5)


#perform Dunn's test
dunn("C00147")
dunn("C00242")
dunn("C00262")
dunn("C02350")
dunn("C00062")
dunn("C00327")
dunn("C00064")
dunn("C00025")

dunn <- function(kegg_id) {
  # Subset data for the specified KEGG ID
  subset_data <- combined[combined$kegg_id == kegg_id, ]
  
  subset_data <- as.data.frame(t(subset_data))
  subset_data <- subset_data[-1, , drop = FALSE]
  colnames(subset_data) <- "value"
  subset_data$canopy <- c("open","open","close","close")
  subset_data$plant <- c("DICA","ERLE","DICA","ERLE")

list_NMR <- NMR_chem[which(NMR_chem$KEGG.Compound.ID == kegg_id),]$metabolite
selected_NMR <- NMR_final[, names(NMR_final) %in% list_NMR, drop = FALSE]

list_LCMS <- c(LCMS_kegg[which(LCMS_kegg$KEGG1 == kegg_id),]$compoundID, 
               LCMS_kegg[which(LCMS_kegg$KEGG2 == kegg_id),]$compoundID,  
               LCMS_kegg[which(LCMS_kegg$KEGG3 == kegg_id),]$compoundID,
               LCMS_kegg[which(LCMS_kegg$KEGG4 == kegg_id),]$compoundID)
selected_rppos <- Root_rppos_final[, names(Root_rppos_final) %in% list_LCMS, drop = FALSE]
selected_rpneg <- Root_rpneg_final[, names(Root_rpneg_final) %in% list_LCMS, drop = FALSE]
selected_hpos <- Root_hpos_final[, names(Root_hpos_final) %in% list_LCMS, drop = FALSE]
selected_hneg <- Root_hneg_final[, names(Root_hneg_final) %in% list_LCMS, drop = FALSE]

dataframes <- list(selected_NMR, selected_rppos, selected_rpneg, selected_hpos, selected_hneg)

# Define combinations of plant and canopy
combinations <- expand.grid(plant = c("DICA", "ERLE"), canopy = c("open", "close"))

# Initialize an empty list to store results
result_list <- list()

# Loop over each dataframe and combination of plant and canopy
for (i in 1:nrow(combinations)) {
  plant <- combinations$plant[i]
  canopy <- combinations$canopy[i]
  combined_result <- NULL
  # Loop over each dataframe
  for (df in dataframes) {
    # Subset the dataframe based on the current combination of plant and canopy
    result <- df[which(metadata$plant == plant & metadata$canopy == canopy), ]
    # If combined_result is NULL, assign the result directly
    if (is.null(combined_result)) {
      combined_result <- result
    } else {
      # Otherwise, combine the result with the existing combined_result using rbind
      combined_result <- rbind(combined_result, result)
    }
  }
  # Assign the combined result to the key (combination of plant and canopy)
  result_list[[paste0(plant, ".", canopy)]] <- combined_result
}


result <- data.frame(
  value = unlist(result_list),  
  group = rep(names(result_list), times = sapply(result_list, length)) 
)

return(dunn_test(value~group, data = result))}

