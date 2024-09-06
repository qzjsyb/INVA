chemRICHFix <- function(x){
  #Data Check:
  if(sum(colnames(x) %in% c('compound_name', 'effect_size', 'pvalue', 'set')) != 4){
    stop('Data Incorrect: Make sure column names are compound_name, effect_size, pvalue, and set')
  }
  
  x_out <- x %>%
    #Decide the direction that metabolites are moving 
    mutate(direction = ifelse(effect_size > 0, 'up', 'down')) %>%
    #Non-significant metabolites are not moving
    mutate(direction = ifelse(pvalue > 0.05, 'no_change', direction)) %>%
    #Conver effect size for downward moving metabolites into 1/efs
    mutate(efs = ifelse(effect_size < 0, 1/abs(effect_size), abs(effect_size))) %>%
    #Non-moving metabolites have no effect size
    mutate(efs = ifelse(pvalue > 0.05, 1, efs)) %>%
    #Drop the NAs from the dataset
    filter(!is.na(set)) %>%
    group_by(set) %>%
    #Remove any compound classes which have no significant differences
    filter(!all(pvalue > 0.05)) %>%
    #Remove any classes with less than 2 compounds
    filter(n() > 2) %>%
    nest() %>%
    #Run the KS-Test essentially asking the question, what is the likelihood the pvalue distribution is uniform for this metabolite?
    mutate(ptest = map_dbl(data, function(x) suppressWarnings(ks.test(x$pvalue, 'punif', alternative = 'greater')$p.value))) %>%
    #Correct p-values rounded to 0
    modify_at('ptest', ~ifelse(.x == 0, 2.2e-20, .x)) %>%
    #unnest('data') 
    #Calculate ration of significant to all
    mutate(altrat = map_dbl(data, function(x) sum(x$pvalue < 0.05)/length(x$pvalue))) %>%
    #Calculate the ratio of metabolites which were increased
    mutate(uprat = map_dbl(data, function(x) sum(x$pvalue < 0.05 & x$direction == 'up')/sum(x$pvalue < 0.05))) %>% 
    #Calculate the number of metaboltes in the group
    mutate(size = map_dbl(data, function(x) nrow(x))) %>%
    #Calculate the number of significant metabolties
    mutate(num_altered = map_dbl(data, function(x) sum(x$pvalue < 0.05))) %>%
    #How many went up?
    mutate(num_up = map_dbl(data, function(x) sum(x$direction == 'up' ))) %>%
    #How many went down?
    mutate(num_down = map_dbl(data, function(x) sum(x$direction == 'down' ))) %>%
    #Ungroup to prevent weird things from happengin
    ungroup() %>%
    #Now adjust p-values
    modify_at('ptest', ~p.adjust(.x, method = 'fdr')) %>%
    #Create a random order for chemical classes
    mutate(order = sample(1:nrow(.), replace = F))
  return(x_out)
}

suppressMW <- function(x){
  suppressMessages(suppressWarnings(x))
}

plotChemRich <- function(CR_output, interactive = FALSE, colors = color, nset = 10){
  p <- suppressMW(ggplot(CR_output, aes(x = order, y = -log(ptest), size = size, color = uprat)) +
                    geom_point(aes(text = paste0('Class: ', set, '\n',
                                                 'K-S P-Val: ', ptest, '\n',
                                                 'Up Expressed: ', num_up, '\n',
                                                 'Down Expressed: ', num_down, '\n',
                                                 'Size: ', size))) +
                    scale_color_gradient(low = colors[1], high = colors[2], limits = c(0,1)) +
                    scale_size(range = c(5,30)) +
                    ggrepel::geom_label_repel(aes(label = set), color = 'gray20', data = subset(CR_output, size > nset), size = 4, force = 5) +
                    theme(legend.position = 'right'))+
                    theme_linedraw()+
                    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank())+ 
                    scale_y_sqrt()
    
  if(interactive){
    suppressMW(plotly::ggplotly(p, tooltip = 'text'))
  }else{
    suppressMW(plot(p))
  }
} 

library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)

setwd("~/INVA/INVA")


performChemRICH <- function(name) {
  file <- paste0(name,"_all.csv")
  df <- read.csv(file, header = TRUE) %>%
    left_join(canopus_result) %>%
    dplyr::select(compoundID, Log2FC, pvalue, ClassyFire.subclass) %>%
    rename(compound_name = compoundID, 
           effect_size = Log2FC, 
           pvalue = pvalue, 
           set = ClassyFire.subclass)
  
  NMR <- df[grepl("NMR_", df$compound_name), ]
  df <- df[!grepl("NMR_", df$compound_name), ]
  NMR$compound_name <- gsub("X(\\d+)\\.", "\\1-", NMR$compound_name)
  NMR$compound_name <- gsub("\\.", " ", NMR$compound_name)
  NMR$compound_name <- gsub(" ", "-", NMR$compound_name)
  
  NMR_chem <- read.csv("NMR_chem.csv", header = TRUE)
  NMR_chem$compound_name <- paste0("NMR_", NMR_chem$metabolite)
  NMR_chem$compound_name <- gsub(" ", "-", NMR_chem$compound_name)
  
  NMR <- left_join(NMR, NMR_chem)
  NMR$set <- NMR$Subclass
  NMR <- NMR[, c(1:4)]
  
  df <- rbind(df, NMR)
  df <- df[df$set != "", ]
  df <- df[!is.na(df$set), ]
  
  result <- chemRICHFix(df)
  assign(paste0(name, "_result"), result, envir = .GlobalEnv)
  plotChemRich(result)
}


canopus_result<-read.table("/Users/benyang/INVA/SIRIUS_results_root/canopus_compound_summary.tsv", header = T, sep = "\t", fill = TRUE, quote = "")
canopus_result$ClassyFire.subclass <- gsub('"', '', canopus_result$ClassyFire.subclass)
color <- c("#50124C", "#4DAF4A")

performChemRICH("root_ode")
performChemRICH("root_cde")

color <-c("#E49518", "#009CA3")
ggsave("enrich_root_doc.pdf",performChemRICH("root_doc"), width = 8, height = 4)
ggsave("enrich_root_eoc.pdf",performChemRICH("root_eoc"), width = 8, height = 4)


canopus_result<-read.table("/Users/benyang/INVA/SIRIUS_results_rhizo/canopus_compound_summary.tsv", header = T, sep = "\t", fill = TRUE, quote = "")
canopus_result$ClassyFire.subclass <- gsub('"', '', canopus_result$ClassyFire.subclass)
color <- c("#50124C", "#4DAF4A")
performChemRICH("rhizo_ode")
performChemRICH("rhizo_cde")

color <-c("#E49518", "#009CA3")
performChemRICH("rhizo_doc")
performChemRICH("rhizo_eoc")

# try to plot volcano plot on shared features only.
shared <- read.csv("shared_Rhizo_compounds.csv", header = T)

#doc
name<- "rhizo_doc"
file <- paste0(name,"_all.csv")
df <- read.csv(file, header = TRUE) 

NMR <- df[grepl("NMR_", df$compoundID), ]
NMR$compoundID <- gsub("X(\\d+)\\.", "\\1-", NMR$compoundID)
NMR$compoundID <- gsub("\\.", " ", NMR$compoundID)
NMR$compoundID <- gsub(" ", "-", NMR$compoundID)

NMR_chem <- read.csv("NMR_chem.csv")
NMR_all_root <- read.table("NMR_root_clean.txt", header = TRUE, sep = "\t") %>%
  .$metabolite
NMR_all_rhizo <- read.table("NMR_rhizo_clean.txt", header = TRUE, sep = "\t") %>%
  .$metabolite

NMR.shared <- subset(NMR_all_rhizo, NMR_all_rhizo %in% NMR_all_root)
NMR_chem <- subset(NMR_chem, NMR_chem$metabolite %in% NMR.shared)

NMR_chem$compoundID <- paste0("NMR_", NMR_chem$metabolite)
NMR_chem$compoundID <- gsub(" ", "-", NMR_chem$compoundID)

NMR <- subset(NMR, NMR$compoundID %in% NMR_chem$compoundID)
NMR <- left_join(NMR, NMR_chem, by = "compoundID")

NMR <- NMR[, c(1,2,5,11,18,19)]
df.shared <- subset(df, df$compoundID %in% shared$compoundID)
df.shared <- df.shared[, c(1,2,5,11,9,10)]

colnames(NMR) <- colnames(df.shared)
df <- rbind(df.shared, NMR)


color_palette <- ggpubr::get_palette(palette = 'Dark2', 6)
maingroups<-c("Amino acids, peptides, and analogues","Carbohydrates and carbohydrate conjugates","Carboxylic acids and derivatives","Fatty acids and conjugates","Lineolic acids and derivatives","Hydroxy acids and derivatives")
colors <- setNames(color_palette,maingroups)

carboacid <- unique(df[which(df$Class == "Carboxylic acids and derivatives"),]$Subclass)
hyrdoacid <- unique(df[which(df$Class == "Hydroxy acids and derivatives"),]$Subclass)

carboacid <- carboacid[-1]

# Create a copy of the input dataframe

df <- df %>%
  mutate(Subclass = ifelse(Subclass %in% carboacid,
                        "Carboxylic acids and derivatives",
                        Subclass))

df <- df %>%
  mutate(Subclass = ifelse(Subclass %in% hyrdoacid,
                        "Hydroxy acids and derivatives",
                        Subclass))

df <- df %>%
  mutate(group = case_when(
    Subclass %in% maingroups ~ Subclass,
    TRUE ~ "other"))

df[which(df$pvalue > 0.05 | (df$Log2FC < 1.5 & df$Log2FC > -1.5)),]$group <- "NO"
df$label <- ifelse(df$Tags == "Confirmed ID (HIgh Confidence)" & df$group != "NO", df$Name, NA)
df$label <- gsub("DL-","", df$label)


# Create volcano plot
plot_volcano_doc <- ggplot(data = df, aes(x = Log2FC, y = -log10(pvalue), col = group)) + 
  geom_vline(xintercept=c(-1.5, 1.5), col="dimgrey") +
  geom_hline(yintercept=-log10(0.05), col="dimgrey") +
  geom_point(size = 1) + 
  theme_minimal() +
  theme(legend.position = "none")+
  scale_color_manual(values = c(colors, "other" = "bisque2", "NO" = "grey"))

##eoc
name<- "rhizo_eoc"
file <- paste0(name,"_all.csv")
df <- read.csv(file, header = TRUE) 

NMR <- df[grepl("NMR_", df$compoundID), ]
NMR$compoundID <- gsub("X(\\d+)\\.", "\\1-", NMR$compoundID)
NMR$compoundID <- gsub("\\.", " ", NMR$compoundID)
NMR$compoundID <- gsub(" ", "-", NMR$compoundID)

NMR_chem <- read.csv("NMR_chem.csv")
NMR_all_root <- read.table("NMR_root_clean.txt", header = TRUE, sep = "\t") %>%
  .$metabolite
NMR_all_rhizo <- read.table("NMR_rhizo_clean.txt", header = TRUE, sep = "\t") %>%
  .$metabolite

NMR.shared <- subset(NMR_all_rhizo, NMR_all_rhizo %in% NMR_all_root)
NMR_chem <- subset(NMR_chem, NMR_chem$metabolite %in% NMR.shared)

NMR_chem$compoundID <- paste0("NMR_", NMR_chem$metabolite)
NMR_chem$compoundID <- gsub(" ", "-", NMR_chem$compoundID)

NMR <- subset(NMR, NMR$compoundID %in% NMR_chem$compoundID)
NMR <- left_join(NMR, NMR_chem, by = "compoundID")

NMR <- NMR[, c(1,2,5,11,18,19)]
df.shared <- subset(df, df$compoundID %in% shared$compoundID)
df.shared <- df.shared[, c(1,2,5,11,9,10)]

colnames(NMR) <- colnames(df.shared)
df <- rbind(df.shared, NMR)


color_palette <- ggpubr::get_palette(palette = 'Dark2', 6)
maingroups<-c("Amino acids, peptides, and analogues","Carbohydrates and carbohydrate conjugates","Carboxylic acids and derivatives","Fatty acids and conjugates","Lineolic acids and derivatives","Hydroxy acids and derivatives")
colors <- setNames(color_palette,maingroups)

carboacid <- unique(df[which(df$Class == "Carboxylic acids and derivatives"),]$Subclass)
hyrdoacid <- unique(df[which(df$Class == "Hydroxy acids and derivatives"),]$Subclass)

carboacid <- carboacid[-1]

# Create a copy of the input dataframe

df <- df %>%
  mutate(Subclass = ifelse(Subclass %in% carboacid,
                           "Carboxylic acids and derivatives",
                           Subclass))

df <- df %>%
  mutate(Subclass = ifelse(Subclass %in% hyrdoacid,
                           "Hydroxy acids and derivatives",
                           Subclass))

df <- df %>%
  mutate(group = case_when(
    Subclass %in% maingroups ~ Subclass,
    TRUE ~ "other"))

df[which(df$pvalue > 0.05 | (df$Log2FC < 1.5 & df$Log2FC > -1.5)),]$group <- "NO"
df$label <- ifelse(df$Tags == "Confirmed ID (HIgh Confidence)" & df$group != "NO", df$Name, NA)
df$label <- gsub("DL-","", df$label)

# Create volcano plot
plot_volcano_eoc <- ggplot(data = df, aes(x = Log2FC, y = -log10(pvalue), col = group)) + 
  geom_vline(xintercept=c(-1.5, 1.5), col="dimgrey") +
  geom_hline(yintercept=-log10(0.05), col="dimgrey") +
  geom_point(size = 1) + 
  theme_minimal() +
  theme(legend.position = "none")+
  scale_color_manual(values = c(colors, "other" = "bisque2", "NO" = "grey"))

##plant_all
name<- "rhizo_plant"
file <- paste0(name,"_all.csv")
df <- read.csv(file, header = TRUE) 

NMR <- df[grepl("NMR_", df$compoundID), ]
NMR$compoundID <- gsub("X(\\d+)\\.", "\\1-", NMR$compoundID)
NMR$compoundID <- gsub("\\.", " ", NMR$compoundID)
NMR$compoundID <- gsub(" ", "-", NMR$compoundID)

NMR_chem <- read.csv("NMR_chem.csv")
NMR_all_root <- read.table("NMR_root_clean.txt", header = TRUE, sep = "\t") %>%
  .$metabolite
NMR_all_rhizo <- read.table("NMR_rhizo_clean.txt", header = TRUE, sep = "\t") %>%
  .$metabolite

NMR.shared <- subset(NMR_all_rhizo, NMR_all_rhizo %in% NMR_all_root)
NMR_chem <- subset(NMR_chem, NMR_chem$metabolite %in% NMR.shared)

NMR_chem$compoundID <- paste0("NMR_", NMR_chem$metabolite)
NMR_chem$compoundID <- gsub(" ", "-", NMR_chem$compoundID)

NMR <- subset(NMR, NMR$compoundID %in% NMR_chem$compoundID)
NMR <- left_join(NMR, NMR_chem, by = "compoundID")

NMR <- NMR[, c(1,2,5,11,18,19)]
df.shared <- subset(df, df$compoundID %in% shared$compoundID)
df.shared <- df.shared[, c(1,2,5,11,9,10)]

colnames(NMR) <- colnames(df.shared)
df <- rbind(df.shared, NMR)


color_palette <- ggpubr::get_palette(palette = 'Dark2', 6)
maingroups<-c("Amino acids, peptides, and analogues","Carbohydrates and carbohydrate conjugates","Carboxylic acids and derivatives","Fatty acids and conjugates","Lineolic acids and derivatives","Hydroxy acids and derivatives")
colors <- setNames(color_palette,maingroups)

carboacid <- unique(df[which(df$Class == "Carboxylic acids and derivatives"),]$Subclass)
hyrdoacid <- unique(df[which(df$Class == "Hydroxy acids and derivatives"),]$Subclass)

carboacid <- carboacid[-1]

# Create a copy of the input dataframe

df <- df %>%
  mutate(Subclass = ifelse(Subclass %in% carboacid,
                           "Carboxylic acids and derivatives",
                           Subclass))

df <- df %>%
  mutate(Subclass = ifelse(Subclass %in% hyrdoacid,
                           "Hydroxy acids and derivatives",
                           Subclass))

df <- df %>%
  mutate(group = case_when(
    Subclass %in% maingroups ~ Subclass,
    TRUE ~ "other"))

df[which(df$pvalue > 0.05 | (df$Log2FC < 1.5 & df$Log2FC > -1.5)),]$group <- "NO"
df$label <- ifelse(df$Tags == "Confirmed ID (HIgh Confidence)" & df$group != "NO", df$Name, NA)
df$label <- gsub("DL-","", df$label)

library(ggrepel)
# Create volcano plot
plot_volcano_all <- ggplot(data = df, aes(x = Log2FC, y = -log10(pvalue), col = group)) + 
  geom_vline(xintercept=c(-1.5, 1.5), col="dimgrey") +
  geom_hline(yintercept=-log10(0.05), col="dimgrey") +
  geom_point(size = 1) + 
  theme_minimal() +
  theme(legend.position = "none")+
  scale_color_manual(values = c(colors, "other" = "bisque2", "NO" = "grey"))

ggsave("volcano_exudates_plant.pdf",plot_volcano_all,width =5, height = 4)
ggsave("volcano_exudates_eoc.pdf",plot_volcano_eoc,width =5, height = 4)
ggsave("volcano_exudates_doc.pdf",plot_volcano_doc,width =5, height = 4)


# Try chemrich on rhizo seprately.

canopus_result<-read.table("/Users/benyang/INVA/SIRIUS_results_rhizo/canopus_compound_summary.tsv", header = T, sep = "\t", fill = TRUE, quote = "")
canopus_result$ClassyFire.subclass <- gsub('"', '', canopus_result$ClassyFire.subclass)
color <- c("#50124C", "#4DAF4A")
# this function only give results df but not plotting
performChemRICH <- function(name) {
  file <- paste0(name,"_all.csv")
  df <- read.csv(file, header = TRUE) %>%
    left_join(canopus_result) %>%
    dplyr::select(compoundID, Log2FC, pvalue, ClassyFire.subclass) %>%
    rename(compound_name = compoundID, 
           effect_size = Log2FC, 
           pvalue = pvalue, 
           set = ClassyFire.subclass)
  
  NMR <- df[grepl("NMR_", df$compound_name), ]
  df <- df[!grepl("NMR_", df$compound_name), ]
  NMR$compound_name <- gsub("X(\\d+)\\.", "\\1-", NMR$compound_name)
  NMR$compound_name <- gsub("\\.", " ", NMR$compound_name)
  NMR$compound_name <- gsub(" ", "-", NMR$compound_name)
  
  NMR_chem <- read.csv("NMR_chem.csv", header = TRUE)
  NMR_chem$compound_name <- paste0("NMR_", NMR_chem$metabolite)
  NMR_chem$compound_name <- gsub(" ", "-", NMR_chem$compound_name)
  
  NMR <- left_join(NMR, NMR_chem)
  NMR$set <- NMR$Subclass
  NMR <- NMR[, c(1:4)]
  
  df <- rbind(df, NMR)
  df <- df[df$set != "", ]
  df <- df[!is.na(df$set), ]
  
  result <- chemRICHFix(df)
  assign(paste0(name, "_result"), result, envir = .GlobalEnv)
}
performChemRICH("rhizo_ode")
performChemRICH("rhizo_cde")

plotChemRich(rhizo_cde_result)
plotChemRich(rhizo_ode_result)

library(reshape2)
cde_20 <- rhizo_cde_result[order(rhizo_cde_result$size, decreasing = TRUE), ][1:20, ]
cde_20 <- cde_20[order(cde_20$set), ]
rank <- cde_20[order(cde_20$size),]$set
cde_20$un_altered <- cde_20$size-cde_20$num_altered
cde_20_long <- cde_20[,c(1,6, 8,9,11)]
cde_20_long<- melt(cde_20_long, id.vars = c("set","size"), variable.name = "identity")
cde_20_long$identity <- factor(cde_20_long$identity, levels = c("un_altered", "num_up", "num_down"))
cde_20_long$set <- factor(cde_20_long$set, levels = rank)

ode_20 <- rhizo_ode_result[order(rhizo_ode_result$size, decreasing = TRUE), ][1:20, ]
ode_20 <- ode_20[order(ode_20$set), ]
rank <- ode_20[order(ode_20$size),]$set
ode_20$un_altered <- ode_20$size-ode_20$num_altered
ode_20_long <- ode_20[,c(1,6, 8,9,11)]
ode_20_long<- melt(ode_20_long, id.vars = c("set", "size"), variable.name = "identity")
ode_20_long$identity <- factor(ode_20_long$identity, levels = c("un_altered", "num_up", "num_down"))
ode_20_long$set <- factor(ode_20_long$set, levels = rank)

plotChemRich <- function(CR_output){
  p <- suppressMW(ggplot(CR_output, aes(x = value, y = set, fill = identity)) +
                    geom_bar(position = "stack", stat = "identity") +
                    scale_fill_manual(values = c("num_down" = "#50124C", "num_up" = "#4DAF4A",  "un_altered" = "grey")) +
                    scale_x_continuous(limits = c(0, 401), breaks = seq(0, 400, by = 50)) + 
                    theme_minimal()+
                    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank())+
                    theme(axis.title = element_blank(), legend.position = "none"))
  return(p)
} 


ode_p10 <- subset(rhizo_ode_result, rhizo_ode_result$size > 10)
ode_p10 <- ode_p10[order(-ode_p10$ptest), ]
rank <- ode_p10[order(-ode_p10$ptest),]$set
ode_p10$set <- factor(ode_p10$set, levels = rank)

cde_p10 <- subset(rhizo_cde_result, rhizo_cde_result$size > 10)
cde_p10 <- cde_p10[order(-cde_p10$ptest), ]
cde_p10$set <- factor(cde_p10$set, levels = rank)


 plotChemRich <- function(CR_output){
  p <- suppressMW(ggplot(CR_output, aes(y = set, x = -log(ptest), size = size, color = uprat)) +
                    geom_point(aes(text = paste0('Class: ', set, '\n',
                                                 'K-S P-Val: ', ptest, '\n',
                                                 'Up Expressed: ', num_up, '\n',
                                                 'Down Expressed: ', num_down, '\n',
                                                 'Size: ', size))) +
                    scale_color_gradient(low = "#50124C", high = "#4DAF4A", limits = c(0,1)) +
                    scale_size(range = c(5,30)) +
                    xlab("-log(P)")+
                    ylab("Subclasses")+
                    theme_minimal()+
                    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank())+
                    theme(legend.position = "right"))
  return(p)
}

plotChemRich(cde_p10)
plotChemRich(ode_p10)

ggsave("rhizo_cde_enrich.pdf",plotChemRich(cde_20_long),height = 14, width = 6)
ggsave("rhizo_ode_enrich.pdf", plotChemRich(ode_20_long),height = 14, width = 6)


ggsave("rhizo_cde_enrich.pdf",plotChemRich(cde_p10),height = 6, width = 5)
ggsave("rhizo_ode_enrich.pdf", plotChemRich(ode_p10),height = 6, width = 5)

#MALDI
MALDI.open <- read.csv("Open_rings.csv")
MALDI.close <- read.csv("Canopy_rings.csv")

MALDI.open <- MALDI.open %>%
  dplyr::rename(
    compound_name = ID, 
    effect_size = Log2FC, 
    pvalue = pvalue, 
    set = Class
  )

MALDI.close <- MALDI.close %>%
  dplyr::rename(
    compound_name = ID, 
    effect_size = Log2FC, 
    pvalue = pvalue, 
    set = Class
  )

MALDI.open.chemrich<-chemRICHFix(MALDI.open)
color <- c("#50124C", "#4DAF4A")
MALDI.open.chemrich<-subset(MALDI.open.chemrich,MALDI.open.chemrich$ptest < 0.05)

rank <- MALDI.open.chemrich[order(-MALDI.open.chemrich$ptest),]$set
MALDI.open.chemrich$set <- factor(MALDI.open.chemrich$set, levels = rank)
ggsave("MALDI_open_enrich.pdf", plotChemRich(MALDI.open.chemrich), height = 8, width = 8)

filtered_data <- MALDI.open %>%
  filter(diffexpressed != "NO")

# Identify unique compound_name entries
unique_compounds <- filtered_data %>%
  group_by(compound_name) %>%
  filter(n() == 1) %>%
  ungroup()

MALDI.close.chemrich<-chemRICHFix(MALDI.close)
MALDI.close.chemrich<-subset(MALDI.close.chemrich,MALDI.close.chemrich$ptest < 0.05)
rank <- MALDI.close.chemrich[order(-MALDI.close.chemrich$ptest),]$set
MALDI.close.chemrich$set <- factor(MALDI.close.chemrich$set, levels = rank)
ggsave("MALDI_close_enrich.pdf", plotChemRich(MALDI.close.chemrich), height = 5.5, width = 8)
