########## library #########
library(dplyr)
library(vegan)
library(ggplot2)
library(ggrepel)
library(stats)
library(factoextra)
library(gridExtra)
library(tibble)
library(performance)
library(lme4)

########## Load Functions ########
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

#Coefficient of Variation
cv <- function(x){
  sd(x)/mean(x)
}


########## data prep #######
setwd("~/INVA/INVA")
metadata <- read.table("metadata.txt", header = TRUE, sep = "\t", row.names = 1)
metadata$block <- as.factor(metadata$block)
NMR_root <- as.data.frame(t(read.table("NMR_root_clean.txt", header = TRUE, sep = "\t", row.names =1)))
#normalization
boxplot(NMR_root)

NMR_pareto <- NMR_root %>%
  log_transform() %>%
  pareto_scale()

NMR_log <- NMR_root %>%
  log_transform()

NMR_auto <- NMR_root %>%
  log_transform()%>%
  auto_scale

boxplot(NMR_pareto)
boxplot(NMR_log)
boxplot(NMR_auto)

rownames(NMR_pareto) <- sub("\\.", "-", rownames(NMR_pareto))
rownames(NMR_pareto) <- metadata[rownames(NMR_pareto),] $SampleID_root
NMR_pareto <- NMR_pareto[order(rownames(NMR_pareto)),]

######### Q1: how does open/close * inva/noninva differ in metabolites profile?#######

#use adonis and NMDS for quantitative results. - pareto
NMR_root_dist <- vegdist(NMR_pareto,method="euclidean")
saveRDS(NMR_root_dist, "Root_NMR_dis.rds", )
a <- readRDS("Root_NMR_dis.rds") #save distance matrix for correlation 

NMR_pareto <- NMR_pareto[metadata$SampleID_root,]
NMR_root_dist <- vegdist(NMR_pareto,method="euclidean")
NMR_root.nmds <- metaMDS(NMR_root_dist)

metadata$NMDS1 <- NMR_root.nmds$points[,1]
metadata$NMDS2 <- NMR_root.nmds$points[,2]

Plot_NMR_root_pareto<-ggplot(metadata, aes(x = NMDS1, y = NMDS2, group = block))+
  geom_point(aes(shape = canopy, color = plant, fill = plant), size = 3)+
  theme_bw(base_size = 15)+
  scale_color_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A"))+
  scale_fill_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A"))+
  scale_shape_manual(values = c("close" = 16, "open" = 1))+
  labs(x = NULL, y = NULL) + 
  annotate("text", x = Inf, y = -Inf, vjust = 0, hjust = 1, label = paste("Stress =", round(as.numeric(NMR_root$stress, 3))))+
  theme(legend.position="none") 

permanova_NMR_root_pareto <- adonis2(NMR_root_dist ~ plant * canopy , data = metadata)
#plant 0.001 (0.15615)； canopy 0.045 (0.04178); plant: canopy 0.012 (0.05609)
NMR_root.nmds$stress #0.1701989

permanova_table <- as.data.frame(permanova_NMR_root_pareto)
permanova_table[, sapply(permanova_table, is.numeric)] <- round(permanova_table[, sapply(permanova_table, is.numeric)], 3)
table_plot <- tableGrob(permanova_table)
ggsave("NMDS_Root_NMR.pdf", grid.arrange(Plot_NMR_root_pareto + ggtitle("Root_NMR"), table_plot, ncol = 1, heights = c(3, 1)), width = 5, height = 7)



#From NMDS, it is clear now DICA/ERLE and open/canopy for DICA is the driving difference
#Use log2fold change to assess those differences.
# I can do paried analysis actually and make volcano plots.


NMR_final <- NMR_pareto %>%
  rownames_to_column('SampleID_root') %>%
  right_join(metadata) %>%
  column_to_rownames('SampleID_root')

lm_list <- apply(NMR_final[, 1:ncol(NMR_pareto)], 2, function(x) {
  lm(x ~ plant * canopy, data = NMR_final, na.action = na.exclude)
})

#calculate p value
p_values_list <- lapply(lm_list, function(model) {
  anova_results <- anova(model)[c(1:3), 5]
  return(anova_results)
})
#export p results
root_NMR_results<- data.frame(
  p_values_plant = sapply(p_values_list, `[`, 1),
  p_values_canopy = sapply(p_values_list, `[`, 2),
  p_values_interaction = sapply(p_values_list, `[`, 3),
  stringsAsFactors = FALSE
)


# calculate log2fold change
# Convert the row names of one dataframe to match the other
rownames(NMR_root) <- gsub("\\.", "-", rownames(NMR_root))  # Convert "-" to "." in row names
NMR_root$RowNames <- rownames(NMR_root)
metadata$RowNames <- rownames(metadata)
# Merge the dataframes by row names
Root_NMR_LC2 <- merge(NMR_root, metadata, by = "RowNames")
Root_NMR_LC2 <- column_to_rownames(Root_NMR_LC2, var = "RowNames")

Root_NMR_LC2 %>%
  ggplot(aes(x = plant, y = Asparagine)) +
  geom_jitter(aes(color = plant, shape = canopy), position = position_jitterdodge(dodge.width = 0.75)) +
  geom_boxplot(aes(color = plant, fill = plant, shape = canopy, alpha = canopy), outlier.shape = NA) +
  scale_fill_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_color_manual(values = c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A")) +
  scale_shape_manual(values = c("close" = 16, "open" = 1)) +
  scale_alpha_manual(values = c("close" = 0.7, "open" = 0)) +
  labs(title = "Root_asparagine") +
  xlab(NULL) +
  ylab("NMR_intensity") +
  theme_bw(base_size = 15) +
  theme(legend.position = "none")




ERLE <- Root_NMR_LC2[which(Root_NMR_LC2$plant == "ERLE"),1:(ncol(NMR_root)-1)]
DICA <- Root_NMR_LC2[which(Root_NMR_LC2$plant == "DICA"),1:(ncol(NMR_root)-1)]
OPEN <- Root_NMR_LC2[which(Root_NMR_LC2$canopy == "open"),1:(ncol(NMR_root)-1)]
CLOSE <- Root_NMR_LC2[which(Root_NMR_LC2$canopy == "close"),1:(ncol(NMR_root)-1)]


LC2_NMR_oc <- log2(colMeans(CLOSE)/colMeans(OPEN))
LC2_NMR <- log2(colMeans(DICA)/colMeans(ERLE))
root_volcano_oc <- data.frame(pvalue = root_NMR_results$p_values_canopy,
                               Log2FC = LC2_NMR_oc)  
root_volcano_NMR <- data.frame(pvalue = root_NMR_results$p_values_plant,
                                Log2FC = LC2_NMR)  


root_volcano_NMR$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "DICA" 
root_volcano_NMR$diffexpressed[root_volcano_NMR$Log2FC > 1.5 & root_volcano_NMR$pvalue < 0.05] <- "DICA"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
root_volcano_NMR$diffexpressed[root_volcano_NMR$Log2FC < -1.5 & root_volcano_NMR$pvalue < 0.05] <- "ERLE"
root_volcano_NMR$delabel <- NA
root_volcano_NMR$delabel[root_volcano_NMR$diffexpressed != "NO"] <- rownames(root_volcano_NMR)[root_volcano_NMR$diffexpressed != "NO"]

volcano_root <- ggplot(data=root_volcano_NMR, aes(x=Log2FC, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("ERLE" = "#E41A1C", "DICA" = "#4DAF4A", "NO" = "grey")) +
  geom_vline(xintercept=c(-1.5, 1.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")


ggsave("root_volcano_NMR.pdf", volcano_root + ggtitle("NMR plant") + theme(legend.position = "none"), width = 6, height = 4)


#### look at different responses of root metabolites under close and open.
t.eo <- NMR_final[NMR_final$plant == "ERLE" & NMR_final$canopy == "open", 1:(ncol(NMR_root)-1)]
t.ec <- NMR_final[NMR_final$plant == "ERLE" & NMR_final$canopy == "close", 1:(ncol(NMR_root)-1)]
t.do <- NMR_final[NMR_final$plant == "DICA" & NMR_final$canopy == "open", 1:(ncol(NMR_root)-1)]
t.dc <- NMR_final[NMR_final$plant == "DICA" & NMR_final$canopy == "close", 1:(ncol(NMR_root)-1)]

p_eoc <- sapply(1:ncol(t.eo), function(i) {
  t.test(t.eo[, i], t.ec[, i])$p.value
})
p_doc <- sapply(1:ncol(t.do), function(i) {
  t.test(t.do[, i], t.dc[, i])$p.value
})

LC2.eo <- Root_NMR_LC2[Root_NMR_LC2$plant == "ERLE" & Root_NMR_LC2$canopy == "open", 1:(ncol(NMR_root)-1)]
LC2.ec <- Root_NMR_LC2[Root_NMR_LC2$plant == "ERLE" & Root_NMR_LC2$canopy == "close", 1:(ncol(NMR_root)-1)]
LC2.do <- Root_NMR_LC2[Root_NMR_LC2$plant == "DICA" & Root_NMR_LC2$canopy == "open", 1:(ncol(NMR_root)-1)]
LC2.dc <- Root_NMR_LC2[Root_NMR_LC2$plant == "DICA" & Root_NMR_LC2$canopy == "close", 1:(ncol(NMR_root)-1)]

LC2.eoc <- log2(colMeans(LC2.ec)/colMeans(LC2.eo))
LC2.doc <- log2(colMeans(LC2.dc)/colMeans(LC2.do))

volcano_eoc <- data.frame(pvalue = p_eoc,
                          Log2FC = LC2.eoc) 
volcano_doc <- data.frame(pvalue = p_doc,
                          Log2FC = LC2.doc) 


volcano_eoc$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "close" 
volcano_eoc$diffexpressed[volcano_eoc$Log2FC > 1.5 & volcano_eoc$pvalue < 0.05] <- "CLOSE"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "open"
volcano_eoc$diffexpressed[volcano_eoc$Log2FC < -1.5 & volcano_eoc$pvalue < 0.05] <- "OPEN"
volcano_eoc$delabel <- NA
volcano_eoc$delabel[volcano_eoc$diffexpressed != "NO"] <- rownames(volcano_eoc)[volcano_eoc$diffexpressed != "NO"]
# eoc has nothing

volcano_ERLE_oc <- ggplot(data=volcano_eoc, aes(x=Log2FC, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("OPEN" = "#1B9E77", "CLOSE" = "#7570B3", "NO" = "grey")) +
  geom_vline(xintercept=c(-1.5, 1.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")



volcano_doc$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "close" 
volcano_doc$diffexpressed[volcano_doc$Log2FC > 1.5 & volcano_doc$pvalue < 0.05] <- "CLOSE"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "open"
volcano_doc$diffexpressed[volcano_doc$Log2FC < -1.5 & volcano_doc$pvalue < 0.05] <- "OPEN"
volcano_doc$delabel <- NA
volcano_doc$delabel[volcano_doc$diffexpressed != "NO"] <- rownames(volcano_doc)[volcano_doc$diffexpressed != "NO"]


volcano_DICA_oc <- ggplot(data=volcano_doc, aes(x=Log2FC, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("OPEN" = "#1B9E77", "CLOSE" = "#7570B3", "NO" = "grey")) +
  geom_vline(xintercept=c(-1.5, 1.5), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")


ggsave("root_volcano_eoc_NMR.pdf", volcano_ERLE_oc + ggtitle("ERLE open vs close") + theme(legend.position = "none"), width = 6, height = 4)
ggsave("root_volcano_doc_NMR.pdf", volcano_DICA_oc + ggtitle("DICA open vs close") + theme(legend.position = "none"), width = 6, height = 4)

