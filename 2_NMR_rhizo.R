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
library(KEGGREST)

########## Load Functions ########
#Pareto Scaling:
PS_helper <- function(x) {
  (x - mean(x)) / sqrt(sd(x, na.rm = T))
}	

pareto_scale <- function(x){
  mtb_scaled <- data.frame(apply(x, 2, PS_helper))
  return(mtb_scaled)
}

### can be sensitive to large range. The difference between Z-score is it squar rhizo the sd again.

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
# fancy! 

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
NMR_rhizo <- as.data.frame(t(read.table("NMR_rhizo_clean.txt", header = TRUE, sep = "\t", row.names =1)))
NMR_chem <- read.csv("NMR_chem.csv", header = TRUE)
NMR_chem <- NMR_chem %>% 
  arrange(Class)

#normalization
boxplot(NMR_rhizo)

NMR_pareto <- NMR_rhizo %>%
  log_transform() %>%
  pareto_scale()

NMR_log <- NMR_rhizo %>%
  log_transform()

NMR_auto <- NMR_rhizo %>%
  log_transform()%>%
  auto_scale

boxplot(NMR_pareto)
boxplot(NMR_log)
boxplot(NMR_auto)

rownames(NMR_pareto) <- sub("\\.", "-", rownames(NMR_pareto))
rownames(NMR_pareto) <- metadata[rownames(NMR_pareto),] $SampleID_rhizo
NMR_pareto <- NMR_pareto[order(rownames(NMR_pareto)),]

######### Q1: how does open/close * inva/noninva differ in metabolites profile?#######

#use adonis and NMDS for quantitative results. - pareto
NMR_rhizo_dist <- vegdist(NMR_pareto,method="euclidean")
saveRDS(NMR_rhizo_dist, "rhizo_NMR_dis.rds", )

NMR_pareto <- NMR_pareto[metadata$SampleID_rhizo,]
NMR_rhizo_dist <- vegdist(NMR_pareto,method="euclidean")
NMR_rhizo.nmds <- metaMDS(NMR_rhizo_dist)

metadata$NMDS1 <- NMR_rhizo.nmds$points[,1]
metadata$NMDS2 <- NMR_rhizo.nmds$points[,2]

Plot_NMR_rhizo_pareto<-ggplot(metadata, aes(x = NMDS1, y = NMDS2, group = block))+
  geom_point(aes(shape = canopy, color = plant, fill = plant), size = 3)+
  theme_bw(base_size = 15)+
  scale_color_manual(values = c("ERLE" = "#50124C", "DICA" = "#4DAF4A"))+
  scale_fill_manual(values = c("ERLE" = "#50124C", "DICA" = "#4DAF4A"))+
  scale_shape_manual(values = c("close" = 16, "open" = 1))+
  labs(x = NULL, y = NULL) +
  theme(legend.position="none") 

permanova_NMR_rhizo_pareto <- adonis2(NMR_rhizo_dist ~ plant * canopy , data = metadata)
#plant 0.001 (0.30129)； canopy 0.001 (0.10554); plant: canopy 0.163 (0.02781)
NMR_rhizo.nmds$stress #0.1089611

permanova_table <- as.data.frame(permanova_NMR_rhizo_pareto)
permanova_table[, sapply(permanova_table, is.numeric)] <- round(permanova_table[, sapply(permanova_table, is.numeric)], 3)
table_plot <- tableGrob(permanova_table)
ggsave("NMDS_rhizo_NMR.pdf", Plot_NMR_rhizo_pareto, width = 4, height = 4)

DICA <- NMR_pareto[metadata[metadata$plant == "DICA",]$SampleID_rhizo,]
ERLE <- NMR_pareto[metadata[metadata$plant == "ERLE",]$SampleID_rhizo,]

DICA.dist <- vegdist(DICA,method="euclidean")
ERLE.dist <- vegdist(ERLE,method="euclidean")

set.seed(123)
adonis2(DICA.dist ~ canopy , data = metadata[metadata$plant=="DICA",])
adonis2(ERLE.dist ~ canopy , data = metadata[metadata$plant=="ERLE",])


########differential analysis -------
#From NMDS, it is clear now DICA/ERLE and open/canopy for DICA is the driving difference
#Use log2fold change to assess those differences.
# I can do paried analysis actually and make volcano plots.


NMR_final <- NMR_pareto %>%
  rownames_to_column('SampleID_rhizo') %>%
  right_join(metadata) %>%
  column_to_rownames('SampleID_rhizo')

p.plant <- apply(NMR_final[,1:ncol(NMR_rhizo)], 2, function(x){kruskal.test(x~plant, data = NMR_final)$p.value})

# calculate log2fold change
# Convert the row names of one dataframe to match the other
rownames(NMR_rhizo) <- gsub("\\.", "-", rownames(NMR_rhizo))  # Convert "-" to "." in row names
NMR_rhizo$RowNames <- rownames(NMR_rhizo)
metadata$RowNames <- rownames(metadata)
# Merge the dataframes by row names
rhizo_NMR_LC2 <- merge(NMR_rhizo, metadata, by = "RowNames")
rhizo_NMR_LC2 <- column_to_rownames(rhizo_NMR_LC2, var = "RowNames")


ERLE <- rhizo_NMR_LC2[which(rhizo_NMR_LC2$plant == "ERLE"),1:(ncol(NMR_rhizo)-1)]
DICA <- rhizo_NMR_LC2[which(rhizo_NMR_LC2$plant == "DICA"),1:(ncol(NMR_rhizo)-1)]
LC2_NMR <- log2(colMeans(DICA)/colMeans(ERLE))
rhizo_volcano_NMR <- data.frame(pvalue = p.plant,
                          Log2FC = LC2_NMR)  


rhizo_volcano_NMR$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "DICA" 
rhizo_volcano_NMR$diffexpressed[rhizo_volcano_NMR$Log2FC > 1.5 & rhizo_volcano_NMR$pvalue < 0.05] <- "DICA"
rhizo_volcano_NMR$diffexpressed[rhizo_volcano_NMR$Log2FC < -1.5 & rhizo_volcano_NMR$pvalue < 0.05] <- "ERLE"
rhizo_volcano_NMR$delabel <- NA
rhizo_volcano_NMR$delabel[rhizo_volcano_NMR$diffexpressed != "NO"] <- rownames(rhizo_volcano_NMR)[rhizo_volcano_NMR$diffexpressed != "NO"]

#assign taxonomy

rhizo_volcano_NMR$metabolite <- rownames(rhizo_volcano_NMR)
rhizo_volcano_NMR$metabolite <- gsub("X", "", rhizo_volcano_NMR$metabolite)
rhizo_volcano_NMR$metabolite <- gsub(".", "-", rhizo_volcano_NMR$metabolite, fixed = TRUE)
NMR_chem$metabolite <- gsub(" ", "-", NMR_chem$metabolite, fixed = TRUE)
rhizo_volcano_NMR <- left_join(rhizo_volcano_NMR, NMR_chem)

class_colors <- get_palette(palette = 'Set1', k = length(unique(NMR_chem$Subclass)))
names(class_colors) <- c(unique(NMR_chem$Subclass))

rhizo_volcano_NMR$group <- rhizo_volcano_NMR$Subclass
rhizo_volcano_NMR[which(rhizo_volcano_NMR$diffexpressed == "NO"),]$group <- "NO"


volcano_rhizo <- ggplot(data=rhizo_volcano_NMR, aes(x=Log2FC, y=-log10(pvalue), col=group, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = Inf) +
  scale_color_manual(values=c(class_colors, "NO" = "grey")) +
  geom_vline(xintercept=c(-1.5, 1.5), col="dimgrey") +
  geom_hline(yintercept=-log10(0.05), col="dimgrey")


ggsave("rhizo_volcano_NMR.pdf", volcano_rhizo + ggtitle("NMR plant") + theme(legend.position = "none"), width = 4, height = 4)
write.table(rhizo_volcano_NMR, file = "Diff_Rhizo_plant_NMR.csv", sep = ",", row.names = FALSE)

#### look at different responses of rhizo metabolites under close and open.
e <- NMR_final[NMR_final$plant == "ERLE" , ] 
d <- NMR_final[NMR_final$plant == "DICA" , ]

p_eoc <- apply(e[,c(1:length(NMR_rhizo)-1)],2, function(x){kruskal.test(x~canopy, data = e)$p.value})
p_doc <- apply(d[,c(1:length(NMR_rhizo)-1)],2, function(x){kruskal.test(x~canopy, data = d)$p.value})

LC2.eo <- rhizo_NMR_LC2[rhizo_NMR_LC2$plant == "ERLE" & rhizo_NMR_LC2$canopy == "open", 1:(ncol(NMR_rhizo)-1)]
LC2.ec <- rhizo_NMR_LC2[rhizo_NMR_LC2$plant == "ERLE" & rhizo_NMR_LC2$canopy == "close", 1:(ncol(NMR_rhizo)-1)]
LC2.do <- rhizo_NMR_LC2[rhizo_NMR_LC2$plant == "DICA" & rhizo_NMR_LC2$canopy == "open", 1:(ncol(NMR_rhizo)-1)]
LC2.dc <- rhizo_NMR_LC2[rhizo_NMR_LC2$plant == "DICA" & rhizo_NMR_LC2$canopy == "close", 1:(ncol(NMR_rhizo)-1)]

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

#assign taxonomy
volcano_eoc$metabolite <- rownames(volcano_eoc)
volcano_eoc$metabolite <- gsub("X", "", volcano_eoc$metabolite)
volcano_eoc$metabolite <- gsub(".", "-", volcano_eoc$metabolite, fixed = TRUE)
NMR_chem$metabolite <- gsub(" ", "-", NMR_chem$metabolite, fixed = TRUE)
volcano_eoc <- left_join(volcano_eoc, NMR_chem)
volcano_eoc$group <- volcano_eoc$Subclass
volcano_eoc[which(volcano_eoc$diffexpressed == "NO"),]$group <- "NO"


volcano_ERLE_oc <- ggplot(data=volcano_eoc, aes(x=Log2FC, y=-log10(pvalue), col=group, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = Inf) +
  scale_color_manual(values=c(class_colors, "NO" = "grey")) +
  geom_vline(xintercept=c(-1.5, 1.5), col="dimgrey") +
  geom_hline(yintercept=-log10(0.05), col="dimgrey")


volcano_doc$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "close" 
volcano_doc$diffexpressed[volcano_doc$Log2FC > 1.5 & volcano_doc$pvalue < 0.05] <- "CLOSE"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "open"
volcano_doc$diffexpressed[volcano_doc$Log2FC < -1.5 & volcano_doc$pvalue < 0.05] <- "OPEN"
volcano_doc$delabel <- NA
volcano_doc$delabel[volcano_doc$diffexpressed != "NO"] <- rownames(volcano_doc)[volcano_doc$diffexpressed != "NO"]
# doc has xylose in close

volcano_doc$metabolite <- rownames(volcano_doc)
volcano_doc$metabolite <- gsub("X", "", volcano_doc$metabolite)
volcano_doc$metabolite <- gsub(".", "-", volcano_doc$metabolite, fixed = TRUE)
NMR_chem$metabolite <- gsub(" ", "-", NMR_chem$metabolite, fixed = TRUE)
volcano_doc <- left_join(volcano_doc, NMR_chem)
volcano_doc$group <- volcano_doc$Subclass
volcano_doc[which(volcano_doc$diffexpressed == "NO"),]$group <- "NO"

volcano_DICA_oc <- ggplot(data=volcano_doc, aes(x=Log2FC, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = Inf) +
  scale_color_manual(values=c(class_colors, "NO" = "grey")) +
  geom_vline(xintercept=c(-1.5, 1.5), col="dimgrey") +
  geom_hline(yintercept=-log10(0.05), col="dimgrey")


ggsave("rhizo_volcano_eoc_NMR.pdf", volcano_ERLE_oc + ggtitle("ERLE open vs close") + theme(legend.position = "none"), width = 4, height = 4)
ggsave("rhizo_volcano_doc_NMR.pdf", volcano_DICA_oc + ggtitle("DICA open vs close") + theme(legend.position = "none"), width = 4, height = 4)

write.table(volcano_doc, file = "Diff_Rhizo_DICA_CvO_NMR.csv", sep = ",", row.names = FALSE)
write.table(volcano_eoc, file = "Diff_Rhizo_ERLE_CvO_NMR.csv", sep = ",", row.names = FALSE)


#### look at different responses of rhizo metabolites under open ERLE vs DICA close ERLE vs DICA
o <- NMR_final[NMR_final$canopy == "open" , ] 
c <- NMR_final[NMR_final$canopy == "close" , ]

p_ode <- apply(o[,c(1:length(NMR_rhizo)-1)],2, function(x){kruskal.test(x~plant, data = o)$p.value})
p_cde<- apply(c[,c(1:length(NMR_rhizo)-1)],2, function(x){kruskal.test(x~plant, data = c)$p.value})

LC2.ode <- log2(colMeans(LC2.do)/colMeans(LC2.eo))
LC2.cde <- log2(colMeans(LC2.dc)/colMeans(LC2.ec))

volcano_ode <- data.frame(pvalue = p_ode,
                          Log2FC = LC2.ode) 
volcano_cde <- data.frame(pvalue = p_cde,
                          Log2FC = LC2.cde) 


volcano_ode$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "close" 
volcano_ode$diffexpressed[volcano_ode$Log2FC > 1.5 & volcano_ode$pvalue < 0.05] <- "DICA"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "open"
volcano_ode$diffexpressed[volcano_ode$Log2FC < -1.5 & volcano_ode$pvalue < 0.05] <- "ERLE"
volcano_ode$delabel <- NA
volcano_ode$delabel[volcano_ode$diffexpressed != "NO"] <- rownames(volcano_ode)[volcano_ode$diffexpressed != "NO"]
# ode has nothing

#assign taxonomy
volcano_ode$delabel <- gsub("X", "", volcano_ode$delabel)
volcano_ode$delabel <- gsub(".", "-", volcano_ode$delabel, fixed = TRUE)

volcano_ode$metabolite <- rownames(volcano_ode)
volcano_ode <- left_join(volcano_ode, NMR_chem)

volcano_open_de <- ggplot(data=volcano_ode, aes(x=Log2FC, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = Inf) +
  scale_color_manual(values=c(class_colors, "NO" = "grey")) +
  geom_vline(xintercept=c(-1.5, 1.5), col="dimgrey") +
  geom_hline(yintercept=-log10(0.05), col="dimgrey")


volcano_cde$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "close" 
volcano_cde$diffexpressed[volcano_cde$Log2FC > 1.5 & volcano_cde$pvalue < 0.05] <-  "DICA"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "open"
volcano_cde$diffexpressed[volcano_cde$Log2FC < -1.5 & volcano_cde$pvalue < 0.05] <- "ERLE"
volcano_cde$delabel <- NA
volcano_cde$delabel[volcano_cde$diffexpressed != "NO"] <- rownames(volcano_cde)[volcano_cde$diffexpressed != "NO"]
# cde has xylose in close

volcano_cde$delabel <- gsub("X", "", volcano_cde$delabel)
volcano_cde$delabel <- gsub(".", "-", volcano_cde$delabel, fixed = TRUE)

volcano_cde$metabolite <- rownames(volcano_cde)
volcano_cde <- left_join(volcano_cde, NMR_chem)


volcano_close_de <- ggplot(data=volcano_cde, aes(x=Log2FC, y=-log10(pvalue), col=diffexpressed, label=delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_text_repel(max.overlaps = Inf) +
  scale_color_manual(values=c(class_colors, "NO" = "grey")) +
  geom_vline(xintercept=c(-1.5, 1.5), col="dimgrey") +
  geom_hline(yintercept=-log10(0.05), col="dimgrey")


ggsave("rhizo_volcano_eoc_NMR.pdf", volcano_ERLE_oc + ggtitle("ERLE open vs close") + theme(legend.position = "none"), width = 4, height = 4)
ggsave("rhizo_volcano_doc_NMR.pdf", volcano_DICA_oc + ggtitle("DICA open vs close") + theme(legend.position = "none"), width = 4, height = 4)

write.table(volcano_ode, file = "Diff_Rhizo_OPEN_DvE_NMR.csv", sep = ",", row.names = FALSE)
write.table(volcano_cde, file = "Diff_Rhizo_CLOSE_DvE_NMR.csv", sep = ",", row.names = FALSE)



NMR_raw <- left_join(NMR_rhizo, metadata)



P <- function(compound_name) {
  compound_name <- substitute(compound_name)
  ggplot(NMR_raw, aes(x = plant, y = !!compound_name)) +
    geom_jitter(aes(color = plant, shape = canopy), position = position_jitterdodge(dodge.width = 0.75)) +
    geom_boxplot(aes(color = plant, fill = plant, shape = canopy, alpha = canopy), outlier.shape = NA) +
    scale_fill_manual(values = c("ERLE" = "#3070AD", "DICA" = "#4DAF4A")) +
    scale_color_manual(values = c("ERLE" = "#3070AD", "DICA" = "#4DAF4A")) +
    scale_shape_manual(values = c("close" = 16, "open" = 1)) +
    scale_alpha_manual(values = c("close" = 0.7, "open" = 0)) +
    labs(title = paste("NMR", "rhizo", as.character(compound_name))) +
    xlab(NULL) +
    ylab("NMR_intensity") +
    theme_bw(base_size = 15) +
    theme(legend.position = "none")
}


P(Aspartate)

#Asparatate, Lactate, threonine get upregulated in DICA in both root and rhizo

