# load library  -------------------------------------------------------------
library(tidyverse)
library(mixOmics)
library(gridExtra)
library(lme4)
library(lmerTest)
library(performance)
library(ggrepel)
library(vegan)
library(ggpubr)
setwd('~/INVA/INVA')
metadata <- read.table("metadata.txt", header = TRUE, sep = "\t", row.names = 1)
metadata$block <- as.factor(metadata$block)
Rhizo_compounds <- read.csv("Rhizo_compound_link.csv", header = TRUE)

#Load Functions ------------------------------------
source(file.path('~/INVA/functions.R'))

## Rp positive---------------------
rppos_raw <- read_csv("Rhizo_LCMS_RPPOS.csv") %>%
  #Make the compound IDs the row names
  column_to_rownames('Compound ID') %>%
  #Transpose the matrix to make it so samples are rows and compounds are columns
  t()  %>%
  as.data.frame()

rppos_pareto <- rppos_raw %>%
  log_transform() %>%
  pareto_scale()

#Combine metadata with metabolite data
Rhizo_rppos_final  <- rppos_pareto %>%
  rownames_to_column('SampleID_rhizo') %>%
  mutate(SampleID_rhizo = sub("_RP_Pos$", "", SampleID_rhizo)) %>% 
  right_join(metadata) %>%
  column_to_rownames('SampleID_rhizo')

#NMDS plot
Plot_nmds(Rhizo_rppos_final)


Rhizo_rppos_final_ERLE <- Rhizo_rppos_final[Rhizo_rppos_final$plant== "ERLE",]
Rhizo_rppos_final_DICA <- Rhizo_rppos_final[Rhizo_rppos_final$plant== "DICA",]
Plot_nmds_oc(Rhizo_rppos_final_ERLE)
Plot_nmds_oc(Rhizo_rppos_final_DICA)


## Rp negative (select and replace rppos to rpneg; RP_Pos to RP_Neg)---------------------
rpneg_raw <- read_csv("Rhizo_LCMS_RPNEG.csv") %>%
  #Make the compound IDs the row names
  column_to_rownames('Compound ID') %>%
  #Transnege the matrix to make it so samples are rows and compounds are columns
  t()  %>%
  as.data.frame()

rpneg_pareto <- rpneg_raw %>%
  log_transform() %>%
  pareto_scale()

#Combine metadata with metabolite data
Rhizo_rpneg_final  <- rpneg_pareto %>%
  rownames_to_column('SampleID_rhizo') %>%
  mutate(SampleID_rhizo = sub("_RP_Neg$", "", SampleID_rhizo)) %>% 
  right_join(metadata)%>%
  column_to_rownames('SampleID_rhizo')

#NMDS plot
Plot_nmds(Rhizo_rpneg_final)

Rhizo_rpneg_final_ERLE <- Rhizo_rpneg_final[Rhizo_rpneg_final$plant== "ERLE",]
Rhizo_rpneg_final_DICA <- Rhizo_rpneg_final[Rhizo_rpneg_final$plant== "DICA",]
Plot_nmds_oc(Rhizo_rpneg_final_ERLE)
Plot_nmds_oc(Rhizo_rpneg_final_DICA)



## HILIC positive-----------------
hpos_raw <- read_csv("Rhizo_LCMS_HPOS.csv") %>%
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
Rhizo_hpos_final  <- hpos_pareto %>%
  rownames_to_column('SampleID_rhizo') %>%
  mutate(SampleID_rhizo = sub("_HILIC_Po$", "", SampleID_rhizo)) %>% 
  right_join(metadata)%>%
  column_to_rownames('SampleID_rhizo')

#NMDS plot
Plot_nmds(Rhizo_hpos_final)
Rhizo_hpos_final_ERLE <- Rhizo_hpos_final[Rhizo_hpos_final$plant== "ERLE",]
Rhizo_hpos_final_DICA <- Rhizo_hpos_final[Rhizo_hpos_final$plant== "DICA",]
Plot_nmds_oc(Rhizo_hpos_final_ERLE)
Plot_nmds_oc(Rhizo_hpos_final_DICA)

# HILIC negative-----------------
hneg_raw <- read_csv("Rhizo_LCMS_HNEG.csv") %>%
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
Rhizo_hneg_final  <- hneg_pareto %>%
  rownames_to_column('SampleID_rhizo') %>%
  mutate(SampleID_rhizo = sub("_HILIC_Neg$", "", SampleID_rhizo)) %>% 
  left_join(metadata)%>%
  column_to_rownames('SampleID_rhizo')

Plot_nmds(Rhizo_hneg_final)
Rhizo_hneg_final_ERLE <- Rhizo_hneg_final[Rhizo_hneg_final$plant== "ERLE",]
Rhizo_hneg_final_DICA <- Rhizo_hneg_final[Rhizo_hneg_final$plant== "DICA",]
Plot_nmds_oc(Rhizo_hneg_final_ERLE)
Plot_nmds_oc(Rhizo_hneg_final_DICA)


# differential analysis ------------
##rppos########-------
p.plant <- apply(Rhizo_rppos_final[,1:ncol(rppos_raw)], 2, function(x){kruskal.test(x~plant, data = Rhizo_rppos_final)$p.value})  
Rhizo_rppos_LC2  <- rppos_raw %>%
  rownames_to_column('SampleID_rhizo') %>%
  mutate(SampleID_rhizo = sub("_RP_Pos$", "", SampleID_rhizo)) %>% 
  right_join(metadata)

DICA <- Rhizo_rppos_LC2[which(Rhizo_rppos_LC2$plant == "DICA"),1:ncol(rppos_raw)+1]
ELRE <- Rhizo_rppos_LC2[which(Rhizo_rppos_LC2$plant == "ERLE"),1:ncol(rppos_raw)+1]
LC2_rppos <- log2(colMeans(DICA)/colMeans(ELRE))
volcano_rppos <- data.frame(pvalue = p.plant,
                            Log2FC = LC2_rppos)  

Plot_Rhizo_volcano_rppos<-Plot_volcano(volcano_rppos)

##rpneg-----------
p.plant <- apply(Rhizo_rpneg_final[,1:ncol(rpneg_raw)], 2, function(x){kruskal.test(x~plant, data = Rhizo_rpneg_final)$p.value})  

Rhizo_rpneg_LC2  <- rpneg_raw %>%
  rownames_to_column('SampleID_rhizo') %>%
  mutate(SampleID_rhizo = sub("_RP_Neg$", "", SampleID_rhizo)) %>% 
  right_join(metadata)

DICA <- Rhizo_rpneg_LC2[which(Rhizo_rpneg_LC2$plant == "DICA"),1:ncol(rpneg_raw)+1]
ELRE <- Rhizo_rpneg_LC2[which(Rhizo_rpneg_LC2$plant == "ERLE"),1:ncol(rpneg_raw)+1]
LC2_rpneg <- log2(colMeans(DICA)/colMeans(ELRE))
volcano_rpneg <- data.frame(pvalue = p.plant,
                            Log2FC = LC2_rpneg)  

Plot_Rhizo_volcano_rpneg<-Plot_volcano(volcano_rpneg) 



##hpos-----------
p.plant <- apply(Rhizo_hpos_final[,1:ncol(hpos_raw)], 2, function(x){kruskal.test(x~plant, data = Rhizo_hpos_final)$p.value})  

Rhizo_hpos_LC2  <- hpos_raw %>%
  rownames_to_column('SampleID_rhizo') %>%
  mutate(SampleID_rhizo = sub("_HILIC_Po$", "", SampleID_rhizo)) %>% 
  right_join(metadata)

DICA <- Rhizo_hpos_LC2[which(Rhizo_hpos_LC2$plant == "DICA"),1:ncol(hpos_raw)+1]
ELRE <- Rhizo_hpos_LC2[which(Rhizo_hpos_LC2$plant == "ERLE"),1:ncol(hpos_raw)+1]
LC2_hpos <- log2(colMeans(DICA)/colMeans(ELRE))
volcano_hpos <- data.frame(pvalue = p.plant,
                            Log2FC = LC2_hpos)  

Plot_Rhizo_volcano_hpos <- Plot_volcano(volcano_hpos)


##hneg-----------
p.plant <- apply(Rhizo_hneg_final[,1:ncol(hneg_raw)], 2, function(x){kruskal.test(x~plant, data = Rhizo_hneg_final)$p.value})  

Rhizo_hneg_LC2  <- hneg_raw %>%
  rownames_to_column('SampleID_rhizo') %>%
  mutate(SampleID_rhizo = sub("_HILIC_Neg$", "", SampleID_rhizo)) %>% 
  right_join(metadata)

DICA <- Rhizo_hneg_LC2[which(Rhizo_hneg_LC2$plant == "DICA"),1:ncol(hneg_raw)+1]
ELRE <- Rhizo_hneg_LC2[which(Rhizo_hneg_LC2$plant == "ERLE"),1:ncol(hneg_raw)+1]
LC2_hneg <- log2(colMeans(DICA)/colMeans(ELRE))
volcano_hneg <- data.frame(pvalue = p.plant,
                           Log2FC = LC2_hneg)  
Plot_Rhizo_volcano_hneg<-Plot_volcano(volcano_hneg)


########### output ----------

Diff_Rhizo_plant <- rbind(volcano_rppos, volcano_hpos, volcano_rpneg, volcano_rpneg)
# Remove duplicates in "Name" column, prioritizing rows with non-NA values in "delabel" column


write.table(Diff_Rhizo_plant, file = "Diff_Rhizo_plant.csv", sep = ",", row.names = FALSE)

# ERLE open vs close and DICA open vs close -----
##rppos-------
e <- Rhizo_rppos_final[Rhizo_rppos_final$plant == "ERLE", ]
d <- Rhizo_rppos_final[Rhizo_rppos_final$plant == "DICA", ]

p_eoc <- apply(e[,c(1:ncol(rppos_raw))],2, function(x){kruskal.test(x~canopy, data = e)$p.value})
p_doc <- apply(d[,c(1:ncol(rppos_raw))],2, function(x){kruskal.test(x~canopy, data = d)$p.value})

LC2.eo <- Rhizo_rppos_LC2[Rhizo_rppos_LC2$plant == "ERLE" & Rhizo_rppos_LC2$canopy == "open", 1:ncol(rppos_raw)+1]
LC2.ec <- Rhizo_rppos_LC2[Rhizo_rppos_LC2$plant == "ERLE" & Rhizo_rppos_LC2$canopy == "close", 1:ncol(rppos_raw)+1]
LC2.do <- Rhizo_rppos_LC2[Rhizo_rppos_LC2$plant == "DICA" & Rhizo_rppos_LC2$canopy == "open", 1:ncol(rppos_raw)+1]
LC2.dc <- Rhizo_rppos_LC2[Rhizo_rppos_LC2$plant == "DICA" & Rhizo_rppos_LC2$canopy == "close", 1:ncol(rppos_raw)+1]

LC2.eoc <- log2(colMeans(LC2.ec)/colMeans(LC2.eo))
LC2.doc <- log2(colMeans(LC2.dc)/colMeans(LC2.do))

volcano_rppos_eoc <- data.frame(pvalue = p_eoc,
                                Log2FC = LC2.eoc) 
volcano_rppos_doc <- data.frame(pvalue = p_doc,
                                Log2FC = LC2.doc) 

Plot_Rhizo_volcano_rppos_eoc <- Plot_volcano_oc(volcano_rppos_eoc)
Plot_Rhizo_volcano_rppos_doc <- Plot_volcano_oc(volcano_rppos_doc)
ggsave("Plot_Rhizo_volcano_rppos_eoc.pdf", Plot_Rhizo_volcano_rppos_eoc + ggtitle("rppos ERLE by canopy") + theme(legend.position = "none"), width = 6, height = 4)
ggsave("Plot_Rhizo_volcano_rppos_doc.pdf", Plot_Rhizo_volcano_rppos_doc + ggtitle("rppos DICA by canopy") + theme(legend.position = "none"), width = 6, height = 4)


##rpneg-------
e <- Rhizo_rpneg_final[Rhizo_rpneg_final$plant == "ERLE", ]
d <- Rhizo_rpneg_final[Rhizo_rpneg_final$plant == "DICA", ]

p_eoc <- apply(e[,c(1:ncol(rpneg_raw))],2, function(x){kruskal.test(x~canopy, data = e)$p.value})
p_doc <- apply(d[,c(1:ncol(rpneg_raw))],2, function(x){kruskal.test(x~canopy, data = d)$p.value})


LC2.eo <- Rhizo_rpneg_LC2[Rhizo_rpneg_LC2$plant == "ERLE" & Rhizo_rpneg_LC2$canopy == "open", 1:ncol(rpneg_raw)+1]
LC2.ec <- Rhizo_rpneg_LC2[Rhizo_rpneg_LC2$plant == "ERLE" & Rhizo_rpneg_LC2$canopy == "close", 1:ncol(rpneg_raw)+1]
LC2.do <- Rhizo_rpneg_LC2[Rhizo_rpneg_LC2$plant == "DICA" & Rhizo_rpneg_LC2$canopy == "open", 1:ncol(rpneg_raw)+1]
LC2.dc <- Rhizo_rpneg_LC2[Rhizo_rpneg_LC2$plant == "DICA" & Rhizo_rpneg_LC2$canopy == "close", 1:ncol(rpneg_raw)+1]

LC2.eoc <- log2(colMeans(LC2.ec)/colMeans(LC2.eo))
LC2.doc <- log2(colMeans(LC2.dc)/colMeans(LC2.do))

volcano_rpneg_eoc <- data.frame(pvalue = p_eoc,
                                Log2FC = LC2.eoc) 
volcano_rpneg_doc <- data.frame(pvalue = p_doc,
                                Log2FC = LC2.doc) 

Plot_Rhizo_volcano_rpneg_eoc <- Plot_volcano_oc(volcano_rpneg_eoc)
Plot_Rhizo_volcano_rpneg_doc <- Plot_volcano_oc(volcano_rpneg_doc)
ggsave("Plot_Rhizo_volcano_rpneg_eoc.pdf", Plot_Rhizo_volcano_rpneg_eoc + ggtitle("rpneg ERLE by canopy") + theme(legend.position = "none"), width = 6, height = 4)
ggsave("Plot_Rhizo_volcano_rpneg_doc.pdf", Plot_Rhizo_volcano_rpneg_doc + ggtitle("rpneg DICA by canopy") + theme(legend.position = "none"), width = 6, height = 4)


##hpos-------
e <- Rhizo_hpos_final[Rhizo_hpos_final$plant == "ERLE", ]
d <- Rhizo_hpos_final[Rhizo_hpos_final$plant == "DICA", ]

p_eoc <- apply(e[,c(1:ncol(hpos_raw))],2, function(x){kruskal.test(x~canopy, data = e)$p.value})
p_doc <- apply(d[,c(1:ncol(hpos_raw))],2, function(x){kruskal.test(x~canopy, data = d)$p.value})

LC2.eo <- Rhizo_hpos_LC2[Rhizo_hpos_LC2$plant == "ERLE" & Rhizo_hpos_LC2$canopy == "open", 1:ncol(hpos_raw)+1]
LC2.ec <- Rhizo_hpos_LC2[Rhizo_hpos_LC2$plant == "ERLE" & Rhizo_hpos_LC2$canopy == "close", 1:ncol(hpos_raw)+1]
LC2.do <- Rhizo_hpos_LC2[Rhizo_hpos_LC2$plant == "DICA" & Rhizo_hpos_LC2$canopy == "open", 1:ncol(hpos_raw)+1]
LC2.dc <- Rhizo_hpos_LC2[Rhizo_hpos_LC2$plant == "DICA" & Rhizo_hpos_LC2$canopy == "close", 1:ncol(hpos_raw)+1]

LC2.eoc <- log2(colMeans(LC2.ec)/colMeans(LC2.eo))
LC2.doc <- log2(colMeans(LC2.dc)/colMeans(LC2.do))

volcano_hpos_eoc <- data.frame(pvalue = p_eoc,
                               Log2FC = LC2.eoc) 
volcano_hpos_doc <- data.frame(pvalue = p_doc,
                               Log2FC = LC2.doc) 

Plot_Rhizo_volcano_hpos_eoc <- Plot_volcano_oc(volcano_hpos_eoc)
Plot_Rhizo_volcano_hpos_doc <- Plot_volcano_oc(volcano_hpos_doc)
ggsave("Plot_Rhizo_volcano_hpos_eoc.pdf", Plot_Rhizo_volcano_hpos_eoc + ggtitle("hpos ERLE by canopy") + theme(legend.position = "none"), width = 6, height = 4)
ggsave("Plot_Rhizo_volcano_hpos_doc.pdf", Plot_Rhizo_volcano_hpos_doc + ggtitle("hpos DICA by canopy") + theme(legend.position = "none"), width = 6, height = 4)



##hneg-------
e <- Rhizo_hneg_final[Rhizo_hneg_final$plant == "ERLE", ]
d <- Rhizo_hneg_final[Rhizo_hneg_final$plant == "DICA", ]

p_eoc <- apply(e[,c(1:ncol(hneg_raw))],2, function(x){kruskal.test(x~canopy, data = e)$p.value})
p_doc <- apply(d[,c(1:ncol(hneg_raw))],2, function(x){kruskal.test(x~canopy, data = d)$p.value})


LC2.eo <- Rhizo_hneg_LC2[Rhizo_hneg_LC2$plant == "ERLE" & Rhizo_hneg_LC2$canopy == "open", 1:ncol(hneg_raw)+1]
LC2.ec <- Rhizo_hneg_LC2[Rhizo_hneg_LC2$plant == "ERLE" & Rhizo_hneg_LC2$canopy == "close", 1:ncol(hneg_raw)+1]
LC2.do <- Rhizo_hneg_LC2[Rhizo_hneg_LC2$plant == "DICA" & Rhizo_hneg_LC2$canopy == "open", 1:ncol(hneg_raw)+1]
LC2.dc <- Rhizo_hneg_LC2[Rhizo_hneg_LC2$plant == "DICA" & Rhizo_hneg_LC2$canopy == "close", 1:ncol(hneg_raw)+1]

LC2.eoc <- log2(colMeans(LC2.ec)/colMeans(LC2.eo))
LC2.doc <- log2(colMeans(LC2.dc)/colMeans(LC2.do))

volcano_hneg_eoc <- data.frame(pvalue = p_eoc,
                               Log2FC = LC2.eoc) 
volcano_hneg_doc <- data.frame(pvalue = p_doc,
                               Log2FC = LC2.doc) 

Plot_Rhizo_volcano_hneg_eoc <- Plot_volcano_oc(volcano_hneg_eoc)
Plot_Rhizo_volcano_hneg_doc <- Plot_volcano_oc(volcano_hneg_doc)
ggsave("Plot_Rhizo_volcano_hneg_eoc.pdf", Plot_Rhizo_volcano_hneg_eoc + ggtitle("hneg ERLE by canopy") + theme(legend.position = "none"), width = 6, height = 4)
ggsave("Plot_Rhizo_volcano_hneg_doc.pdf", Plot_Rhizo_volcano_hneg_doc + ggtitle("hneg DICA by canopy") + theme(legend.position = "none"), width = 6, height = 4)


# output -------

Diff_Rhizo_ERLE_canopy <- rbind(volcano_rppos_eoc, volcano_hpos_eoc, volcano_rpneg_eoc, volcano_hneg_eoc)
# Remove duplicates in "Name" column, prioritizing rows with non-NA values in "delabel" column
Diff_Rhizo_DICA_canopy <- rbind(volcano_rppos_doc, volcano_hpos_doc, volcano_rpneg_doc, volcano_hneg_doc)


write.table(Diff_Rhizo_ERLE_canopy, file = "Diff_Rhizo_ERLE_canopy.csv", sep = ",", row.names = FALSE)
write.table(Diff_Rhizo_DICA_canopy, file = "Diff_Rhizo_DICA_canopy.csv", sep = ",", row.names = FALSE)



# now compare open ERLE vs DICA and close ERLE vs DICA -----
##rppos-------
o <- Rhizo_rppos_final[Rhizo_rppos_final$canopy == "open", ]
c <- Rhizo_rppos_final[Rhizo_rppos_final$canopy == "close", ]

p_oed <- apply(o[,c(1:ncol(rppos_raw))],2, function(x){kruskal.test(x~plant, data = o)$p.value})
p_ced <- apply(c[,c(1:ncol(rppos_raw))],2, function(x){kruskal.test(x~plant, data = c)$p.value})

LC2.eo <- Rhizo_rppos_LC2[Rhizo_rppos_LC2$plant == "ERLE" & Rhizo_rppos_LC2$canopy == "open", 1:ncol(rppos_raw)+1]
LC2.ec <- Rhizo_rppos_LC2[Rhizo_rppos_LC2$plant == "ERLE" & Rhizo_rppos_LC2$canopy == "close", 1:ncol(rppos_raw)+1]
LC2.do <- Rhizo_rppos_LC2[Rhizo_rppos_LC2$plant == "DICA" & Rhizo_rppos_LC2$canopy == "open", 1:ncol(rppos_raw)+1]
LC2.dc <- Rhizo_rppos_LC2[Rhizo_rppos_LC2$plant == "DICA" & Rhizo_rppos_LC2$canopy == "close", 1:ncol(rppos_raw)+1]

LC2.oed <- log2(colMeans(LC2.do)/colMeans(LC2.eo))
LC2.ced <- log2(colMeans(LC2.dc)/colMeans(LC2.ec))

volcano_rppos_oed <- data.frame(pvalue = p_oed,
                                Log2FC = LC2.oed) 
volcano_rppos_ced <- data.frame(pvalue = p_ced,
                                Log2FC = LC2.ced) 

Plot_Rhizo_volcano_rppos_oed <- Plot_volcano(volcano_rppos_oed)
Plot_Rhizo_volcano_rppos_ced <- Plot_volcano(volcano_rppos_ced)
ggsave("Plot_Rhizo_volcano_rppos_oed.pdf", Plot_Rhizo_volcano_rppos_oed + ggtitle("rppos open ERLE vs DICA") + theme(legend.position = "none"), width = 6, height = 4)
ggsave("Plot_Rhizo_volcano_rppos_ced.pdf", Plot_Rhizo_volcano_rppos_ced + ggtitle("rppos close ERLE vs DICA") + theme(legend.position = "none"), width = 6, height = 4)


##rpneg-------
o <- Rhizo_rpneg_final[Rhizo_rpneg_final$canopy == "open", ]
c <- Rhizo_rpneg_final[Rhizo_rpneg_final$canopy == "close", ]

p_oed <- apply(o[,c(1:ncol(rpneg_raw))],2, function(x){kruskal.test(x~plant, data = o)$p.value})
p_ced <- apply(c[,c(1:ncol(rpneg_raw))],2, function(x){kruskal.test(x~plant, data = c)$p.value})

LC2.eo <- Rhizo_rpneg_LC2[Rhizo_rpneg_LC2$plant == "ERLE" & Rhizo_rpneg_LC2$canopy == "open", 1:ncol(rpneg_raw)+1]
LC2.ec <- Rhizo_rpneg_LC2[Rhizo_rpneg_LC2$plant == "ERLE" & Rhizo_rpneg_LC2$canopy == "close", 1:ncol(rpneg_raw)+1]
LC2.do <- Rhizo_rpneg_LC2[Rhizo_rpneg_LC2$plant == "DICA" & Rhizo_rpneg_LC2$canopy == "open", 1:ncol(rpneg_raw)+1]
LC2.dc <- Rhizo_rpneg_LC2[Rhizo_rpneg_LC2$plant == "DICA" & Rhizo_rpneg_LC2$canopy == "close", 1:ncol(rpneg_raw)+1]

LC2.oed <- log2(colMeans(LC2.do)/colMeans(LC2.eo))
LC2.ced <- log2(colMeans(LC2.dc)/colMeans(LC2.ec))

volcano_rpneg_oed <- data.frame(pvalue = p_oed,
                                Log2FC = LC2.oed) 
volcano_rpneg_ced <- data.frame(pvalue = p_ced,
                                Log2FC = LC2.ced) 

Plot_Rhizo_volcano_rpneg_oed <- Plot_volcano(volcano_rpneg_oed)
Plot_Rhizo_volcano_rpneg_ced <- Plot_volcano(volcano_rpneg_ced)
ggsave("Plot_Rhizo_volcano_rpneg_oed.pdf", Plot_Rhizo_volcano_rpneg_oed + ggtitle("rpneg open ERLE vs DICA") + theme(legend.position = "none"), width = 6, height = 4)
ggsave("Plot_Rhizo_volcano_rpneg_ced.pdf", Plot_Rhizo_volcano_rpneg_ced + ggtitle("rpneg close ERLE vs DICA") + theme(legend.position = "none"), width = 6, height = 4)

##hpos-------
o <- Rhizo_hpos_final[Rhizo_hpos_final$canopy == "open", ]
c <- Rhizo_hpos_final[Rhizo_hpos_final$canopy == "close", ]

p_oed <- apply(o[,c(1:ncol(hpos_raw))],2, function(x){kruskal.test(x~plant, data = o)$p.value})
p_ced <- apply(c[,c(1:ncol(hpos_raw))],2, function(x){kruskal.test(x~plant, data = c)$p.value})

LC2.eo <- Rhizo_hpos_LC2[Rhizo_hpos_LC2$plant == "ERLE" & Rhizo_hpos_LC2$canopy == "open", 1:ncol(hpos_raw)+1]
LC2.ec <- Rhizo_hpos_LC2[Rhizo_hpos_LC2$plant == "ERLE" & Rhizo_hpos_LC2$canopy == "close", 1:ncol(hpos_raw)+1]
LC2.do <- Rhizo_hpos_LC2[Rhizo_hpos_LC2$plant == "DICA" & Rhizo_hpos_LC2$canopy == "open", 1:ncol(hpos_raw)+1]
LC2.dc <- Rhizo_hpos_LC2[Rhizo_hpos_LC2$plant == "DICA" & Rhizo_hpos_LC2$canopy == "close", 1:ncol(hpos_raw)+1]

LC2.oed <- log2(colMeans(LC2.do)/colMeans(LC2.eo))
LC2.ced <- log2(colMeans(LC2.dc)/colMeans(LC2.ec))

volcano_hpos_oed <- data.frame(pvalue = p_oed,
                               Log2FC = LC2.oed) 
volcano_hpos_ced <- data.frame(pvalue = p_ced,
                               Log2FC = LC2.ced) 

Plot_Rhizo_volcano_hpos_oed <- Plot_volcano(volcano_hpos_oed)
Plot_Rhizo_volcano_hpos_ced <- Plot_volcano(volcano_hpos_ced)
ggsave("Plot_Rhizo_volcano_hpos_oed.pdf", Plot_Rhizo_volcano_hpos_oed + ggtitle("hpos open ERLE vs DICA") + theme(legend.position = "none"), width = 6, height = 4)
ggsave("Plot_Rhizo_volcano_hpos_ced.pdf", Plot_Rhizo_volcano_hpos_ced + ggtitle("hpos close ERLE vs DICA") + theme(legend.position = "none"), width = 6, height = 4)


##hneg-------
o <- Rhizo_hneg_final[Rhizo_hneg_final$canopy == "open", ]
c <- Rhizo_hneg_final[Rhizo_hneg_final$canopy == "close", ]

p_oed <- apply(o[,c(1:ncol(hneg_raw))],2, function(x){kruskal.test(x~plant, data = o)$p.value})
p_ced <- apply(c[,c(1:ncol(hneg_raw))],2, function(x){kruskal.test(x~plant, data = c)$p.value})

LC2.eo <- Rhizo_hneg_LC2[Rhizo_hneg_LC2$plant == "ERLE" & Rhizo_hneg_LC2$canopy == "open", 1:ncol(hneg_raw)+1]
LC2.ec <- Rhizo_hneg_LC2[Rhizo_hneg_LC2$plant == "ERLE" & Rhizo_hneg_LC2$canopy == "close", 1:ncol(hneg_raw)+1]
LC2.do <- Rhizo_hneg_LC2[Rhizo_hneg_LC2$plant == "DICA" & Rhizo_hneg_LC2$canopy == "open", 1:ncol(hneg_raw)+1]
LC2.dc <- Rhizo_hneg_LC2[Rhizo_hneg_LC2$plant == "DICA" & Rhizo_hneg_LC2$canopy == "close", 1:ncol(hneg_raw)+1]

LC2.oed <- log2(colMeans(LC2.do)/colMeans(LC2.eo))
LC2.ced <- log2(colMeans(LC2.dc)/colMeans(LC2.ec))

volcano_hneg_oed <- data.frame(pvalue = p_oed,
                               Log2FC = LC2.oed) 
volcano_hneg_ced <- data.frame(pvalue = p_ced,
                               Log2FC = LC2.ced) 

Plot_Rhizo_volcano_hneg_oed <- Plot_volcano(volcano_hneg_oed)
Plot_Rhizo_volcano_hneg_ced <- Plot_volcano(volcano_hneg_ced)
ggsave("Plot_Rhizo_volcano_hneg_oed.pdf", Plot_Rhizo_volcano_hneg_oed + ggtitle("hneg open ERLE vs DICA") + theme(legend.position = "none"), width = 6, height = 4)
ggsave("Plot_Rhizo_volcano_hneg_ced.pdf", Plot_Rhizo_volcano_hneg_ced + ggtitle("hneg close ERLE vs DICA") + theme(legend.position = "none"), width = 6, height = 4)



# make table -------
Diff_Rhizo_open_plant <- rbind(volcano_rppos_oed, volcano_hpos_oed, volcano_rpneg_oed, volcano_hneg_oed)
# Remove duplicates in "Name" column, prioritizing rows with non-NA values in "delabel" column

Diff_Rhizo_close_plant <- rbind(volcano_rppos_ced, volcano_hpos_ced, volcano_rpneg_ced, volcano_hneg_ced)


write.table(Diff_Rhizo_open_plant, file = "Diff_Rhizo_open_plant.csv", sep = ",", row.names = FALSE)
write.table(Diff_Rhizo_close_plant, file = "Diff_Rhizo_close_plant.csv", sep = ",", row.names = FALSE)





