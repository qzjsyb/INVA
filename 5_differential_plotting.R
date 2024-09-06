# Environment -------------------------------------------------------------
library(tidyverse)
library(mixOmics)
library(gridExtra)
library(rstatix)
library(agricolae)
setwd('~/INVA/INVA')
metadata <- read.table("metadata.txt", header = TRUE, sep = "\t", row.names = 1)
metadata$block <- as.factor(metadata$block)
Root_compounds <- read.table("Root_compound_link.txt", header = TRUE, sep = "\t")
Rhizo_compounds <- read.table("Rhizo_compound_link.txt", header = TRUE, sep = "\t")
metadata$canopy <- factor(metadata$canopy, level = c("open", "close"))


plot_compound_rhizo <- function(compound_id) {
  compound_id <- "FT_205.09723_2.132_rppos"
  # Extract the last part of the compound ID to identify the dataframe
  df_suffix <- sub("^.*_", "", compound_id)
  df_name <- paste0("Rhizo_", tolower(df_suffix), "_final")
  
  # Get the dataframe corresponding to the compound ID
  df <- get(df_name)
  values <- df[, compound_id]
  # Extract the compound name from the dataframe plantcanopy
  compound_name <- Rhizo_compounds$Name[Rhizo_compounds$Compound.ID == compound_id]
  
  #calculate 
  canopy_formula <- reformulate("canopy", response = compound_id)
  canopy <- kruskal.test(canopy_formula, data = df)
  plant_formula <- reformulate("plant", response = compound_id)
  plant <- kruskal.test(plant_formula, data = df)
  
  df$group <- paste(df$plant, df$canopy, sep = "_")
  y <- as.data.frame(values, df$group)%>%
    rownames_to_column("group")
  
  
  replace_outliers <- function(x, ...) {
    qnt <- quantile(x, probs = c(.25, .75), ...)
    val <- 1.5 * IQR(x)
    y <- x
    y[x < (qnt[1] - val)] <- NA
    y[x > (qnt[2] + val)] <- NA
    y
  }
  
  dc <- quantile(replace_outliers(y[y$group == "DICA_close",]$values), probs = 1,na.rm = TRUE)
  do <- quantile(replace_outliers(y[y$group == "DICA_open",]$values), probs = 1,na.rm = TRUE)
  ec <- quantile(replace_outliers(y[y$group == "ERLE_close",]$values), probs = 1,na.rm = TRUE)
  eo <- quantile(replace_outliers(y[y$group == "ERLE_open",]$values), probs = 1,na.rm = TRUE)
  
  #set compact letter display
  cld<- kruskal(values, df$group)$groups
  cld <- cld[c("DICA_close","DICA_open","ERLE_close","ERLE_open" ),]
  
  cld$y <- c(dc, do, ec, eo)
  
  # Plotting
  gg<- df %>%
    ggplot(aes(x = plant, y = !!sym(compound_id))) +
    geom_boxplot(aes(fill = plant, shape = canopy), outlier.shape = NA) +
    scale_fill_manual(values = c("ERLE" = "#9765A8", "DICA" = "#4DAF4A")) +
    scale_color_manual(values = c("ERLE" = "#9765A8", "DICA" = "#4DAF4A")) +
    labs(title = paste0(compound_name)) +
    xlab(NULL) +
    ylab("LCMS_intensity") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    geom_text(data = cld, aes(label = groups, x = 1, y = y), vjust = -0.5)+
    theme(legend.position = "none") + 
    facet_wrap(~canopy)+  
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  # Add cld annotations on top of each box
  
  # Add p-values to the bottom
  gg <- gg + 
    labs(caption = paste("P(plant) =", round(plant$p.value, 3), "; P(canopy) =", round(canopy$p.value, 3)))
  
  return(gg)
}

plot_compound_root <- function(compound_id) {
  
  compound_id <- "FT_188.06982_5.454_HPOS"
  # Extract the last part of the compound ID to identify the dataframe
  df_suffix <- sub("^.*_", "", compound_id)
  df_name <- paste0("Root_", tolower(df_suffix), "_final")
  
  # Get the dataframe corresponding to the compound ID
  df <- get(df_name)
  values <- df[, compound_id]
  # Extract the compound name from the dataframe plantcanopy
  compound_name <- Root_compounds$Name[Root_compounds$Compound.ID == compound_id]
  
  #calculate 
  canopy_formula <- reformulate("canopy", response = compound_id)
  canopy <- kruskal.test(canopy_formula, data = df)
  plant_formula <- reformulate("plant", response = compound_id)
  plant <- kruskal.test(plant_formula, data = df)
  
  df$group <- paste(df$plant, df$canopy, sep = "_")
  y <- as.data.frame(values, df$group)%>%
    rownames_to_column("group")
  
  
  replace_outliers <- function(x, ...) {
    qnt <- quantile(x, probs = c(.25, .75), ...)
    val <- 1.5 * IQR(x)
    y <- x
    y[x < (qnt[1] - val)] <- NA
    y[x > (qnt[2] + val)] <- NA
    y
  }
  
  dc <- quantile(replace_outliers(y[y$group == "DICA_close",]$values), probs = 1,na.rm = TRUE)
  do <- quantile(replace_outliers(y[y$group == "DICA_open",]$values), probs = 1,na.rm = TRUE)
  ec <- quantile(replace_outliers(y[y$group == "ERLE_close",]$values), probs = 1,na.rm = TRUE)
  eo <- quantile(replace_outliers(y[y$group == "ERLE_open",]$values), probs = 1,na.rm = TRUE)
  
  #set compact letter display
  cld<- kruskal(values, df$group)$groups
  cld <- cld[c("DICA_close","DICA_open","ERLE_close","ERLE_open" ),]
  
  cld$y <- c(dc, do, ec, eo)

  
  
  # Plotting
  gg<- df %>%
    ggplot(aes(x = plant, y = !!sym(compound_id))) +
    geom_boxplot(aes(fill = plant), outlier.shape = NA) +
    scale_fill_manual(values = c("ERLE" = "#9765A8", "DICA" = "#4DAF4A")) +
    scale_color_manual(values = c("ERLE" = "#9765A8", "DICA" = "#4DAF4A")) +
    labs(title = paste0(compound_name)) +
    xlab(NULL) +
    ylab("LCMS_intensity") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    geom_text(data = cld, aes(label = groups, x = 1, y = y), vjust = -0.5)+
    theme(legend.position = "none") + 
    facet_wrap(~canopy)+  
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  # Add cld annotations on top of each box

  # Add p-values to the bottom
  gg <- gg + 
    labs(caption = paste("P(plant) =", round(plant$p.value, 3), "; P(canopy) =", round(canopy$p.value, 3)))
  
  return(gg)
}

dunn <- function(compound_id) {
  # Extract the last part of the compound ID to identify the dataframe
  df_suffix <- sub("^.*_", "", compound_id)
  df_name <- paste0("Root_", tolower(df_suffix), "_final")
  
  # Get the dataframe corresponding to the compound ID
  df <- get(df_name)
  values <- df[, compound_id]
  # Extract the compound name from the dataframe plantcanopy
  compound_name <- Root_compounds$Name[Root_compounds$Compound.ID == compound_id]
  
  #calculate 
  canopy_formula <- reformulate("canopy", response = compound_id)
  canopy <- kruskal.test(canopy_formula, data = df)
  plant_formula <- reformulate("plant", response = compound_id)
  plant <- kruskal.test(plant_formula, data = df)
  
  df$group <- paste(df$plant, df$canopy, sep = "_")
  y <- as.data.frame(values, df$group)%>%
    rownames_to_column("group")
return(dunn_test(values~group, data = y))
  }


#ggsave("malditest.pdf",plot_compound_maldi("FT_303.08647_5.798_RPPOS",Root_rppos_final, c("R06","R21")), width = 3, height =4)

## load data ---------------------
Root_rppos_raw <- read_csv("Root_LCMS_RPPOS.csv") %>%
  #Make the compound IDs the row names
  column_to_rownames('Compound ID') %>%
  #Transpose the matrix to make it so samples are rows and compounds are columns
  t()  %>%
  as.data.frame()

Root_rpneg_raw <- read_csv("Root_LCMS_RPNEG.csv") %>%
  column_to_rownames('Compound ID') %>%
  t()  %>%
  as.data.frame()

Root_hpos_raw <- read_csv("Root_LCMS_HPOS.csv") %>%
  column_to_rownames('Compound ID') %>%
  t()  %>%
  as.data.frame()

Root_hneg_raw <- read_csv("Root_LCMS_HNEG.csv") %>%
  column_to_rownames('Compound ID') %>%
  t()  %>%
  as.data.frame()

Rhizo_rppos_raw <- read_csv("Rhizo_LCMS_RPPOS.csv") %>%
  column_to_rownames('Compound ID') %>%
  t()  %>%
  as.data.frame()

Rhizo_rpneg_raw <- read_csv("Rhizo_LCMS_RPNEG.csv") %>%
  column_to_rownames('Compound ID') %>%
  t()  %>%
  as.data.frame()

Rhizo_hpos_raw <- read_csv("Rhizo_LCMS_HPOS.csv") %>%
  column_to_rownames('Compound ID') %>%
  t()  %>%
  as.data.frame()

Rhizo_hneg_raw <- read_csv("Rhizo_LCMS_HNEG.csv") %>%
  column_to_rownames('Compound ID') %>%
  t()  %>%
  as.data.frame()

Root_rppos_final  <- Root_rppos_raw %>%
  rownames_to_column('SampleID_root') %>%
  mutate(SampleID_root = sub("_RP_Pos$", "", SampleID_root)) %>% 
  right_join(metadata)

Root_rpneg_final  <- Root_rpneg_raw %>%
  rownames_to_column('SampleID_root') %>%
  mutate(SampleID_root = sub("_RP_Neg$", "", SampleID_root)) %>% 
  right_join(metadata)

Root_hpos_final  <- Root_hpos_raw %>%
  rownames_to_column('SampleID_root') %>%
  mutate(SampleID_root = sub("_HILIC_Pos$", "", SampleID_root)) %>% 
  right_join(metadata)

Root_hneg_final  <- Root_hneg_raw %>%
  rownames_to_column('SampleID_root') %>%
  mutate(SampleID_root = sub("_HILIC_Neg$", "", SampleID_root)) %>% 
  right_join(metadata)

Rhizo_rppos_final  <- Rhizo_rppos_raw %>%
  rownames_to_column('SampleID_rhizo') %>%
  mutate(SampleID_rhizo = sub("_RP_Pos$", "", SampleID_rhizo)) %>% 
  right_join(metadata)

Rhizo_rpneg_final  <- Rhizo_rpneg_raw %>%
  rownames_to_column('SampleID_rhizo') %>%
  mutate(SampleID_rhizo = sub("_RP_Neg$", "", SampleID_rhizo)) %>% 
  right_join(metadata)

Rhizo_hpos_final  <- Rhizo_hpos_raw %>%
  rownames_to_column('SampleID_rhizo') %>%
  mutate(SampleID_rhizo = sub("_HILIC_Po$", "", SampleID_rhizo)) %>% 
  right_join(metadata)

Rhizo_hneg_final  <- Rhizo_hneg_raw %>%
  rownames_to_column('SampleID_rhizo') %>%
  mutate(SampleID_rhizo = sub("_HILIC_Neg$", "", SampleID_rhizo)) %>% 
  right_join(metadata)

### plotting -----



plot_compound_root("FT_157.03599_4.756_HNEG") #allantoin higher in close, particularly with DICA
plot_compound_root("FT_157.03596_1.167_RPNEG") #allantoin higher in close


ggsave("allantoin_rpneg_root.pdf", plot_compound_root("FT_157.03596_1.167_RPNEG"), width =4, height = 3.167)
ggsave("asparagine_hneg.pdf",plot_compound_root("FT_131.04558_7.326_HNEG"), width =4, height = 3.167)
dunn("FT_131.04558_7.326_HNEG")


dunn("FT_145.06119_1.103_RPNEG")

plot_compound_root("FT_145.06113_7.297_HNEG") 
ggsave("glutamine_rppos.pdf", plot_compound_root("FT_145.06119_1.103_RPNEG"), width =4, height = 3.167) #good results
plot_compound_root("FT_147.07614_0.965_RPPOS") #good results
plot_compound_root("FT_132.11263_7.852_HPOS")


plot_compound_root("FT_137.04572_1.035_RPPOS") #Hypoxanthine, higher in open for both, lower in ERLE close
ggsave("Hypoxanthine_rppos_root.png", plot_compound_root("FT_137.04572_1.035_RPPOS"),  width =4, height = 3.5, dpi = 300)

plot_compound_root("FT_134.04463_0.973_RPPOS") # aspartate, high outliers in ERLE open, on average low in Erle close compare to DICA
plot_compound_root("FT_132.02959_1.122_RPNEG") #aspartate, similar trend.
plot_compound_root("FT_134.04414_7.266_HPOS") # higher in DICA, ERLE has large spread.
plot_compound_root("FT_132.0296_7.268_HNEG")

plot_compound_root("FT_131.0456_1.1_RPNEG")

plot_compound_root("FT_152.05594_5.761_HPOS")

plot_compound_root("FT_185.00817_1.186_RPPOS")

plot_compound_root("FT_136.06169_0.956_RPPOS")
ggsave("adenine_rppos_root.png", plot_compound_root("FT_136.06169_0.956_RPPOS"), width =4, height = 3.5, dpi = 300)

plot_compound_root("FT_268.1026_5.762_HPOS")
plot_compound_rhizo("FT_268.1029_5.756_hpos")
ggsave("deoxyguanosine_rppos.png", plot_compound_rhizo("FT_268.1029_5.756_hpos"), width =4, height = 3.5, dpi = 300)
ggsave("deoxyguanosine_rppos_root.png", plot_compound_root("FT_268.1026_5.762_HPOS"), width =4, height = 3.5, dpi = 300)

plot_compound_root("FT_152.05594_5.761_HPOS")
ggsave("guanine_hpos_root.png", plot_compound_root("FT_152.05594_5.761_HPOS"), width =4, height = 3.5, dpi = 300)

ggsave("guanosine_hpos_root.png", plot_compound_root("FT_284.09744_6.195_HPOS"), width =4, height = 3.5, dpi = 300)

plot_compound_rhizo("FT_157.03575_1.156_rpneg")
ggsave("allantoin_rpneg.png", plot_compound_rhizo("FT_157.03575_1.156_rpneg"), width =4, height = 3.5, dpi = 300)

plot_compound_rhizo("FT_134.04461_0.977_rppos")
plot_compound_rhizo("FT_134.04436_7.272_hpos")

plot_compound_rhizo("FT_152.05644_1.022_rppos")
ggsave("guanine_rppos.png", plot_compound_rhizo("FT_152.05644_1.022_rppos"), width =4, height = 3.5, dpi = 300)

plot_compound_rhizo("FT_152.05608_5.753_hpos")
ggsave("guanine_hpos.png", plot_compound_rhizo("FT_152.05608_5.753_hpos"), width =4, height = 3.5, dpi = 300)




plot_compound_rhizo("FT_137.04563_1.026_rppos") 
ggsave("hypoxanthine_rppos.png", plot_compound_rhizo("FT_137.04563_1.026_rppos"), width =4, height = 3.5, dpi = 300)


plot_compound_rhizo("FT_136.06151_0.964_rppos")
ggsave("adenine_rppos.png", plot_compound_rhizo("FT_136.06151_0.964_rppos"), width =4, height = 3.5, dpi = 300)

plot_compound_rhizo("FT_153.04043_1.034_rppos")
ggsave("xanthine_rppos.png", plot_compound_rhizo("FT_153.04043_1.034_rppos"), width =4, height = 3.5, dpi = 300)

plot_compound_rhizo("FT_284.09772_6.192_hpos")
ggsave("guanosine_hpos.png", plot_compound_rhizo("FT_284.09772_6.192_hpos"), width =4, height = 3.5, dpi = 300)

ggsave("12,13-DiHOME.pdf", plot_compound_rhizo("FT_315.25166_1.022_hpos"), width =3, height = 3.5)
ggsave("DiHOME 5.pdf", plot_compound_rhizo("FT_313.23765_10.814_rpneg"), width =3, height = 3.5)
ggsave("HODE 2.pdf", plot_compound_rhizo("FT_297.24219_9.542_rppos"), width =3, height = 3.5)
ggsave("HOME 1.pdf", plot_compound_rhizo("FT_297.24265_1.025_hneg"),width =3, height = 3.5)
ggsave("HOME 2.pdf", plot_compound_rhizo("FT_297.24268_1.328_hneg"), width =3, height = 3.5)
ggsave("Trihydroxyoctadecenoic acid.pdf", plot_compound_rhizo("FT_329.23206_1.03_hneg"), width =3, height = 3.5)
ggsave("12-hydroxy-9-cis-octadecenoic acid.pdf", plot_compound_rhizo("FT_297.24248_11.079_rpneg"), width =3, height = 3.5)


ggsave("Hpode2_rhizo.pdf", plot_compound_rhizo("FT_311.22191_9.66_rpneg"), width =3, height = 3.5)

#enemy release

ggsave("Gentisic acid.pdf", plot_compound_root("FT_153.01872_4.393_HNEG"), width =4, height = 3.167)
ggsave("indoleacylic acid.pdf",plot_compound_root("FT_188.06982_5.454_HPOS"), width =4, height = 3.167)
ggsave("Hesperetin.pdf", plot_compound_root("FT_303.08647_5.798_RPPOS"), width =4, height = 3.167)

ggsave("hpode2.pdf", plot_compound_root("FT_311.22236_10.256_RPNEG"), width =4, height = 3.167)

ggsave("hpode.pdf", plot_compound_root("FT_311.22234_1.018_HNEG"), width =4, height = 3.167)

ggsave("er1.pdf", plot_compound_root("FT_189.04002_4.426_HNEG"), width =3, height = 3.5)
ggsave("er2.pdf",plot_compound_root("FT_153.01872_4.393_HNEG"), width =3, height = 3.5)
ggsave("er3.pdf",plot_compound_root("FT_121.02895_1.02_HNEG"), width =3, height = 3.5)
ggsave("er4.pdf",plot_compound_root("FT_191.05542_6.256_HNEG"), width =3, height = 3.5)

ggsave("er5.pdf",plot_compound_root("FT_195.06577_0.838_HNEG"), width =3, height = 3.5)
ggsave("hpode_rhizo.pdf", plot_compound_rhizo("FT_311.22192_1.016_hneg"), width =4, height = 3.167)


ggsave("I3A_rhizo.pdf", plot_compound_rhizo("FT_146.05996_5.835_rppos"), width =4, height = 3.167)
ggsave("indoleacylic acid_rhizo.pdf", plot_compound_rhizo("FT_205.09723_2.132_rppos"), width = 4, height = 3.167)

