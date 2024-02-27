#load library
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)


### peak overlapping --------
root_comp <- read.csv("Root_peak_comp.csv")
#root_comp is a dataframe with the weights (could be m/z or neutrualized weight) of all features by each method, with method as columns

maldi_peak <- read.csv("MALDI_peaks.csv")
maldi_peak_clean <- maldi_peak[!duplicated(subset(maldi_peak, select = c(ion, mz, moleculeNames))), ]
#remove duplicates occur at different samples

MALDI_mz <- data.frame(
  Method = rep("MALDI_mz", length(maldi_peak_clean$mz)),  # Fill Method column with "MALDI_mz"
  MolecularWeight = maldi_peak_clean$mz  # Fill MoleculeWeight column with maldi_peak$mz
)

### root
root_comp_long <- tidyr::gather(root_comp, key = "Method", value = "MolecularWeight", na.rm = TRUE)

root_comp_long <- rbind(root_comp_long, MALDI_mz)

root_comp_long <- root_comp_long %>%
  mutate(Method = case_when(
    Method %in% c("RPPOS_mz", "RPNEG_mz", "HPOS_mz", "HNEG_mz") ~ "LCMS_mz",
    TRUE ~ Method
  ))
#I didn't separete the LCMS modes because too much methods will be hard to visualize the overlapping.

counts <- root_comp_long %>% 
  group_by(Method) %>% 
  summarize(n = n())
#count the total features for include that in legend.

# Merge counts with root_comp_long
root_comp_long <- left_join(root_comp_long, counts, by = "Method")

root_comp_long$Method <- paste0(root_comp_long$Method, " (n = ", root_comp_long$n, ")")
legend_order <- c("NMR_weight (n = 60)", "LCMS_mz (n = 901)", "MALDI_mz (n = 3089)", "FTICR_NeutralMass (n = 31939)")
# there must be a way to automatic it but I am just being lazy.

plot_root_comp <- ggplot(root_comp_long, aes(x = MolecularWeight, fill = Method)) +
  geom_histogram(data = filter(root_comp_long, Method %in% c("LCMS_mz (n = 901)")),
                 position = "identity", bins = 83, alpha = 0.5, aes(y = after_stat(count))) +
  geom_histogram(data = filter(root_comp_long, Method %in% c("FTICR_NeutralMass (n = 31939)")),
                 position = "identity", bins = 191, alpha = 0.5, aes(y = after_stat(count) / 10)) +
  geom_histogram(data = filter(root_comp_long, Method %in% c("MALDI_mz (n = 3089)")),
                 position = "identity", bins = 31, alpha = 0.5, aes(y = after_stat(count) / 10)) +
  geom_histogram(data = filter(root_comp_long, Method %in% c("NMR_weight (n = 60)")),
                 position = "identity", bins = 17, alpha = 0.5, aes(y = after_stat(count))) +
  #need to plot separately for double y axis, be aware of that "aes(y = ..count.. / 10)". The number need to match with the "trans" in sec.axis
  scale_fill_manual(values = c("LCMS_mz (n = 901)" = "blue", "NMR_weight (n = 60)" = "green","MALDI_mz (n = 3089)" = "yellow", "FTICR_NeutralMass (n = 31939)" = "red"), limits = legend_order) +  
  # Specify colors, RGB works good for overlapped color. 
  labs(title = "Histogram of Root metabolites",
       x = "Molecular Weight",
       y = "Frequency") +
  cowplot::theme_minimal_grid() +
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  scale_y_continuous(
    name = "NMR and LCMS frequency",
    sec.axis = sec_axis(trans = ~.*10, name = "MALDI and FTICR Frequency"),  # Specify labels for secondary axis
  ) 


rhizo_comp <- read.csv("Rhizo_peak_comp.csv")
### rhizo
rhizo_comp_long <- tidyr::gather(rhizo_comp, key = "Method", value = "MolecularWeight", na.rm = TRUE)

rhizo_comp_long <- rhizo_comp_long %>%
  mutate(Method = case_when(
    Method %in% c("RPPOS_mz", "RPNEG_mz", "HPOS_mz", "HNEG_mz") ~ "LCMS_mz",
    TRUE ~ Method
  ))

counts <- rhizo_comp_long %>% 
  group_by(Method) %>% 
  summarize(n = n())

# Merge counts with rhizo_comp_long
rhizo_comp_long <- left_join(rhizo_comp_long, counts, by = "Method")

rhizo_comp_long$Method <- paste0(rhizo_comp_long$Method, " (n = ", rhizo_comp_long$n, ")")
unique(rhizo_comp_long$Method)

legend_order <- c("NMR_weight (n = 48)", "LCMS_mz (n = 1313)", "FTICR_NeutralMass (n = 30171)")

plot_rhizo_comp <- ggplot(rhizo_comp_long, aes(x = MolecularWeight, fill = Method)) +
  geom_histogram(data = filter(rhizo_comp_long, Method == "LCMS_mz (n = 1313)"),
                 position = "identity", bins = 61, alpha = 0.5, aes(y = ..count..)) +
  geom_histogram(data = filter(rhizo_comp_long, Method == "FTICR_NeutralMass (n = 30171)"),
                 position = "identity", bins = 83, alpha = 0.5, aes(y = ..count.. / 10)) +
  geom_histogram(data = filter(rhizo_comp_long, Method == "NMR_weight (n = 48)"),
                 position = "identity", bins = 17, alpha = 0.5, aes(y = ..count..)) +
  scale_fill_manual(values = c("NMR_weight (n = 48)" = "green","LCMS_mz (n = 1313)" = "blue", "FTICR_NeutralMass (n = 30171)" = "red"),limits = legend_order) +  # Specify colors
  labs(title = "Histogram of rhizo metabolites",
       x = "Molecular Weight",
       y = "Frequency") +
  cowplot::theme_minimal_grid() +
  #theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  scale_y_continuous(
    name = "NMR and LCMS frequency",
    sec.axis = sec_axis(trans = ~.*10, name = "FTICR Frequency"),  # Specify labels for secondary axis  # Set limits for the primary axis
  ) 

ggsave("rhizo_comp.pdf",plot_rhizo_comp, height =4, width = 8)
ggsave("root_comp.pdf",plot_root_comp, height =4, width = 8)
