#load library
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(DescTools)
setwd("~/INVA")
### peak overlapping --------
root_comp <- read.csv("peak_comp_root_neutral.csv")
rhizo_comp <- read.csv("peak_comp_rhizo_neutral.csv")
#root_comp is a dataframe with the weights of all features by each method, with method as columns


### root_neg -------
root_comp_long <- tidyr::gather(root_comp, key = "Method", value = "neutral mass", na.rm = TRUE)

root_comp_long <- root_comp_long %>%
  mutate(Method = case_when(
    Method %in% c("RPNEG", "HNEG", "RPPOS", "HPOS") ~ "LCMS",
    TRUE ~ Method
  ))
#I didn't separete the LCMS modes because too much methods will be hard to visualize the overlapping.

root_comp_long <- root_comp_long %>%
  mutate(Method = case_when(
    Method %in% c("DHB", "NEDC", "APEBA") ~ "MALDI",
    TRUE ~ Method
  ))

counts <- root_comp_long %>% 
  group_by(Method) %>% 
  summarize(n = n())
#count the total features for include that in legend.

range_and_median <- root_comp_long %>%
  group_by(Method) %>%
  summarize(
    Range = paste0("[", min(`neutral mass`), ", ", max(`neutral mass`), "]"),
    Median = median(`neutral mass`)
  )


# Merge counts with root_comp_long
legend_order <- c("NMR", "LCMS", "MALDI", "FTICR")
# there must be a way to automatic it but I am just being lazy.

plot_root_comp<- ggplot(root_comp_long, aes(x = `neutral mass`, fill = Method)) +
  geom_histogram(data = filter(root_comp_long, Method %in% c("FTICR")),
                 position = "identity", bins = 97, alpha = 0.5, aes(y = after_stat(count))) + 
  geom_histogram(data = filter(root_comp_long, Method %in% c("LCMS")),
                 position = "identity", bins = 53, alpha = 0.5, aes(y = after_stat(count))) +
  geom_histogram(data = filter(root_comp_long, Method %in% c("MALDI")),
                 position = "identity", bins = 41, alpha = 0.5, aes(y = after_stat(count))) +
  geom_histogram(data = filter(root_comp_long, Method %in% c("NMR")),
                 position = "identity", bins = 37, alpha = 0.5, aes(y = after_stat(count))) +
  #need to plot separately for double y axis, be aware of that "aes(y = ..count.. / 10)". The number need to match with the "trans" in sec.axis
  scale_fill_manual(values = c("LCMS" = "blue", "NMR" = "green","MALDI" = "yellow", "FTICR" = "red"), limits = legend_order) +  
  # Specify colors, RGB works good for overlapped color. 
  labs(title = "Histogram of Root metabolites by Neutral mass",
       x = "Neutral mass (Da)",
       y = "Frequency") +
  cowplot::theme_minimal_grid() +
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  scale_y_continuous(
    name = "Frequency"
  ) 

### root_pos -------
root_comp_long <- tidyr::gather(root_comp_pos, key = "Method", value = "neutral mass", na.rm = TRUE)

root_comp_long <- root_comp_long %>%
  mutate(Method = case_when(
    Method %in% c("RPPOS", "HPOS") ~ "LCMS (pos)",
    TRUE ~ Method
  ))
#I didn't separete the LCMS modes because too much methods will be hard to visualize the overlapping.

root_comp_long <- root_comp_long %>%
  mutate(Method = case_when(
    Method %in% c("DHB", "NEDC", "X4.APEBA") ~ "MALDI (pos)",
    TRUE ~ Method
  ))

counts <- root_comp_long %>% 
  group_by(Method) %>% 
  summarize(n = n())
#count the total features for include that in legend.

range_and_median <- root_comp_long %>%
  group_by(Method) %>%
  summarize(
    Range = paste0("[", min(`neutral mass`), ", ", max(`neutral mass`), "]"),
    Median = median(`neutral mass`)
  )


# Merge counts with root_comp_long
legend_order <- c("NMR", "LCMS (pos)", "MALDI (pos)")
# there must be a way to automatic it but I am just being lazy.

plot_root_comp_pos <- ggplot(root_comp_long, aes(x = `neutral mass`, fill = Method)) +

  geom_histogram(data = filter(root_comp_long, Method %in% c("LCMS (pos)")),
                 position = "identity", bins = 57, alpha = 0.5, aes(y = after_stat(count))) +
  geom_histogram(data = filter(root_comp_long, Method %in% c("MALDI (pos)")),
                 position = "identity", bins = 37, alpha = 0.5, aes(y = after_stat(count) / 10)) +
  geom_histogram(data = filter(root_comp_long, Method %in% c("X4.APEBA")),
                 position = "identity", bins = 41, alpha = 0.5, aes(y = after_stat(count) / 10)) + 
  geom_histogram(data = filter(root_comp_long, Method %in% c("NMR")),
                 position = "identity", bins = 31, alpha = 0.5, aes(y = after_stat(count))) +
  #need to plot separately for double y axis, be aware of that "aes(y = ..count.. / 10)". The number need to match with the "trans" in sec.axis
  scale_fill_manual(values = c("LCMS (pos)" = "blue", "NMR" = "green","MALDI (pos)" = "yellow", "X4.APEBA" = "orange"), limits = legend_order) +  
  # Specify colors, RGB works good for overlapped color. 
  labs(title = "Histogram of Root metabolites postive mode",
       x = "m/z ratio",
       y = "Frequency") +
  cowplot::theme_minimal_grid() +
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  scale_y_continuous(
    name = "NMR and LCMS frequency",
    sec.axis = sec_axis(trans = ~.*10, name = "MALDI Frequency"),  # Specify labels for secondary axis
  ) 


### rhizo neg -------

rhizo_comp_long <- tidyr::gather(rhizo_comp, key = "Method", value = "neutral mass", na.rm = TRUE)

rhizo_comp_long <- rhizo_comp_long %>%
  mutate(Method = case_when(
    Method %in% c("RPNEG", "HNEG", "RPPOS", "HPOS") ~ "LCMS",
    TRUE ~ Method
  ))

range_and_median <- rhizo_comp_long %>%
  group_by(Method) %>%
  summarize(
    Range = paste0("[", min(`neutral mass`), ", ", max(`neutral mass`), "]"),
    Median = median(`neutral mass`)
  )


counts <- rhizo_comp_long %>% 
  group_by(Method) %>% 
  summarize(n = n())

# Merge counts with rhizo_comp_long
legend_order <- c("NMR", "LCMS", "FTICR")

plot_rhizo_comp <- ggplot(rhizo_comp_long, aes(x = `neutral mass`, fill = Method)) +

  geom_histogram(data = filter(rhizo_comp_long, Method == "FTICR"),
                 position = "identity", bins = 97, alpha = 0.5, aes(y = ..count..)) +
  geom_histogram(data = filter(rhizo_comp_long, Method == "LCMS"),
                 position = "identity", bins = 53, alpha = 0.5, aes(y = ..count..)) +
  geom_histogram(data = filter(rhizo_comp_long, Method == "NMR"),
                 position = "identity", bins = 37, alpha = 0.5, aes(y = ..count..)) +
  scale_fill_manual(values = c("NMR" = "green","LCMS" = "blue", "FTICR" = "red"),limits = legend_order) +  # Specify colors
  labs(title = "Histogram of rhizo metabolites by neutral mass",
       x = "Neutral mass (Da)",
       y = "Frequency") +
  cowplot::theme_minimal_grid() +
  #theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  scale_y_continuous(
    name = "Frequency",  # Specify labels for secondary axis  # Set limits for the primary axis
  ) 

# Merge counts with rhizo_comp_long
legend_order <- c("NMR_pos", "LCMS (pos)")

plot_rhizo_comp_pos <- ggplot(rhizo_comp_long, aes(x = `m/z`, fill = Method)) +
  geom_histogram(data = filter(rhizo_comp_long, Method == "LCMS (pos)"),
                 position = "identity", bins = 61, alpha = 0.5, aes(y = ..count..)) +
  geom_histogram(data = filter(rhizo_comp_long, Method == "NMR_pos"),
                 position = "identity", bins = 17, alpha = 0.5, aes(y = ..count..)) +
  scale_fill_manual(values = c("NMR_pos" = "green","LCMS (pos)" = "blue"),limits = legend_order) +  # Specify colors
  labs(title = "Histogram of rhizo metabolites positive mode",
       x = "m/z ratio",
       y = "Frequency") +
  cowplot::theme_minimal_grid() +
  #theme_minimal() +
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  scale_y_continuous(
    name = "NMR and LCMS frequency"  # Specify labels for secondary axis  # Set limits for the primary axis
  ) 

ggsave("rhizo_comp_neg.pdf",plot_rhizo_comp_neg, height =5, width = 12)
ggsave("root_comp_neg.pdf",plot_root_comp_neg, height =5, width = 12)

ggsave("rhizo_comp_pos.pdf",plot_rhizo_comp_pos, height =5, width = 12)
ggsave("root_comp_pos.pdf",plot_root_comp_pos, height =5, width = 12)

ggsave("rhizo_comp_neu.pdf",plot_rhizo_comp, height =5, width = 12)
ggsave("root_comp_neu.pdf",plot_root_comp, height =5, width = 12)
