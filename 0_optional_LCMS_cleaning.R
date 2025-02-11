library(tidyverse)

find_ppm <- function(mo, me){
  ppm <- (mo - me)/me * 10^6
  return(ppm)
}

read_mgf <- function(path){
  #Read in the mgf file
  l <- readLines(path)
  #Find the start of each feature
  ti <- which(str_detect(l, 'FEATURE_ID'))
  #Find the end of each feature
  tt <- which(str_detect(l, 'END IONS'))
  ck <- data.frame()
  #Pull out the necessary information
  for(i in 1:length(ti)){
    tmp <- tibble(Compound_ID = str_split(l[ti[i]], '=')[[1]][2],
                  ions = list(as.numeric(sapply(str_split(l[(ti[i]+8):(tt[i]-1)], ' '), function(x) x[1]))))
    ck <- rbind(ck, tmp)
  }
  # #Check to see if "FT_" is at the start of feature names
  # if(str_detect(ck$Compound_ID[1], 'FT_', negate = TRUE)){
  #   #If not, add it on
  # frags <- ck %>%
  #   modify_at('Compound_ID', ~paste0('FT_', .x))
  # return(frags)
  # }else{
  #   return(frags)
  # }
  return(ck)
}

setwd("~/PERMA")
tm <- read_mgf('PERMA_RPPOS.mgf')

tm <- read_mgf('hneg_rhizo.mgf')

duplicate_ids <- tm %>%
  dplyr::group_by(Compound_ID) %>%
  dplyr::filter(dplyr::n() > 1) %>%
  dplyr::arrange(Compound_ID)


find_ISFrags <- function(link_data, intensity_data, mgf_file,
                         time_tol = 0.005, ppm_tol = 5, cor_cutoff = 0.9){
  
  #Check data integrity:
  if(sum(colnames(link_data) %in% c('Compound_ID', 'm/z', 'RT [min]')) != 3){
    stop('Check Link Data, it must contain the following columns named exactly: Compound_ID, RT [min], m/z')
  }
  
  #Format the fragment table to be correct
  fragtab <- link_data %>%
    dplyr::select(Compound_ID, `m/z`, `RT [min]`) %>%
    dplyr::rename('mz' = `m/z`, 'rt' = `RT [min]`) %>%
     tidyr::drop_na()
  
  #Calculate the difference between each rt and all others:
  #Output the results into a list
  olist <- list()
  for(i in 1:nrow(fragtab)){
    ovec <- fragtab$rt - fragtab$rt[i]
    names(ovec) <- fragtab$Compound_ID
    olist[[i]] <- ovec
  }
  
  #Now trim each vector to only places where the rt dif <= 1 and pull the names
  olist_clean <- lapply(olist, function(x) names(x[abs(x) <= time_tol]))
  
  #Remove features that are a length of 1 - These have no matches
  olist_trimmed <- olist_clean[sapply(olist_clean, length) != 1]
  #And each set will show up twice so trim it down
  olist_final <- olist_trimmed[!duplicated(olist_trimmed)]
  
  #Now turn it into a dataframe
  #Making a unique group for each set
  dfo <- tibble()
  for(i in 1:length(olist_final)){
    df <- tibble(group = paste0('G', i),
                 matches = list(olist_final[[i]]))
    dfo <- rbind(dfo, df)
  }
  
  #Now make a table with the groups
  group_tab <- dfo %>%
    unnest('matches') %>%
    dplyr::rename('Compound_ID' = matches)
  
  #Make a table of intensities
  int_table <- intensity_data %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column('Compound_ID')
  
  #Make a table with a correlation matrix of each group
  big_group <- group_tab %>%
    #Join on teh intensity data
    left_join(int_table) %>%
    #Make long
    pivot_longer(cols = c(-group, -Compound_ID), names_to = 'sample') %>%
    #Group
    group_by(group) %>%
    nest() %>%
    #Make a correlation matrix for each set of candidates
    mutate(cor_tab = purrr::map(data, function(x){
      ct <- x %>%
        pivot_wider(names_from = Compound_ID) %>%
        column_to_rownames('sample') %>%
        dplyr::select_if(is.numeric) %>%
        cor()
      return(ct)
    })) 
  
  #Filter the candidates based on correlation info:
  IS_cand <- big_group %>%
    mutate(cands = purrr::map_chr(cor_tab, function(x){
      x %>%
        as.data.frame() %>%
        rownames_to_column('var1') %>%
        #Make the correlation table long for ease
        pivot_longer(starts_with('FT_'), names_to = 'var2') %>%
        #Remove the diag
        filter(var1 != var2) %>%
        #Grab only ones above the cutoff
        filter(value >= cor_cutoff) %>%
        #Pull out var1 since it will include both sides
        pull(var1) %>%
        #Remove duplicates
        unique() %>%
        paste0(collapse = '; ')
    })) %>%
    ungroup() %>%
    #Remove candidates where nothing passed
    filter(cands != '') %>%
    #Clear out the duplicates
    filter(!duplicated(cands))
  
  #Now using the grouping information prep for MS2 data
  grp <- IS_cand %>%
    #Clean up the fragment data
    mutate(split = str_split(cands, pattern = '; ')) %>%
    dplyr::select(group, split) %>%
    unnest('split') %>%
    #Rename for merging
    dplyr::rename('Compound_ID' = split) %>%
    #Merge in the fragment data
    left_join(fragtab) %>%
    group_by(group) %>%
    #Keep the highest mz as the parent
    mutate(type = ifelse(max(mz) == mz, 'parent', 'frag'))
  
  frags <- read_mgf(mgf_file)
  
  #Merge in the MS2 data
  grp_frag <- grp %>%
    left_join(frags) %>%
    #Grab the necessary information from the fragment tab
    reframe(parent_mz = mz[type == 'parent'],
            parent_ft = Compound_ID[type == 'parent'],
            frag_mz = list(mz[type == 'frag']),
            frag_fts = list(Compound_ID[type == 'frag']),
            ions_parent = ions[type == 'parent']) %>%
    #Calculate PPMs of each fragment candidate against the MS2 data of parent fragment
    mutate(ppm_vals = map2(frag_mz, ions_parent, function(x,y){
      purrr::map(x, function(a){
        find_ppm(a,y)
      })
    })) %>%
    mutate(match = purrr::map2(ppm_vals, frag_fts, function(x,y){
      lgls <- sapply(x, function(x) any(abs(x) <= ppm_tol))
      data.frame(frag_cand = y,
                 match = lgls)
    })) %>%
    dplyr::select(parent_ft, match) %>%
    unnest('match') %>%
    filter(match)
  
  final_output <- grp_frag %>%
    dplyr::select(-match) %>%
    group_by(parent_ft) %>%
    summarise(fragments = paste0(unique(frag_cand), collapse = '; '))
  
  return(final_output)
}

#link_data contains columns with m/z, RT [min.], and Compound_ID
#Intensity data in a numeric matrix of RAW (or adjusted for things like sample weight) of feature intensities
#Each column is a feature and each row is a sample with sample names as rownames ONLY
#mgf_file is the path of the mgf file you gave to SIRIUS
#time_tol is the tolerance for finding initial candidates - Defaults to 0.005
#cor_cutoff is the cutoff for correlations between two features to be considered a fragment - Default is 0.9
#ppm_tol is the ppm tolerance between candidate fragments (which pass the correlation test) and 
#the parent (heaviest) ion MS2 fragments - Default is 5

#Example use:


setwd("~/INVA/cleaning")

sample_types <- c("rhizo", "root")
modes <- c("hneg", "hpos", "rppos", "rpneg")

# Loop over each combination
for (sample_type in sample_types) {
  for (mode in modes) {
    # Generate the appropriate file name based on current loop iteration
    link_file <- paste0(mode, "_", sample_type, ".csv")
    mgf_file <- paste0(mode, "_", sample_type, ".mgf")
    
    # Read the link data and intensity data
    link <- read_csv(link_file)
    intensity <- read.csv(paste0(sample_type, "_LCMS_", mode, ".csv"))
    
    # Rename column in intensity data
    names(intensity)[names(intensity) == "Compound.ID"] <- "Compound_ID"
    
    # Identify duplicate IDs in link data
    duplicate_ids <- link %>%
      dplyr::group_by(Compound_ID) %>%
      dplyr::filter(dplyr::n() > 1) %>%
      dplyr::arrange(Compound_ID)
    
    # Prepare intensity data for the function
    intensity_new <- intensity %>%
      column_to_rownames('Compound_ID') %>%
      t()
    
    # Run the find_ISFrags function and assign to dynamically named variable
    frag_var_name <- paste0(sample_type, "_", mode, "_frag")
    assign(frag_var_name, find_ISFrags(link_data = link, 
                                       intensity_data = intensity_new, 
                                       mgf_file = mgf_file,
                                       time_tol = 0.005, 
                                       cor_cutoff = 0.98))
    
    # Print the name of the generated fragment variable (optional)
    print(paste("Created", frag_var_name))
  }
}

setwd('~/INVA/INVA')
metadata <- read.table("metadata.txt", header = TRUE, sep = "\t", row.names = 1)
metadata$block <- as.factor(metadata$block)
Rhizo_compounds <- read.csv("Rhizo_compound_link.csv", header = TRUE)
Root_compounds <- read.csv("Root_compound_link.csv", header = TRUE)

all_root_fragments <- c()

# Loop through each mode for root samples
for (mode in modes) {
  # Construct the variable name dynamically
  frag_var_name <- paste0("root_", mode, "_frag")
  
  # Extract fragments from the current root_[mode]_frag and clean up
  fragments <- get(frag_var_name) %>%
    pull(fragments) %>%
    str_split(";") %>%
    unlist() %>%
    str_trim()
  
  # Combine fragments from each mode into a single vector
  all_root_fragments <- c(all_root_fragments, fragments)
}

all_rhizo_fragments <- c()

# Loop through each mode for rhizo samples
for (mode in modes) {
  # Construct the variable name dynamically
  frag_var_name <- paste0("rhizo_", mode, "_frag")
  
  # Extract fragments from the current rhizo_[mode]_frag and clean up
  fragments <- get(frag_var_name) %>%
    pull(fragments) %>%
    str_split(";") %>%
    unlist() %>%
    str_trim()
  
  # Combine fragments from each mode into a single vector
  all_rhizo_fragments <- c(all_rhizo_fragments, fragments)
}

Root_compounds_filtered <- Root_compounds %>%
  filter(!CompoundID %in% all_root_fragments)

Rhizo_compounds_filtered <- Rhizo_compounds %>%
  filter(!compoundID %in% all_rhizo_fragments)

write.csv(Rhizo_compounds_filtered,"Rhizo_compounds_filtered.csv")
write.csv(Root_compounds_filtered, "Root_compounds_filtered.csv")


###remove gap filling
setwd('~/INVA/cleaning')
for (sample_type in sample_types) {
  for (mode in modes) {
    
    intensity <- read.csv(paste0(sample_type, "_LCMS_", mode, ".csv"))
    gapfill <- read.csv(paste0(sample_type, "_LCMS_", mode, "_gapfill.csv"))%>%
      select(-`m.z`, -`RT..min.`)
    
    colnames(gapfill) <- gsub(".*(Z\\d{2}.*)", "\\1", colnames(gapfill))
    colnames(gapfill) <- gsub(".*(R\\d{2}.*)", "\\1", colnames(gapfill))
    colnames(gapfill) <- gsub("^(([^_]*_){2}[^_]*).*", "\\1", colnames(gapfill))
  
    common_ids <-sort(intersect(gapfill$CompoundID, intensity$Compound.ID))
    common_columns <- intersect(names(gapfill), names(intensity))
    
    gapfill_aligned <- gapfill %>%
      filter(CompoundID %in% common_ids) %>%
      arrange(match(CompoundID, common_ids))%>%
      column_to_rownames("CompoundID")%>%
      select(all_of(common_columns))
    
    intensity_aligned <- intensity %>%
      filter(Compound.ID %in% common_ids) %>%
      arrange(match(Compound.ID, common_ids))%>%
      column_to_rownames("Compound.ID")%>%
      select(all_of(common_columns))
    
    gapfill_aligned[] <- apply(gapfill_aligned, 2, function(col) {
      # Replace values "8" and "32" with NA for later reference
      col[col == "8" | col == "32"] <- NA
      return(col)
    })
    
    # Now set the corresponding positions in intensity_aligned to 0
    for (i in 1:nrow(gapfill_aligned)) {
      for (j in 1:ncol(gapfill_aligned)) {
        if (is.na(gapfill_aligned[i, j])) {
          intensity_aligned[i, j] <- 0
        }
      }
    }
    if (sample_type == "root") {
      intensity_aligned_filtered <- intensity_aligned[rownames(intensity_aligned) %in% Root_compounds_filtered$CompoundID, ]
    } else if (sample_type == "rhizo") {
      intensity_aligned_filtered <- intensity_aligned[rownames(intensity_aligned) %in% Rhizo_compounds_filtered$compoundID, ]
    }
    final_df_name <- paste0(sample_type, "_LCMS_", mode, "_clean")
    assign(final_df_name, intensity_aligned_filtered)
    write.csv(intensity_aligned_filtered, paste0(final_df_name, ".csv"), row.names = TRUE)
  }
}
