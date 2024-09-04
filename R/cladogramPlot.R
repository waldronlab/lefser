library(tidyverse)
library(proxy)
suppressPackageStartupMessages(library(dendextend))
library(circlize)
data(zeller14)
zeller14 <- zeller14[, zeller14$study_condition != "adenoma"]
res <- lefser(zeller14, groupCol = "study_condition", blockCol = "age_category")

# format organism's taxa for use in cladogram() function
# results in a string with organism's taxa names separated by single space
format_full_organism_name <- function(full_organism_name, separator= NULL){
  names_split_tovector <- unlist(strsplit(full_organism_name, separator))
  sample_name_formatted <- paste(names_split_tovector, collapse = " ")
  sample_name <- tail(names_split_tovector, n=1)
  names(sample_name_formatted) <- sample_name
  return(sample_name_formatted)
}
# format_full_organism_name() function is used with lapply to process all   
get_all_names <- function(names_vector, separator= NULL){
  formated_names <- unlist(lapply(names_vector, format_full_organism_name, separator))
  df <- data.frame(organism_name=names(formated_names), full_organism_name=formated_names)
  return(df)
}

# taxonomic dendrogram 
calc_evolutionary_distances <- function(df){
  distances <- 
    df %>%
    separate_rows(full_organism_name, sep=" ") %>%
    mutate(has_ID = 1) %>%
    pivot_wider(names_from = full_organism_name, values_from = has_ID, values_fill = list(has_ID = 0)) %>%
    column_to_rownames("organism_name") %>%
    proxy::dist(by_rows = TRUE, method = "Jaccard")
  return(distances)
}

#cladogram
cladogram <- function(res, separator=NULL){
  
  if(is.null(separator)){
    df_names <- get_all_names(res$Names, "\\|")
  } else{
    df_names <- get_all_names(res$Names,separator)
  }
  
  evo_dist <- calc_evolutionary_distances(df_names)
  dend=as.dendrogram(hclust(evo_dist))
  labels = labels(dend)  # name of birds
  n = length(labels)  # number of bird species
  
  dend=dend%>%
    set("branches_lwd",1) %>%
    set("branches_lty",1) %>%
    set("leaves_pch", c(19)) %>% 
    set("leaves_cex", abs(res$scores[match(labels(dend),a$organism_name)])/2) %>% 
    set("leaves_col", ifelse(sign(res$scores[match(labels(dend),a$organism_name)])==-1, "red","darkgreen"))%>% 
    set("labels_cex", c(.4))%>%
    set("labels_colors", ifelse(sign(res$scores[match(labels(dend),a$organism_name)])==-1, "red","darkgreen"))
  
  circos.initialize("a", xlim = c(0, n)) # only one sector
  circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.3, 
               panel.fun = function(x, y) {
                 for(i in seq_len(n)) {
                   circos.text(i-0.5, 0, labels[i], adj = c(0, 0.5), 
                               facing = "clockwise", niceFacing = TRUE,
                               col = ifelse(sign(res$scores[match(labels[i],a$organism_name)])==-1, "red","darkgreen"),cex = 0.5)
                 }
               })

  dend_height = attr(dend, "height")
  circos.track(ylim = c(0, dend_height), bg.border = NA, 
               track.height = 0.4, panel.fun = function(x, y) {
                 circos.dendrogram(dend)
               })
}
  

  