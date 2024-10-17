library(ggplot2)
library(ggalluvial)
library(dplyr)
library(ggrepel)
library(wesanderson)
library(ComplexHeatmap)

#-----------
# version 1.1
# added type check for clone_id
# version 1.2
# some minor performance optimization
#-----------


clone_expansion_donut <- function(sample_name, clone_exp){
  #sample_name = 'ISAC99_1' # for debug purpose
  if (typeof(clone_exp$clone_id) != 'character'){
    clone_exp <- clone_exp %>% mutate(clone_id = as.character(clone_id))
  }
  df <- clone_exp %>%
    filter(Sample_Name == sample_name) %>%
    mutate(clone_id = case_when(
      clone_count == 1 ~ 'unique',
      clone_count >1 ~ clone_id
    )) %>%
    arrange(clone_count)
  #df <- clone_exp %>% 
  #  filter(Sample_Name == sample_name) %>%
  #  arrange(clone_count)
  df <- rbind(tibble(CDR3_concat='unique', 
                     Sample_Name = sample_name, 
                     clone_count = dim(df%>%filter(clone_id=='unique'))[1], 
                     clone_id = 'unique'),
              df %>% filter(clone_id != 'unique')) %>%
    mutate(label = factor(case_when(
      clone_id == 'unique' ~ 'Unique',
      clone_id != 'unique' & clone_count == 2 ~ '2',
      clone_id != 'unique' & clone_count == 3 ~ '3',
      clone_id != 'unique' & clone_count == 4 ~ '4',
      clone_id != 'unique' & clone_count == 5 ~ '5',
      clone_id != 'unique' & clone_count > 5 & clone_count <=9 ~ '6-9',
      clone_id != 'unique' & clone_count > 9 & clone_count <= 19 ~ '10-19',
      clone_id != 'unique' & clone_count > 19 ~ '>=20'
    ), levels=c('Unique', '2', '3', '4', '5', '6-9', '10-19', '>=20')),
    fraction = clone_count/sum(clone_count),
    ymax = cumsum(fraction),
    ymin = c(0, head(ymax, n=-1)),
    )
  breaks <- df %>% group_by(label) %>% summarise(label_frac=sum(fraction)) %>% 
    mutate(max = cumsum(label_frac), 
           min = c(0,head(max,n=-1)),
           pos = (max+min)/2
    )
  # make a named vector for color
  cols = c('#FFFFFF', wes_palette("Rushmore1", nlevels(breaks$label)-1, type = "continuous"))
  names(cols) = levels(breaks$label)
  g <- ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=2.7, fill=label)) +
    geom_rect(color = 1, size=0.2) +
    annotate("text", x=2, y=0, label = sum(df$clone_count), size = 8) +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    scale_y_continuous(breaks = breaks$pos, labels = breaks$label) +
    scale_fill_manual(values = cols) +
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
          axis.text = element_text(size = 12), 
          legend.position = "none",
          panel.background = element_rect(fill = "white"),
          plot.title = element_text(hjust = 0.5)) +
    labs(title = sample_name)
  return(g)
}

get_expanded_id <- function(clone_exp, thre=1){
  return (
    (clone_exp %>% filter(clone_count > thre) %>% select(clone_id))$clone_id
    )
}

get_top_expansion_id <- function(clone_exp, top_number=10){
  top_clone_id <- clone_exp %>% 
    group_by(clone_id) %>% 
    summarise(max_count = max(clone_count)) %>% 
    filter(max_count > 1) %>%
    slice_max(order_by = max_count, n = top_number) %>% 
    select(clone_id) %>% 
    unlist(use.names = F)
  return(top_clone_id)
}

get_large_expansion_id <- function(clone_exp, clone_thre=10){
  top_clone_id <- clone_exp %>% 
    group_by(clone_id) %>% 
    summarise(max_count = max(clone_count)) %>% 
    filter(max_count >= clone_thre) %>% 
    select(clone_id) %>% 
    unlist(use.names = F)
  return(top_clone_id)
}

get_clone_exp_sub <- function(data_scTCR, sub_name){
  clone_exp_sub <- data_scTCR %>%
    filter(cell_type == sub_name) %>%
    # mutate(Sample_Name = paste0(proj_id, '_', Sample_Name)) %>% # when proj_wise plots are needed
    group_by(CDR3_concat, Sample_Name) %>%
    summarise(clone_count = n()) %>%
    ungroup()
}

clone_exp_alluvium_preparation <- function(individual_name, clone_exp, top_mod){
  if (typeof(clone_exp$clone_id) != 'character'){
    clone_exp <- clone_exp %>% mutate(clone_id = as.character(clone_id))
  }
  clone_exp <- clone_exp  %>%
    tidyr::separate(Sample_Name, into = c('individual', 'time_point')) %>%
    filter(individual == individual_name)
  # QC
  failed_exp <- clone_exp %>% group_by(time_point) %>% tally() %>% filter(n<10) %>% select(time_point) %>% unlist(use.names = F)
  clone_exp <- clone_exp %>% filter(!time_point %in% failed_exp)
  #
  if (top_mod == 'union'){
    top_clone_id <- union(get_top_expansion_id(clone_exp), get_large_expansion_id(clone_exp))
  } else if (top_mod == 'top'){
    top_clone_id <- get_top_expansion_id(clone_exp, n)
  } else if (top_mod == 'large'){
    top_clone_id <- get_large_expansion_id(clone_exp, n)
  } else {
    stop('unkown top mod')
  }
  
  df <- clone_exp %>%
    group_by(time_point) %>%
    mutate(clone_ratio = clone_count/sum(clone_count)) %>%
    filter(clone_id %in% top_clone_id) %>%
    mutate(relative_clone_ratio = clone_ratio/sum(clone_ratio)) %>%
    mutate(time_point = factor(time_point, levels = gtools::mixedsort(unique(clone_exp$time_point)))) %>%
    ungroup()
  return(df)
}


clone_expansion_alluvium<- function(individual_name, clone_exp, clone_label = F, top_mod = 'union', n=10){
  df <- clone_exp_alluvium_preparation(individual_name, clone_exp, top_mod)
  
  g <- ggplot(df, aes(x = time_point, stratum = clone_id, alluvium = clone_id, 
                      #y= clone_ratio, 
                      y = relative_clone_ratio,
                      fill = clone_id, color=clone_id)) +
    geom_alluvium() +
    geom_stratum(size = 0.1) +
    scale_fill_manual(values = wes_palette("Rushmore1", length(unique(df$clone_id)), type = "continuous"))+
    theme(legend.position = "none")+
    labs(title = individual_name)
  if (clone_label) {
    g <- g + geom_text_repel(aes(label = clone_id), stat = 'stratum', size = 4)
  }
  return(g)
}

trend_determination_plot <- function(individual_name, clone_exp, top_mod='union'){

  df <- clone_exp_alluvium_preparation(individual_name, clone_exp, top_mod)

  df_complete <- df%>%select(time_point, clone_id,relative_clone_ratio)%>%
    ungroup()%>%
    tidyr::complete(time_point, clone_id,fill = list(relative_clone_ratio=0))
  tmp <- tidyr::pivot_wider(df_complete,names_from = time_point, values_from = relative_clone_ratio)
  tmp2 <- tmp %>% select(-clone_id) %>% as.data.frame()
  rownames(tmp2) <- tmp$clone_id
  tmp2 <- t(scale(t(tmp2)))
  g <- Heatmap(tmp2, show_column_dend = F, column_order = c(levels(df_complete$time_point)))
  return(g)
}
