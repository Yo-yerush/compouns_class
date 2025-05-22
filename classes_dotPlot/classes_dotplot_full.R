library(dplyr)

#' comp_classes
#'
#' description
#' @param class_table_for_metaboanalyst  full path to "class_table_for_metaboanalyst.csv"
#' @param group1 character, treatment samples name as wright in "table_for_metaboanalyst" file
#' @param group2 character, control samples name as wright in "table_for_metaboanalyst" file
#' @param combine_flavonoids [optional] logical (default: FALSE), if TRUE - combine all flavonoids subgroups in the output dotplot
#' @param print_t.test_table [optional] logical (default: FALSE), if TRUE - print table that contains pathway_classes and t.test results
#' @param print_plot_table [optional] logical (default: FALSE), if TRUE - print table that contains the output dotplot values
#' @return none. side effect : writing merged file - ???
#'
#' @example comp_classes(class_table_for_metaboanalyst = "C:/user/lcms/raw_file.csv",
#'                       group1 = "mut1",
#'                       group2 = "wt")
#'                          

#' @export comp_classes



classes_dotplot <- function(class_table_for_metaboanalyst,
                         group1 = NULL, #treatment
                         group2 = NULL, #control
                         p.value_type = "p.value", # FDR
                         combine_flavonoids = F,
                         print_t.test_table = F,
                         print_plot_table = F)
  {
  
  suppressWarnings(library(dbplyr))
  suppressWarnings(library(ggplot2))

  if (is.data.frame(class_table_for_metaboanalyst)) {
    df <- class_table_for_metaboanalyst
  } else {
    df = read.csv(class_table_for_metaboanalyst)
  }
  names(df) = df[1,]
  names(df)[1] = "Name"
  df = df[-1,]

  if (is.null(group1) & is.null(group2)) {
    group1 = unique(names(df)[-1])[1]
    group2 = unique(names(df)[-1])[2]
  }
  
  pvalue_vec = NA
  fdr_vec = NA
  ttest_vec = NA
  FC_vec = NA
  i=1
  while (i <= nrow(df)) {
    pvalue_vec[i] = t.test(as.numeric(df[i,grep(group1,names(df))]),
                           as.numeric(df[i,grep(group2,names(df))]), var.equal = T)$p.value
    
    ttest_vec[i] = t.test(as.numeric(df[i,grep(group1,names(df))]),
                          as.numeric(df[i,grep(group2,names(df))]), var.equal = T)$statistic
    
    FC_vec[i] = mean(as.numeric(df[i,grep(group1,names(df))]))/mean(as.numeric(df[i,grep(group2,names(df))]))
    
    i=i+1
  }
  log2FC_vec = log2(FC_vec)
  fdr_vec = p.adjust(pvalue_vec, "fdr", length(pvalue_vec))
  
  
  df_2.view = data.frame(Name = df[,"Name"],
                         Fold_Change = FC_vec,
                         log2FC = log2FC_vec,
                         p.value = pvalue_vec,
                         FDR = fdr_vec,
                         t = ttest_vec)
  if (print_t.test_table == T) {
    return(df_2.view = df_2.view)
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
    }
  
  df_2 = df_2.view
  df_2$Name = gsub(" X.*","",df_2$Name)
  
  # change and unite names
#  df_2$Name = gsub("tannin","HTs",df_2$Name)
#  df_2$Name = gsub("gamm.*","amino acid",df_2$Name)
#  df_2$Name = gsub("amino acid zwitterion","amino acid",df_2$Name)
  df_2$Name = gsub("aromatic ketone","aromatic compound",df_2$Name)
  df_2$Name = gsub("benzenoid aromatic compound","aromatic compound",df_2$Name)
  df_2$Name = gsub("phenol","aromatic compound",df_2$Name)
  df_2$Name = gsub("aromatic compounds","aromatic compound",df_2$Name)
  if (combine_flavonoids == T) {
    df_2$Name = gsub("flav.*","flavonoids",df_2$Name)
    df_2$Name = gsub("isoflav.*","flavonoids",df_2$Name)
    df_2$Name = gsub("proanthocyanidin","flavonoids",df_2$Name)
  }

  # p.value type to use
  if (p.value_type == "p.value") {
    p.value_vector = df_2$p.value
  } else if (p.value_type == "fdr") {
    p.value_vector = df_2$FDR
  } #else {
    #errorCondition("change p.value_type to 'p.value' or 'fdr'")
  #}
  df_sig = df_2[p.value_vector < 0.05,]
  df_up = df_sig[df_sig$t > 0,]
  df_down = df_sig[df_sig$t < 0,]
  
  
  plot_df = data.frame(pathway = unique(df_sig$Name),
                       count_total = "",
                       count_total_sig = "",
                       count_up = "",
                       count_down = "",
                       ratio_up = "",
                       ratio_down = "",
                       ratio_sig = "")
  
  i=1
  while (i <= length(plot_df$pathway)) {
    plot_df$count_total[i] = length(df_2[df_2$Name == plot_df$pathway[i],1])
    i=i+1
  }
  
  i=1
  while (i <= length(plot_df$pathway)) {
    plot_df$count_total_sig[i] = length(df_sig[df_sig$Name == plot_df$pathway[i],1])
    i=i+1
  }
  
  i=1
  while (i <= length(plot_df$pathway)) {
    plot_df$count_up[i] = length(df_up[df_up$Name == plot_df$pathway[i],1])
    i=i+1
  }
  
  i=1
  while (i <= length(plot_df$pathway)) {
    plot_df$count_down[i] = as.numeric(length(df_down[df_down$Name == plot_df$pathway[i],1]))
    i=i+1
  }
  
  ##
  i=2
  while (i <= length(names(plot_df)[1:5])) {
    plot_df[,i] = as.numeric(plot_df[,i])
    i=i+1
  }
  ##
  
  i=1
  while (i <= length(plot_df$pathway)) {
    plot_df$ratio_up[i] = plot_df$count_up[i]/plot_df$count_total_sig[i]
    i=i+1
  }
  
  #
  
  i=1
  while (i <= length(plot_df$pathway)) {
    plot_df$ratio_down[i] = plot_df$count_down[i]/plot_df$count_total_sig[i]
    i=i+1
  }
  
  #
  df_name_without_X = gsub(" X[0-9].*","",df$Name)
  i=1
  while (i <= length(plot_df$pathway)) {
    plot_df$ratio_sig[i] = (length(df_sig[df_sig$Name == plot_df$pathway[i],1])/
      length(df_name_without_X[df_name_without_X == plot_df$pathway[i]]))
    i=i+1
  }
  
  i=5
  while (i <= length(names(plot_df))) {
    plot_df[,i] = as.numeric(plot_df[,i])
    i=i+1
  }
  
  #

  
  ###
  
  # remove unknown comp
  plot_df = plot_df[grep("[A-z]", plot_df$pathway),]
  if (print_plot_table == T) {
    return(plot_df = plot_df)
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }
  ###
  
  plot_df %>% 
    ggplot(aes(count_total_sig, reorder(pathway,(count_total_sig)), size = count_total, color = ratio_up))+#, shape = Flavenoid)) + 
    #guides(size = guide_legend(title = "Significants\n   (p<0.05)")) + theme_classic() +
    guides(size = guide_legend(title = "Total\nCompounds")) + theme_classic() +
    #  scale_color_gradient2("ratio_up", low = "#1e5fa3", high = "#c33a3b", mid = "#ecded6", midpoint = 0.4) + 
    scale_color_gradientn(paste0("Ratio to\n",group1),
                          colours = c("#134b87","#2066ab","#3f8ebf","#83bcd9","#d1e4ee","#e9dfd8","#f8bfa2","#ec9374","#d86450","#9b1127","#7e0722"),) + #theme_classic() +
    labs(title = paste0(group1, " vs ", group2), x = "Significant Compounds\n(p<0.05)", y = "") +
    theme(plot.title=element_text(hjust=0.5, face = "bold"),
          legend.key.size = unit(0.45, 'cm'),
          legend.title = element_text(size=9.5),
          axis.text = element_text(face="bold")) + 
    geom_point()
  
  ###
  
}
