classes_corr <- function(classDF, comp_pathways = NULL) {
  
  #classDF = "P:/yonatan/RNAseq_yonatan_2021/lcms/lcms_results_060722/peel/peel_files/files_for_metaboanalyst_18072022/pathway.name.for.metabo.18072022.csv",
  #comp_pathways = c("HTs","Anthocyanidin")
  
  # https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html
  
  if (is.data.frame(classDF)) {
    file1 <- classDF[-1,]
  } else {
    file1 = read.csv(classDF)[-1,]
  }
  class_names = gsub(" X[0-9]+","", file1[,1])
  row.names(file1) = file1[,1]
  file1 = file1[,-1]
  mat1 = apply(as.matrix(file1),1,as.numeric)
  cor_file1 = cor(mat1)
  
  library(ComplexHeatmap)
  #ht1 = Heatmap(cor_file1, name = "corralation",show_row_names = F, show_column_names = F)
  #ht2 = Heatmap(class_names, name = "Classes", width = 20,bottom_annotation = T)
  #ht1 + ht2
  
  # replace empty cell with 'unknown'
  class_names = gsub("^$","unknown",class_names)
  
  # create class2col index
  n = length(unique(class_names))
  color.hex = rainbow(n, s = 1, v = 1, start = 0, end = max(1, n - 1)/n,
                      alpha = NULL, rev = FALSE)
  comp.color.index = data.frame(unique(class_names),color.hex)
  
  # create table of all metabolites with their related pathway and col
  class2col = data.frame(class_names = class_names, class_col = class_names)
  i=1
  while (i <= nrow(comp.color.index)) {
    class2col[,2] = sub(pattern = paste0(".*",comp.color.index[,1][i],".*"),
                         replacement = comp.color.index[,2][i],
                         x = class2col[,2])
    i=i+1
  }
  
  ############################################
  # make only the chose pathway to be colored (all the others will get white color)
  if (is.null(comp_pathways) == F) {
    u=1
    while (u <= nrow(class2col)) {
      if (all(class2col$class_names[u] != comp_pathways) == T) {
        class2col$class_col[u] = "#ffffff" 
      }
      u=u+1
    }
  }
  
  
  
  ############
  # prepare right side class annotation
  #col_fun = circlize::colorRamp2(1:nrow(class2col), class2col$class_col)
  #xx <- class2col$class_col
  #yy <- class2col$class_names
  #zz =list()
  #for(i in 1:length(yy)) {zz[[yy[i]]]<-xx[i]}
  
  #ht3 = Heatmap(class2col$class_names, col = list(class = setNames(class2col$class_col,class2col$class_names)))
  
  #comp_colored_0  = grep(comp_pathways,class2col$class_names)
  #comp_colored = class2col$class_names[comp_colored_0]
  #ha2 = HeatmapAnnotation(class = row.names(cor_file1))
  #row_ha = rowAnnotation(. = class2col$class_names,
  #                       col = list(class = setNames(class2col$class_col,
  #                                                   class2col$class_names)))
  
  #############
  col_fun = list(class = setNames(class2col$class_col,
                                  class2col$class_names))
  
  ha = HeatmapAnnotation(class = class2col$class_names,
                         col = col_fun,
                         show_annotation_name = F)
  
  ht1 = Heatmap(cor_file1, name = "corralation",show_row_names = F, show_column_names = F,
                bottom_annotation = ha)#, right_annotation = row_ha)
  ht1# + ht3
}
