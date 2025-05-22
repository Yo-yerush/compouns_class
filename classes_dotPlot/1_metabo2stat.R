load_csv_op1 <- function(pathway_table_for_metaboanalyst,
                         group1,
                         group2) {
  
  df = read.csv(pathway_table_for_metaboanalyst)
  names(df) = df[1,]
  names(df)[1] = "Name"
  df = df[-1,]
  
  pvalue_vec = NA
  ttest_vec = NA
  FC_vec = NA
  i=1
  while (i <= nrow(df)) {
    pvalue_vec[i] = t.test(as.numeric(df[i,grep(group1,names(df))]),
                           as.numeric(df[i,grep(group2,names(df))]))$p.value
    
    ttest_vec[i] = t.test(as.numeric(df[i,grep(group1,names(df))]),
                          as.numeric(df[i,grep(group2,names(df))]))$statistic
    
    FC_vec[i] = as.numeric(df[i,grep(group1,names(df))]) /
      as.numeric(df[i,grep(group2,names(df))])
    
    i=i+1
  }
  
  
  
  
  df_2.view = data.frame(Name = df[,"Name"],
                         Fold_Change = FC_vec,
                         p.value = pvalue_vec,
                         t = ttest_vec)
  
  return(df_2.view = df_2.view)
}