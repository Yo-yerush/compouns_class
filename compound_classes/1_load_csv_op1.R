load_csv_op1 <- function(dataFile, PeakArea_columns_pattern) {

  raw.df<- read.csv(file=dataFile)
  
  # ------------------------------
  # (1) fix 1st column name
  # ------------------------------
  names(raw.df)[1] = gsub(".*\\.\\.", "", names(raw.df)[1])
  
  # ------------------------------
  # (2) mark tables
  # ------------------------------
  
  #firstColumnsNamePos = grep("Name", names(raw.df))
  #firstColumnsName = names(raw.df)[firstColumnsNamePos]
  #raw.df[raw.df$Name == "",firstColumnsNamePos] = NA
  
  
  #PeakAreaColumnsName<-names(raw.df)[1]
  
  FirstPeakAreaColumnsName<-names(raw.df)[grep(PeakArea_columns_pattern, names(raw.df))][1]
  #raw.df[,1] = na_if(raw.df[,1], '')
  
  
  raw_df_1<-raw.df %>%
    mutate("table_type"=ifelse(is.na(get(FirstPeakAreaColumnsName)), 
                               "lcms_inchi_df","lcms_edit"), 
           .before=1)%>%
    mutate("boolean_mark_table"=(table_type=="lcms_edit"), .before=Name)%>%
    mutate("id"=cumsum(boolean_mark_table),.before=table_type)%>%
    select(-boolean_mark_table)
  # ------------------------------
  # (3) split tables
  # ------------------------------
  # blue lines : the metabolites statistics
  lcms_edit<-raw_df_1%>%
    filter(table_type=="lcms_edit")
  
  # orange lines : the inchi code
  lcms_inchi_df_0<-raw_df_1%>%
    filter(table_type=="lcms_inchi_df")
  
  
  # ------------------------------
  # (4) Managing lcms_inchi_df
  # ------------------------------
  # remove empty columns
  emptyColumns.ix<-which((is.na(lcms_inchi_df_0[1,])) | (lcms_inchi_df_0[1,]==""))
  lcms_inchi_df_1<-lcms_inchi_df_0[,-emptyColumns.ix]
  
  # fix lcms_inchi_df columns names
  names(lcms_inchi_df_1)<-lcms_inchi_df_1[1,]
  
  lcms_inchi_df_2<-lcms_inchi_df_1%>%
    filter(`Name` != "Name")
  
  lcms_edit = lcms_edit[,-c(1,2)]
  lcms_inchi_df = lcms_inchi_df_2[,-c(1,2)]
  
#  inchi_pos = grep("InCh", names(raw.df))
#  inchi_type = "InChI"
#  names(raw.df)[inchi_pos] = inchi_type
  
  inchi_pos = grep("InCh", names(lcms_inchi_df))
  inchi_type = "InChI"
  names(lcms_inchi_df)[inchi_pos] = inchi_type
  
  return(list(lcms_edit = lcms_edit,
              lcms_inchi_df = lcms_inchi_df,
              inchi_type = inchi_type))
}
