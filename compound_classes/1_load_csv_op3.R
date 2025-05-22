load_csv_op3 <- function(dataFile) {
  # ------------------------------
  # upload lcms data file
  # ------------------------------
  raw.df<- read.csv(file=dataFile)
  
  # fix 1st column name
  names(raw.df)[1] = gsub(".*\\.\\.", "", names(raw.df)[1])
  
  inchi_pos = grep("InCh", names(raw.df))
  inchi_type = names(raw.df)[inchi_pos]

  
  lcms_edit<-raw.df
  # remove inchi column from "lcms_edit" for "4_merge_with_RA_library" function
  lcms_edit = lcms_edit[,-grep(inchi_type, names(lcms_edit))]
  
  # ------------------------------
  # split to prepare inchi file
  # ------------------------------
  lcms_inchi_df<-raw.df[,c("Name",inchi_type)]
  
  return(list(lcms_edit = lcms_edit,
              lcms_inchi_df = lcms_inchi_df,
              inchi_type = inchi_type))
}


