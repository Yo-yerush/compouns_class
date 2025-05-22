load_csv_op2 <- function(dataFile) {
  # ------------------------------
  # upload lcms data file
  # ------------------------------
  raw.df<- read.csv(file=dataFile[1])
  
  # fix 1st column name
  names(raw.df)[1] = gsub(".*\\.\\.", "", names(raw.df)[1])
  lcms_edit<-raw.df
  
  # ------------------------------
  # upload lcms inchi file
  # ------------------------------
  lcms_inchi_df<- read.csv(file=dataFile[2])
  
  inchi_pos = grep("InCh", names(raw.df))
  inchi_type = names(raw.df)[inchi_pos]
  
  return(list(lcms_edit = lcms_edit,
              lcms_inchi_df = lcms_inchi_df,
              inchi_type = inchi_type))
}


