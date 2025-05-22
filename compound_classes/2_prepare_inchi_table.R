prepare_inchi_table <- function(lcms_inchi_df, lcms_edit, inchi_type) {
  # ------------------------------
  # prepare data frame
  # ------------------------------
  # upload ChEBI structures_lite file
  if (inchi_type == "InChI.Key") {
    inchi_df = read.csv("data/structures_lite_InChIKey.csv")
  } else {
    inchi_df = read.csv("data/structures_lite_inchi.csv")
    }

  names(inchi_df)[2] = inchi_type
  message("\nuploaded ChEBI structures_lite file")
  
  
  # merge LCMS inchi and data files created by "load_csv" function
  lcms_inchi_df_final = lcms_inchi_df[,c("Name",inchi_type)] %>% filter(!duplicated(Name))
  #names(lcms_inchi_df_final)[2] = "InChI"
  
  lcms_df = merge.data.frame(lcms_edit,lcms_inchi_df_final, by = "Name")
  message("extracted InChI from LCMS data file")
  
  # merge new data frame (contain InChI column) with CHEBI_id
  comp_file_0 = merge.data.frame(inchi_df, lcms_df, by = inchi_type)
  
  # order column positions: inchi and chebi_id (after name) and pathway_column (befor "Name")
  comp_file = comp_file_0[,-c(1,2)] %>% 
    mutate(comp_file_0[,c(2,1)], .after = Name) %>%
    mutate("Comp_Class" = rep(NA,nrow(comp_file_0)), .before = Name)
  
  #  message(paste0("InChI to ChEBI id: found ",nrow(comp_file)," compounds\n"))
  
  return(comp_file = comp_file)
}
    
