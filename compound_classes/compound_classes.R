if (FALSE){
  setwd("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/chebi db/pipline_scripts/package_files_310523/")
  dataFile_1 = "data/test_data/test_data_op1.csv"
  dataFile_2 = c("data/test_data/test_data_op2.csv","data/test_data/test_inchi_op2.csv")
  dataFile_3 = "data/test_data/test_data_op3.csv"
  
  #dataFile
  path_to_save = NULL
  pathway_class = 4
  PeakArea_columns_pattern = "Norm..Area.."
  parent_class_file = "data/parent_class.csv"
  chebi_db_file_uploaded = NULL
  InChi_column_exist = F
  merge_with_RA_library = T

  #  structures_lite_file = "data/structures_lite_inchi.csv"
  #  chebi_db_file = "data/chebi.obo"
  #  RA_lib_file = "data/RA_comp_library_100423.csv"
  }

compound_classes <- function(dataFile,
                            path_to_save = NULL,
                            pathway_class = 4,
                            PeakArea_columns_pattern = "Norm..Area..",
                            parent_class_file = "data/parent_class.csv",
                            chebi_db_file_uploaded = NULL,
                            InChi_column_exist = F,
                            merge_with_RA_library = T) {
  suppressWarnings(library(dplyr))
  setwd("C:/Users/yonatany/Migal/Rachel Amir Team - General/yonatan/chebi db/pipline_scripts/package_files_310523/")
  
  ############ part 1
  ### option 3
  if (InChi_column_exist) {
    source("compound_classes/1_load_csv_op3.R")
    part_1 = load_csv_op3(dataFile)
  } else {
    if (length(dataFile) == 1) {
      ### option 1
      source("compound_classes/1_load_csv_op1.R")
      part_1 = load_csv_op1(dataFile, PeakArea_columns_pattern)
    } else if (length(dataFile) == 2) {
      ### option 2
      source("compound_classes/1_load_csv_op2.R")
      part_1 = load_csv_op2(dataFile)
    }
  }
  
  ############ part 2
  source("compound_classes/2_prepare_inchi_table.R")
  part_2 = prepare_inchi_table(part_1$lcms_inchi_df,
                               part_1$lcms_edit,
                               part_1$inchi_type)
  
  ############ part 3
  source("compound_classes/3_prepare_final_table.R")
  part_3 = prepare_final_table(part_2,
                               chebi_db_file_uploaded,
                               parent_class_file,
                               pathway_class)
  
  ############ part 4
  if (merge_with_RA_library == T) {
    source("compound_classes/4_merge_with_RA_library.R")
    part_4 = merge_RA_lib(part_3,
                          part_1$lcms_edit,
                          PeakArea_columns_pattern,
                          part_1$inchi_type)
  } else {part_4 = part_3}
  
  ############ save file
  if (is.null(path_to_save) == F) {
    write.csv(part_4, paste0(path_to_save,"compound_classes_",pathway_class,".csv"), row.names = F)
    message(paste0("\nfile saved:\n",path_to_save,"compound_classes_",pathway_class,".csv"))
  } else {return(part_4)}
}
