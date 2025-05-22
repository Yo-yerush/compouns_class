merge_RA_lib <- function(comp_file,
                         lcms_edit,
                         PeakArea_columns_pattern,
                         inchi_type = part_1$inchi_type) {
  # ------------------------------
  # merge "comp_file" table with RA library
  # ------------------------------
  
  lib_file = read.csv("data/RA_comp_library.csv")[,1:2]
  lib_merge = merge.data.frame(lib_file, lcms_edit, by = "Name")
  message(paste0("Merged with RA library: found ",nrow(lib_merge)," compounds\n"))
  
  
  # create "RA library data frame" from "comp_file"
  lib_merge_edit = cbind(Comp_Class = lib_merge[,"Comp_Class"],
                         Name = lib_merge[,"Name"],
                         chebi_id = rep(NA, nrow(lib_merge)),
                         InChI_col = rep(NA, nrow(lib_merge)),
                         lib_merge[,3:ncol(lib_merge)])
  names(lib_merge_edit) = gsub("InChI_col", inchi_type, names(lib_merge_edit))
  
  merge_lib_comp = rbind(lib_merge_edit, comp_file)
  
  # find position of the last "norm_area" value to remove duplicate rows
  lib_mark = grep(PeakArea_columns_pattern, names(merge_lib_comp))
  lib_mark = lib_mark[length(lib_mark)]
  
  # merge "RA library data frame" with "comp_file"
  comp_file_final = merge_lib_comp[!duplicated(merge_lib_comp[,lib_mark]),]
  message(paste0("Merged with library: found ",nrow(comp_file_final)," compounds\n"))
  return(comp_file_final = comp_file_final)
}
