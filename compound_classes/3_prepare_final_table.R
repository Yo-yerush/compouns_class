prepare_final_table <- function(comp_file,
                                chebi_db_file_uploaded,
                                parent_class_file,
                                pathway_class) {
  # ------------------------------
  # upload CHEBI ontology data base (obo file)
  # ------------------------------
  
  suppressWarnings(library(ontologyIndex))
  #  suppressWarnings(library(stringr))
  
  if (is.null(chebi_db_file_uploaded) == T) {
    message("upload ChEBI ontology database file...")
    chebi_db = get_OBO("data/chebi.obo")
    message("uploaded ChEBI ontology database")
  } else {chebi_db = chebi_db_file_uploaded}
  
  # prepare ChEBI names and ids table
  chebi_names = data.frame(Name = chebi_db$name)
  chebi_names$chebi_id = row.names(chebi_names)
  
  # ------------------------------
  # upload parents table
  # ------------------------------
  parent_file = read.csv(parent_class_file)
  class_1 = parent_file[parent_file$class == "class_1", 1]
  class_2 = parent_file[parent_file$class == "class_2", 1]
  class_3 = parent_file[parent_file$class == "class_3", 1]
  class_4 = parent_file[parent_file$class == "class_4", 1]
  
  if (pathway_class == 4) {
    all_groups = c(class_1,class_2,class_3,class_4)
  } else if (pathway_class == 3) {
    all_groups = c(class_1,class_2,class_3)
  } else if (pathway_class == 2) {
    all_groups = c(class_1,class_2)
  } else if (pathway_class == 1) {
    all_groups = c(class_1)
  }
  
  parents_df = data.frame(Name = NA, chebi_id = NA)
  i=1
  while (i <= length(all_groups)) {
    parents_df[i,] = chebi_names[chebi_names$Name == all_groups[i],]
    i=i+1
  }
  
  # ------------------------------
  # prepare final table
  # ------------------------------
  i=1
  while (i <= nrow(comp_file)) {
    comp_id = comp_file[i,"chebi_id"]
    ances = get_ancestors(chebi_db, comp_id)
    
    u=1
    while (u <= nrow(parents_df)) {
      parent = ances[grep(parents_df[u,2], ances)]
      if (length(parent) != 0) {
        comp_file[i,"Comp_Class"] = parents_df[u,1]
      } 
      u=u+1
    }
    i=i+1
  }
  comp_file$Comp_Class[is.na(comp_file$Comp_Class)] = ""
  return(comp_file = comp_file)
}
 