library(ontologyIndex)

obo_in  <- "data/chebi.obo"
rds_out <- "chebi_slim.rds"

chebi_slim <- get_OBO(
  obo_in,
  extract_tags = "minimal",vector
  propagate_relationships = c() # skip extra rels, keeps only is_a
)
chebi_slim$definition <- NULL
chebi_slim$synonym    <- NULL
chebi_slim$namespace  <- NULL
# keep $parents / $ancestors – they’re required by get_ancestors()

saveRDS(chebi_slim, rds_out, compress = "xz")

chebi_db <- readRDS("chebi_slim.rds")
