# Script to read/generate database files:
# 1. NCBI taxonomy db
# 2. 16S rRNA copy number
# 3. pathogen info

library(taxonomizr) # converts taxid to genus & species
library(readxl)
library(writexl)

makeTaxonomyDB <- function()  {
  # Get the taxonomy database
  getNamesAndNodes(outDir = ".",
             url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz",
             fileNames = c("names.dmp", "nodes.dmp"))
  read.names.sql('names.dmp')
  read.nodes.sql('nodes.dmp')
  # produces nameNote.sqlite
}

get16Sdb <- function () {
  # get the 16S rRNA copy number here
  # file: rrnDB-5.6_pantaxa_stats_NCBI.tsv
  rrndb <- read.csv("rrnDB-5.6_pantaxa_stats_NCBI.tsv", sep="\t")
  return(rrndb)
}

getPathogenDB <- function() {
  # get pathogen data here
  # pathogensDB.xlsx
  # uses getId from taxonomizr
  redo_table = FALSE # change to TRUE if need to regenerate the genus taxids
  pathDB <- read_xlsx("PathogensDB.xlsx")
  if (redo_table) {
    genus_taxid <- getId(pathDB$Genus)
    pathDB <- cbind(genusTaxid, pathDB)
    write_xlsx(pathDB, "PathogensDB.xlsx")
  }
  return(pathDB)
}


