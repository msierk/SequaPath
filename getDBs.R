# Script to read/generate database files:
# 1. NCBI taxonomy db
# 2. 16S rRNA copy number
# 3. pathogen info

library(taxonomizr) # converts taxid to genus & species
library(readxl) # to read in PathodgensDB.xlsx
library(writexl) # to write out new PathogensDB.xlsx

makeTaxonomyDB <- function()  {
  # Get the taxonomy database
  getNamesAndNodes(outDir = ".",
             url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz",
             fileNames = c("db/names.dmp", "db/nodes.dmp"))
  read.names.sql('db/names.dmp')
  read.nodes.sql('db/nodes.dmp')
  # produces nameNode.sqlite
}

makeAcc2taxidDB <- function() {
  setwd("/Users/msierk/Dropbox/Consulting/Nephros/MinION/data")
  read.accession2taxid("nucl_gb.accession2taxid", "accession2taxid.sqlite")
}

get16Sdb <- function () {
  # get the 16S rRNA copy number here
  # file: rrnDB-5.6_pantaxa_stats_NCBI.tsv
  # from https://rrndb.umms.med.umich.edu/static/download/?C=N;O=D
  temp <- tempfile()
  download.file("https://rrndb.umms.med.umich.edu/static/download/rrnDB-5.6_pantaxa_stats_NCBI.tsv.zip",temp)
  rrndb <- read.csv(unz(temp, "rrnDB-5.6_pantaxa_stats_NCBI.tsv"), sep="\t")
  unlink(temp)
  return(rrndb)
}

getPathogenDB <- function(redo_table) {
  #print(paste0("redo_table: ", redo_table))
  # get pathogen data here
  # default file is pathogensDB.xlsx
  # uses getId from taxonomizr
  #redo_table = TRUE if need to regenerate the genus taxids
  pathDB = data.frame()
  if (redo_table) {
    print("Redoing the pathogens DB...")
    pathDB <- read_xlsx("db/SequaPathPathogenList20200923.xlsx")
    genusTaxid <- getId(pathDB$Genus)
    pathDB <- cbind(genusTaxid, pathDB)
    write_xlsx(pathDB, "db/PathogensDB.xlsx")
  } else {
    pathDB <- read_xlsx("db/PathogensDB.xlsx")
  }
  return(pathDB)
}


