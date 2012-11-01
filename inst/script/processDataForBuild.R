################################################################################
## This is just a set of scripts that will get IDs that match with entrez gene
## IDs and then populate them into a DB.

## Call this script from src_build.sh, and AFTER organism_annotation/script.sh

## 1) Get a set of entrez gene IDs for a species (and tax ID) (
## (extract ALL of this information from the relevant chipsrc DB)

library("RSQLite")
library("UniProt.ws")
##Connect to the DB
drv <- dbDriver("SQLite")

## Path to all the Db's.
dir <- file.path("/mnt/cpb_anno/mcarlson/proj/mcarlson/sqliteGen/annosrc/db")
## ## Use this one (temporarily while the /db rebuilds)
## dir <- file.path("/home/mcarlson/arch/x86_64/R-devel/library/human.db0/extdata")
    


## For now just get data for the ones that we have traditionally supported
## I don't even know if the other species are available...
speciesList = c("chipsrc_human.sqlite",
  "chipsrc_rat.sqlite",
  "chipsrc_chicken.sqlite",
  "chipsrc_zebrafish.sqlite",
  #  "chipsrc_worm.sqlite",
  #  "chipsrc_fly.sqlite",
  "chipsrc_mouse.sqlite",
  "chipsrc_bovine.sqlite",
  #  "chipsrc_arabidopsis.sqlite"  ## this is available and could be "activated"
  ## But to activate arabidopsis, remember you have to pre-add the tables...
  #  "chipsrc_canine.sqlite",
  #  "chipsrc_rhesus.sqlite",
  #  "chipsrc_chimp.sqlite",
  #  "chipsrc_anopheles.sqlite"
  )

## ## temp for testing:
## speciesList = c("chipsrc_human.sqlite")
## speciesList = c("chipsrc_arabidopsis.sqlite")


## 2) use the package accessors to get all the relevant data back (IPI, PFAM
## and Prosite IDs).

## 3) AND: merge all these frames together into two data frames:
## SO        gene_id, ipi_id, pfam_id
## & ALSO    gene_id, ipi_id, prosite_id


getuniProtAndIPIs <- function(genes){
  
  ## get the UniProt IDs.
  ## issue: URI is too large, so I need a helper since I will have to do this
  ## in pieces...  - and actually, helper mapUniprot should do that
  ## automagically.
  ups <- UniProt.ws:::mapUniprot(from='P_ENTREZGENEID',to='ACC',query=genes)
  ## get the IPI IDs
  upKeys <- as.character(t(ups$ACC))
  ips <- UniProt.ws:::mapUniprot(from='ACC',to='P_IPI', upKeys)
  
  ## return as a single frame.
  ## Currently, I use an inner join here b/c DB is gene centric 
  base <- merge(ups, ips, by.x ="ACC", by.y ="ACC")#, all=TRUE)
  base
}


getData <- function(dbFile, db){
  ## look up the tax ID
  taxId <- sqliteQuickSQL(db, "SELECT value FROM metadata WHERE name='TAXID'")
  ## look up the entrez gene IDs
  if(dbFile != "chipsrc_arabidopsis.sqlite"){ ## if its not arabidopsis
    genes <- sqliteQuickSQL(db, "SELECT gene_id FROM genes")
  }else{
    genes <- sqliteQuickSQL(db, "SELECT gene_id FROM entrez_genes")
  }
  genes <- as.character(t(genes))

  ## get the UniProt and IPI Id's (merged into a table)
  base <- getuniProtAndIPIs(genes)
  
  
  ## get the pfam Id's
  pfam <- UniProt.ws:::getOneToMany(taxId, type="PFAM")
  colnames(pfam) <- c("ACC", "PFAM")
  ## and the prosite Id's.
  prosite <- UniProt.ws:::getOneToMany(taxId, type="prosite")
  colnames(prosite) <- c("ACC", "PROSITE")
  
  
  
  ## merge it all together. # 3 above - then return a list of length two
  ## Currently I am using an inner join here b/c the DB is gene centric, so
  ## there is no benefit fo haveing pfam/UniProt accessions that are not
  ## connected to an entrez gene
  lst <- list()
  lst[[1]] <- merge(base, pfam, by.x="ACC", by.y="ACC") #,all=TRUE)
  lst[[2]] <- merge(base, prosite, by.x="ACC", by.y="ACC") #,all=TRUE)
  names(lst) <- c("pfam","prosite")

## finally, be sure to drop the UniProt IDS?  I think we should keep
## em...  ;) Later on I can make use of them to enhance the devel annots.
  lst
}





## 4) For each species, get the data, using getData, and then go
## straight to humansrc.sqlite etc and populate the pfam and prosite
## tables. ALSO, be sure to add entries to metadata about where the data came
## from.  (and remove relevant code from the scripts).  - See the bindb.R
## script in ensembl/script.


## Helper for doing inserts from pfam
doInserts <- function(db, table, data){

  ## make a temp pfam table pfamt
  sqlCreate <- paste0("CREATE TABLE ",table,"t (
                gene_id TEXT,
                ipi_id TEXT,
                ",table,"_id TEXT);")
  sqliteQuickSQL(db, sqlCreate)

  
  ## 1st insert for pfam
  sqlIns <- paste0("INSERT into ",table,"t
             (gene_id, ipi_id, ",table,"_id)
             VALUES ($P_ENTREZGENEID,$P_IPI,$",toupper(table),")")
  dbBeginTransaction(db)
  rset <- dbSendPreparedQuery(db, sqlIns, data)
  dbClearResult(rset)
  dbCommit(db)

  ## then insert into the real pfam table
  sqlIns2 <- paste0("INSERT INTO ",table,"
                SELECT DISTINCT g._id as _id, i.ipi_id, i.",table,"_id
                FROM genes as g, ",table,"t as i
                WHERE g.gene_id=i.gene_id
                ORDER BY _id")
  sqliteQuickSQL(db, sqlIns2)

  
  ## then drop the table
  sqlDrop <- paste0("DROP TABLE ",table,"t")
  sqliteQuickSQL(db, sqlDrop)
  
}


## So: loop, where we call getData and then just populate the tables
require("RSQLite")

for(species in speciesList){

  ## DB connection
  db <- dbConnect(drv,dbname=file.path(dir, species))
  
  message("Getting data for:",species)
  res <- getData(species, db)

  message("Inserting data for:",species)
  ## Now I need to insert the data:
  doInserts(db, "pfam", res[["pfam"]])
  doInserts(db, "prosite", res[["prosite"]])
  
  ## And then I need to add metadata:
  date <- date()
  url <- "http://www.UniProt.org/"
  name <- "Uniprot"
  sqlMeta1 <- paste0("INSERT INTO metadata (name,value) VALUES ('UPSOURCENAME','",name,"')")
  sqliteQuickSQL(db, sqlMeta1)
  sqlMeta2 <- paste0("INSERT INTO metadata (name,value) VALUES ('UPSOURCEURL','",url,"')")
  sqliteQuickSQL(db, sqlMeta2)
  sqlMeta3 <- paste0("INSERT INTO metadata (name,value) VALUES ('UPSOURCEDATE','",date,"')")
  sqliteQuickSQL(db, sqlMeta3)
  sqliteQuickSQL(db,"DELETE FROM metadata WHERE name LIKE 'IPISOURCE%'")

  
  ## And don't forget the map_counts for PROSITE AND PFAM
  sqlmapcnt1 <- "INSERT INTO map_counts
                 SELECT 'PFAM', count(DISTINCT _id)
                 FROM pfam;"
  sqliteQuickSQL(db, sqlmapcnt1)

  sqlmapcnt2 <- "INSERT INTO map_counts
                 SELECT 'PROSITE', count(DISTINCT _id)
                 FROM prosite;"
  sqliteQuickSQL(db, sqlmapcnt2)

  ## ALSO: modify the map_metadata (1st drop the PFAM and prosite entries
  sqliteQuickSQL(db, "DELETE FROM map_metadata where map_name ='PFAM' ") 
  sqliteQuickSQL(db, "DELETE FROM map_metadata where map_name ='PROSITE' ")
  ## then put our own entries in...
  sqlPFMM <- paste0( "INSERT INTO map_metadata (map_name, source_name, ",
                    "source_url, source_date) VALUES ('PFAM','",name,
                    "','",url,"','",date,"')")
  sqliteQuickSQL(db, sqlPFMM)
  sqlPSMM <- paste0( "INSERT INTO map_metadata (map_name, source_name, ",
                    "source_url, source_date) VALUES ('PROSITE','",name,
                    "','",url,"','",date,"')")
  sqliteQuickSQL(db, sqlPSMM)
  
}



## 5) ALSO: Be sure to also add metadata to each DB as we loop!


## 6) ALSO: Be sure to add map counts for PFAM and PROSITE too.
