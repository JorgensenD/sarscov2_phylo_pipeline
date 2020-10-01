## edit gisaid MSA
setwd("~/saudi_sarscov2")
require(ape)
require(seqinr)
#BiocManager::install('Biostrings')
require(Biostrings)
require(lubridate)
require(insect)
require(tidyverse)

## trimmed to posn 382 in current alignment previously -- end 30154(inc. gaps) prior to import


# use biostrings as crashes with read.DNA
full_alignment <- readDNAStringSet("./msa_0820/msa_0820_trim.fasta")
algn <- as.DNAbin(full_alignment)
#rm(full_alignment)

# remove region from names
names(algn) <- sub("\\|[^|]+$", "", names(algn))

#some sequences not full length
lengths <- sapply(algn,length)
algn2 <- algn[unname(which(lengths==29837))]

# decimal dates from these names
dec_dates <- sapply( strsplit( names(algn2), '\\|' ) , function(x){
  decimal_date( ymd( tail(x,1)))})



# rm sequences where the date is ambiguous/not complete
algn3 <- subset.DNAbin(algn2,!is.na(dec_dates))

## write this to file for distance matrix
#writeFASTA(algn3, file='algn3.fas')


# adapt metadata
## find duplicate sequences

## order by date 
new_dec_dates <- sapply( strsplit( names(algn3), '\\|' ) , function(x){
  decimal_date( ymd( tail(x,1)))})
ordin <- algn3
ordin <- ordin[order(new_dec_dates)]

## keep oldest of each duplicate sequence
subs <- unique.DNAbin(ordin)

meta <- readr::read_tsv('metadata.tsv')

#add seqnames to meta
meta$seqnames <- paste0('hCoV-19/', meta$strain, '|', paste(meta$gisaid_epi_isl, meta$date, sep='|') )

# -> innodups column in meta (mark these unique columns 1)
meta$inNoDups<- as.numeric(meta$seqnames %in% names(subs))
meta$sampleDate <- as.Date(meta$date) 

meta <- meta %>% 
  rename(
    Continent = region,
    Country = country,
    RegionOrState = division,
    CityOrCounty = location,
    seqName = seqnames
  )

write.csv(meta, 'meta.csv')

## running the tn93 etc. is very slow on a file this large w/ this many sequences. Would probably be ideal to downsample the inernational
## datta ahead of time in future




  