## Script write by Anat Melamed 2020
## Script for extraction, quantification, abundance estimation and detwinning of
#   Brice Jagado & Magali Naville integration site data. 
## Input: R12 file in sidebyside sam format
# Launch: Rscript --vanilla 02_quantify_and_detwin.R <pathtoinput> <pathtooutputbasename>[optional] 
## Outline: 
#   1. read in file from each sample
#   2. reformat positions, sequences if on the reverse strand
#   3. collapse to list of shearsites. 
#   4. estimate abundance using <sonicLength>
#   5. detwin using sequence similarity
#   6. output list of sites for further analysis. 
## Output:
#   1. list of sites in <pathtooutputbasename>.sites
#   2. report in <pathtooutputbasename>.summary

# install Biostrings :
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  +     install.packages("BiocManager")
#  BiocManager::install("Biostrings")

#commandline on terminal : Rscript --vanilla chemin/02_quantify_and_detwin.R R12_sample.txt 2> err
# In loop : for filename in R12_*.txt 
# > do
# > Rscript --vanilla chemin/02_quantify_and_detwin.R $filename 2> $filename.err
# > done


## Setup ----

# load required libraries
library(Biostrings)
library(tidyverse) 
library(sonicLength)

# # Record package status
# sessionInfo()

# recover information from command line. 
commandline = commandArgs(trailingOnly = TRUE)

# setup paths to files. if output prefix is absent, it will be extracted from input name. 
inputpath <- commandline[1]
outputpath <- commandline[2]
# inputpath <- "~/Dropbox/Projects/collab-ext/jagado/temp/R12_42_bwa_SFV_co.txt"
# outputpath <- "~/Dropbox/Projects/collab-ext/jagado/temp/R12_42_bwa_SFV_co"
if(is.na(outputpath) | nchar(outputpath) == 0){
  outputpath <- gsub("\\..+", "", inputpath)
}

## Functions ----

whichort <- function(x){
  # test imput
  if (!is.vector(x)) {stop("input must be a vector")}
  
  # if single value, process once
  if(!as.character(intToBits(x)[5]) %in% c("01", "00")){stop("incorrect input")}
  res <- ifelse(as.character(intToBits(x)[5]) == "01", "-", "+")

  if(length(x) > 1){
    res <- character(length(x))
    for (i in 1:length(x)){
      res[i] <- whichort(x[i])
    }
  }
  
  return(res)
}

parse_cigar_neg <- function(cigar, seqlen){
  if (length(cigar) == 1) {
    # adjust seqlen based on deletions, insertions
    # inspired by script in http://www.bioinformatics.babraham.ac.uk/projects/bismark/deduplicate_bismark
    ops <- unlist(str_split(cigar, "\\d+")) # store the operations
    lens <- unlist(str_split(cigar, "\\D")) # store the length per operation
    ops <- ops[ops != ""]
    lens <- as.integer(lens[lens != ""]) # remove blanks from ends
    
    longcigar <- tibble(ops, lens)
    
    correctionfactor = seqlen + 
      sum(longcigar$lens[longcigar$ops == "D"]) - 
      sum(longcigar$lens[longcigar$ops == "I"]) - 
      parse_cigar_pos(cigar)
  } else{
    correctionfactor <- numeric(length(cigar))
    for (i in 1:length(cigar)) {
      correctionfactor[i] <- parse_cigar_neg(cigar[i], seqlen[i])
    }
  }
  return(correctionfactor)
}

parse_cigar_pos <- function(cigar){
  # function for calculating the correction factor required for negative stranded reads. 
  # NOTE: This is different from the fusion read version of this function
  # output: the correction factor (numeric) to use.
  if (length(cigar) == 1) {
    if (grepl("^\\d+M", cigar)) {
      correctionfactor <- 0
    } else if (grepl("^\\d+H\\d+M", cigar)) {
      correctionfactor <- 0
    } else if (grepl("^\\d+S\\d+M", cigar)) {
      correctionfactor <- as.numeric(gsub("^(\\d+)S\\d+M.*", "\\1", cigar))
    } else if (grepl("^\\d+I\\d+M", cigar)) {
      correctionfactor <- 0
    } else {
      stop("cigar not planned for ", cigar)
    }
  }
  
  if (length(cigar) > 1) {
    correctionfactor <- numeric(length(cigar))
    for (i in 1:length(cigar)) {
      correctionfactor[i] <- parse_cigar_pos(cigar[i])
    }
  }
  
  return(correctionfactor)
}

correct_reads <- function(x, nslices){
  x <- x %>% select(-starts_with("qual"), -identifier) %>% 
    group_by_all()%>% 
    tally() %>% 
    ungroup()
  
  x$ort1 <- whichort(x$flag1)
  x$ort2 <- whichort(x$flag2)
  x$seqlen1 <- nchar(x$samseq1)
  x$seqlen2 <- nchar(x$samseq2)
  x$cigarcorfact1 <- numeric(nrow(x))
  x$cigarcorfact1[x$ort1 == "-"] <- parse_cigar_neg(x$cigar1[x$ort1 == "-"], x$seqlen1[x$ort1 == "-"])
  x$cigarcorfact1[x$ort1 == "+"] <- parse_cigar_pos(x$cigar1[x$ort1 == "+"])
  x$cigarcorfact2 <- numeric(nrow(x))
  x$cigarcorfact2[x$ort2 == "-"] <- parse_cigar_neg(x$cigar2[x$ort2 == "-"], x$seqlen2[x$ort2 == "-"])
  x$cigarcorfact2[x$ort2 == "+"] <- parse_cigar_pos(x$cigar2[x$ort2 == "+"])
  x$seq1 <- character(nrow(x))
  x$seq1[x$ort1 == "-"] <- as.character(reverseComplement(DNAStringSet(x$samseq1[x$ort1 == "-"])))
  x$seq1[x$ort1 == "+"] <- x$samseq1[x$ort1 == "+"]
  x$pos1 <- numeric(nrow(x))
  x$pos1[x$ort1 == "-"] <- x$sampos1[x$ort1 == "-"]  + x$cigarcorfact1[x$ort1 == "-"]  - 1
  x$pos1[x$ort1 == "+"] <- x$sampos1[x$ort1 == "+"]  - x$cigarcorfact1[x$ort1 == "+"] 
  x$pos2 <- numeric(nrow(x))
  x$pos2[x$ort2 == "-"] <- x$sampos2[x$ort2 == "-"]  + x$cigarcorfact2[x$ort2 == "-"]  - 1
  x$pos2[x$ort2 == "+"] <- x$sampos2[x$ort2 == "+"]  - x$cigarcorfact2[x$ort2 == "+"] 
  x[, c("chr1", "pos1", "ort1", "seq1", "seqlen1", "chr2", "pos2", "ort2", "n")]
  
}


wrapestabun <- function(x) {
  safeestabun <- safely(estAbund)
  y <- safeestabun(x$cloneid, x$length, min.length = min(x$length)) 
  if (is.null(y$error)) {
    aft <- enframe(y$result$theta, name = "cloneid", value = "sisters")
  } else if (is.null(y$result)) {
    cat("\nsonicLength fail, shearsites kept as sisters, total clones", 
        length(unique(x$cloneid)), ", total shear sites", nrow(x), 
        ", total reads", sum(x$totaldupl), "\n\n") 
    aft <- tibble(cloneid = as.character(x$cloneid)) %>%
      count(cloneid) %>%
      rename(sisters = n)
  } else {
    stop("ERROR check wrapestabun")
  }
  return(aft)
}

cleansites <- function(sitestoclean, alltwins = NULL){
  cat("\n")
  
  for (i in 1:max(sitestoclean$rn)) {
    
    if (i %% 100 == 0) {cat(".")}
    if (i %% 10000 == 0) {cat(" ", as.character(Sys.time()), "\n")}
    if (!i %in% sitestoclean$rn[sitestoclean$keep]) {next}
    if (i >= max(sitestoclean$rn[sitestoclean$keep])) {break}
    
    if (is.null(alltwins)) {
      twins <- findtwins(sitestoclean, i) %>%
        mutate(decl = TRUE) %>%
        filter(decl)
    }else{
      twins <- alltwins[[i]]
      twins <- twins[twins$V2 %in% sitestoclean$rn[sitestoclean$keep], ]
    }

      removetwins <- 
      unique(c(twins$V1[twins$totaldupl.x == twins$totaldupl.y & 
                          twins$sisters.x == twins$sisters.y], 
               twins$V2))
    
    if (length(removetwins) != 0) {
      sitestoclean$keep[sitestoclean$rn %in% removetwins] <- FALSE
    }
  }
  
  return(sitestoclean)
}

findtwins <- function(sitestoclean, i){
  twins <- tibble(V1 = i, V2 = sitestoclean$rn[sitestoclean$rn > i] ) %>%
    inner_join(sitestoclean, by = c("V1" = "rn")) %>%
    inner_join(sitestoclean, by = c("V2" = "rn")) %>%
    filter(keep.x, keep.y) %>%
    mutate(subseq.x = str_sub(seqexample.x, 1, 20), 
           subseq.y = str_sub(seqexample.y, 1, 20)) %>%
    mutate(twin = agrepl(unique(subseq.x), subseq.y) | 
             (chr.x == chr.y & 
                pos.x <= pos.y + 2 & pos.x >= pos.y - 2 & 
                ort.x == ort.y)) %>%
    filter(twin) %>% 
    select(V1, V2, subseq.x, subseq.y, totaldupl.x, totaldupl.y, 
           sisters.x, sisters.y) 
  
  return(twins)
}

## Run starts here ----

# load read file 
timestamp()
cat("Reading in file: ", basename(inputpath), "...")
readdata <- read_tsv(file = file.path(inputpath), 
                     col_names = c("identifier", "flag1", "chr1", "sampos1", "mapq1", "cigar1", "samseq1", "qual1", 
                                   "flag2", "chr2", "sampos2", "mapq2", "cigar2", "samseq2", "qual2"), 
                     col_types = c("cicdic---ccicdic---cc-"))

# Test whether there was a problem with the reading of the quality string
if (length(which(
  (nchar(readdata$samseq1) != nchar(readdata$qual1)) |
  (nchar(readdata$samseq2) != nchar(readdata$qual2)))) > 0){
  stop("reads not read in correctly!\n"); quit(status = 123)
} else {
  cat("File read in, total", nrow(readdata), "reads. \n\n")
} 
  

# Correct integration site position, orientation, sequence. 
timestamp()
cat("Processing reads to unique pairs\n\n")
cleaneddata <- correct_reads(readdata)


# collapse to shear sites. 
shearsitedata <- cleaneddata %>%
  mutate(length = abs(pos2 - pos1)) %>%
  group_by(chr1, pos1, ort1, length) %>%
  arrange(chr1, pos1, ort1, length, desc(seqlen1)) %>%
  summarise(totalduplicates = sum(n), seqlen1 = first(seqlen1), seqexample = first(seq1)) %>%
  ungroup()


# Estimate clonal abundance using sonicLength
timestamp()
cat("Estimating clonal abundance\n\n")

intsitedata <- shearsitedata %>%
  group_by(chr1, pos1, ort1) %>%
  arrange(chr1, pos1, ort1, desc(seqlen1)) %>%
  summarise(totaldupl = sum(totalduplicates), shearsites = n(), 
            seqexample = first(seqexample)) %>%
  ungroup() %>%
  mutate(cloneid = paste(chr1, pos1, ort1, sep = "_"))

intsitedata <- shearsitedata %>%
  mutate(cloneid = paste(chr1, pos1, ort1, sep = "_")) %>%
  wrapestabun() %>%
  inner_join(intsitedata, by = "cloneid") %>%
  mutate(relabun = sisters / sum(sisters)) %>%
  arrange(desc(sisters), desc(totaldupl)) %>%
  select(cloneid, chr = chr1, pos = pos1, ort = ort1, totaldupl, shearsites, sisters, relabun, seqexample)
  

# Remove twins based on sequence similarity
timestamp()
cat("Clean twins from site list", "\n")

sitelist <- intsitedata %>%
  mutate(rn = row_number(), keep = TRUE) 

cleanedsites <- cleansites(sitelist) 

timestamp()
cat("Done,  exporting to output files:\n\t", paste0(outputpath, ".sites"), "\n\t", paste0(outputpath, ".summary"), "\n")

cleanedsites %>% 
  group_by(keep) %>% 
  summarise(file = basename(inputpath), 
            clones = n(), sisters = sum(sisters), totaldupl = sum(totaldupl)) %>%
  select(file, keep, clones, sisters, totaldupl) %>%
  write_tsv(paste0(outputpath, ".summary"))

cleanedsites %>%
  filter(keep) %>%
  select(-keep, -rn) %>%
  write_tsv(paste0(outputpath, ".sites"))

# cleanedsites %>%
#   select(-rn) %>%
#   write_tsv(paste0(outputpath, ".rawsites"))


## Session info when setting up script ----

comment(cleanedsites) <- "SessionInfo: 
3.5.0 (2018-04-23)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS  10.15.6

Matrix products: default
BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
  [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
  [1] splines   stats4    parallel  stats     graphics  grDevices utils    
[8] datasets  methods   base     

other attached packages:
  [1] sonicLength_1.4.6   forcats_0.4.0       stringr_1.4.0      
[4] dplyr_1.0.1         purrr_0.3.3         readr_1.3.1        
[7] tidyr_1.0.0         tibble_3.0.1        ggplot2_3.3.0      
[10] tidyverse_1.3.0     Biostrings_2.50.2   XVector_0.22.0     
[13] IRanges_2.16.0      S4Vectors_0.20.1    BiocGenerics_0.28.0

loaded via a namespace (and not attached):
  [1] Rcpp_1.0.3       cellranger_1.1.0 pillar_1.4.3     compiler_3.5.0  
[5] dbplyr_1.4.2     tools_3.5.0      zlibbioc_1.28.0  jsonlite_1.6    
[9] lubridate_1.7.4  lifecycle_0.2.0  gtable_0.3.0     pkgconfig_2.0.3 
[13] rlang_0.4.7      reprex_0.3.0     cli_2.0.1        rstudioapi_0.11 
[17] DBI_1.1.0        haven_2.2.0      withr_2.1.2      xml2_1.2.2      
[21] httr_1.4.1       fs_1.4.1         generics_0.0.2   vctrs_0.3.2     
[25] hms_0.5.3        grid_3.5.0       tidyselect_1.1.0 glue_1.4.0      
[29] R6_2.4.1         fansi_0.4.1      readxl_1.3.1     modelr_0.1.5    
[33] magrittr_1.5     backports_1.1.5  scales_1.1.0     ellipsis_0.3.0  
[37] rvest_0.3.5      assertthat_0.2.1 colorspace_1.4-1 stringi_1.4.5   
[41] munsell_0.5.0    broom_0.7.0.9001 crayon_1.3.4"
