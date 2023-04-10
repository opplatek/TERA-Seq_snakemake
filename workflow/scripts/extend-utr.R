#!/usr/bin/env Rscript
#
# Extend first and last exon. If UTRs are present, extend them as well so they match the exons.
# Adjusts sizes of UTRs, exons, transcripts, and genes
# Very simple, primarily done for one of the yeast projects
#
# IMPORTANT: This might (and will in compact genomes) cause overlapping regions!
# IMPORTANT: As of 06/23/2022 the script can either copy and extend existings UTRs or make new ones. It cannot do both.
#
# TODO: Check existing UTRs and extend them. Make new UTRs only for transcripts with start & stop codon and without UTR in the annotation.
# TODO: Disallow overlapping of UTRs with other UTRs or exons - but how to decide what is correct UTR and what's not?
#
################################################################################

suppressPackageStartupMessages(library("dplyr"))
# suppressPackageStartupMessages(library("GenomicRanges")) # loaded by rtracklayer
suppressPackageStartupMessages(library("rtracklayer"))

################################################################################
# Add order column based on group and sort column (multiple features in a single group will get numerical value based on the sort column)
add_order_col <- function(tab, group_col, sort_col, revs = FALSE) {
  if (revs) {
    tab <- tab %>%
      as.data.frame() %>%
      group_by(!!sym(group_col)) %>%
      arrange(desc(!!sym(sort_col))) %>%
      mutate(exon_order = row_number())
  } else {
    tab <- tab %>%
      as.data.frame() %>%
      group_by(!!sym(group_col)) %>%
      arrange(!!sym(sort_col)) %>%
      mutate(exon_order = row_number())
  }

  return(tab)
}

# Function to extend region (strand specific) by xx bp upstream (5' UTR) or downstream (3' UTR); https://support.bioconductor.org/p/78652/
extend <- function(x, upstream = 0, downstream = 0) {
  if (any(strand(x) == "*")) {
    warning("'*' ranges were treated as '+'")
  }
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}

################################################################################

opt <- NULL
opt$gtf <- snakemake@input[["gtf"]] # Input GTF
opt$fiveutr <- snakemake@params[["extend_fiveutr"]] # Five UTR
opt$threeutr <- snakemake@params[["extend_threeutr"]] # Three UTR
opt$gen_ind <- snakemake@input[["fai"]] # Genome FAI index (chrom\tsize
opt$output <- snakemake@output[[1]] # Output GT

# args <- commandArgs(trailingOnly = TRUE)
# opt$gtf <- args[1] # Input GTF
# opt$fiveutr <- args[2] # Five UTR
# opt$threeutr <- args[3] # Three UTR
# opt$gen_ind <- args[4] # Genome FAI index (chrom_name\tsize)
# opt$output <- args[5] # Output GTF

################################################################################
# ### TESTING VARIABLES ###
# opt <- NULL
# opt$gtf <- "/usr/local/lib/R/site-library/Saccharomyces_cerevisiae.R64-1-1.109.gtf.gz"
# opt$gen_ind <- "/usr/local/lib/R/site-library/genome.fa.fai"
# opt$fiveutr <- 200
# opt$threeutr <- 150
# opt$output <- "/usr/local/lib/R/site-library/Saccharomyces_cerevisiae.R64-1-1.109.extend-utr.gtf"
# ### TESTING VARIABLES ###
################################################################################
if (opt$fiveutr == "NA" || opt$threeutr == "NA") {
  print(paste0("Snakemake params extend_fiveutr or extend_threeutr are NA. Only renaming the input GTF to the output GTF."))
  system(paste0("cat ", opt$gtf, " > ", opt$output))
} else {
  gtf <- rtracklayer::import(opt$gtf, format = "gtf")
  length(gtf)
  gtf$rows <- seq(1:length(gtf)) # add a helper column to copy back all the missing rows at the end
  gtf.bckp <- gtf

  gtf_gene <- gtf[gtf$type == "gene"] # Get genes lines
  gtf_transcript <- gtf[gtf$type == "transcript"] # Get transcript lines
  gtf <- gtf[!is.na(gtf$transcript_id) & gtf$type != "gene" & gtf$type != "transcript"] # Get all but genes, transcript, and transcript_id missing lines

  gen_ind <- read.table(opt$gen_ind, stringsAsFactors = F, header = F) # chrom length is V2
  gen_ind <- gen_ind[, c(1, 2)]
  colnames(gen_ind) <- c("seqnames", "chr_end")

  utr_exist <- FALSE
  if (sum(gtf$type == "three_prime_utr") > 0 & sum(gtf$type == "five_prime_utr") > 0) {
    utr_exist <- TRUE # Remember UTRs exist
    print("I see the annotation contains UTRs. I am going to extend them. I will not make new ones!")

    plus_start <- gtf[strand(gtf) == "+" & gtf$type == "five_prime_utr"]
    minus_start <- gtf[strand(gtf) == "-" & gtf$type == "five_prime_utr"]
    plus_stop <- gtf[strand(gtf) == "+" & gtf$type == "three_prime_utr"]
    minus_stop <- gtf[strand(gtf) == "-" & gtf$type == "three_prime_utr"]

    gtf <- gtf[gtf$type != "five_prime_utr" & gtf$type != "three_prime_utr"] # Remove UTRs from the annotation because we'll extend them and append them later

    # Add order but sort it so we always take "exon" 1 - we sort it reverse and get only row index
    plus_start_ind <- add_order_col(plus_start, "transcript_id", "start") %>%
      filter(exon_order == 1) %>%
      pull(rows)
    plus_stop_ind <- add_order_col(plus_stop, "transcript_id", "end", revs = TRUE) %>%
      filter(exon_order == 1) %>%
      pull(rows)
    minus_start_ind <- add_order_col(minus_start, "transcript_id", "end", revs = TRUE) %>%
      filter(exon_order == 1) %>%
      pull(rows)
    minus_stop_ind <- add_order_col(minus_stop, "transcript_id", "start") %>%
      filter(exon_order == 1) %>%
      pull(rows)

    utr5_plus <- extend(plus_start[plus_start$type == "five_prime_utr" & plus_start$rows %in% plus_start_ind, ], upstream = opt$fiveutr, downstream = 0)
    utr5_minus <- extend(minus_start[minus_start$type == "five_prime_utr" & minus_start$rows %in% minus_start_ind, ], upstream = opt$fiveutr, downstream = 0)
    utr5 <- c(utr5_plus, utr5_minus)

    utr3_plus <- extend(plus_stop[plus_stop$type == "three_prime_utr" & plus_stop$rows %in% plus_stop_ind, ], upstream = 0, downstream = opt$threeutr)
    utr3_minus <- extend(minus_stop[minus_stop$type == "three_prime_utr" & minus_stop$rows %in% minus_stop_ind, ], upstream = 0, downstream = opt$threeutr)
    utr3 <- c(utr3_plus, utr3_minus)
  } else {
    print("I see the annotation doesn't contain UTRs. I am going to make them.")

    plus_start <- gtf[strand(gtf) == "+" & gtf$type == "start_codon"]
    minus_start <- gtf[strand(gtf) == "-" & gtf$type == "start_codon"]
    plus_stop <- gtf[strand(gtf) == "+" & gtf$type == "stop_codon"]
    minus_stop <- gtf[strand(gtf) == "-" & gtf$type == "stop_codon"]

    # Fix gbkey if exists not to confuse us
    if (any(colnames(elementMetadata(gtf)) == "gbkey")) {
      plus_start$gbkey <- "mRNA"
      minus_start$gbkey <- "mRNA"
      plus_stop$gbkey <- "mRNA"
      minus_stop$gbkey <- "mRNA"
    }

    utr5_plus <- extend(plus_start[plus_start$type == "start_codon", ], upstream = opt$fiveutr, downstream = 0)
    end(utr5_plus) <- start(plus_start) - 1
    utr5_minus <- extend(minus_start[minus_start$type == "start_codon", ], upstream = opt$fiveutr, downstream = 0)
    start(utr5_minus) <- end(minus_start) + 1
    utr5 <- c(utr5_plus, utr5_minus)

    utr3_plus <- extend(plus_stop[plus_stop$type == "stop_codon", ], upstream = 0, downstream = opt$threeutr)
    start(utr3_plus) <- end(plus_stop) + 1
    utr3_minus <- extend(minus_stop[minus_stop$type == "stop_codon", ], upstream = 0, downstream = opt$threeutr)
    end(utr3_minus) <- start(minus_stop) - 1
    utr3 <- c(utr3_plus, utr3_minus)

    utr5$type <- "five_prime_utr"
    utr3$type <- "three_prime_utr"
  }

  utrs <- c(utr5, utr3)

  if (!utr_exist) {
    utrs$exon_number <- NA
    utrs$rows <- NA # it inherited row numbers from start/stop codon, we don't want that because we use it to get missing rows at the end
  }

  # Get transcript with annotated start/stop codong for exon extension
  cds <- gtf[gtf$transcript_id %in% c(plus_start$transcript_id, minus_start$transcript_id, plus_stop$transcript_id, minus_stop$transcript_id)]
  cds <- cds[cds$type == "exon"] # get only exons

  # 5 UTR on + & 3 UTR on - (first exons)
  # Extend first exon by opt$fiveutr and last by opt$threeutr
  cds[cds$exon_number == "1"] <- extend(cds[cds$exon_number == "1"], upstream = opt$fiveutr, downstream = 0)

  # 3 UTR on + & 5 UTR on - (last exons)
  mcols(cds) <- mcols(cds) %>%
    as.data.frame() %>%
    group_by(transcript_id) %>%
    mutate(exon_number_max = as.character(max(as.numeric(exon_number)), na.rm = T)) %>% # exon numbers are stored as char, for some reason
    ungroup()

  cds[cds$exon_number == cds$exon_number_max] <- extend(cds[cds$exon_number == cds$exon_number_max], upstream = 0, downstream = opt$threeutr)

  cds$exon_number_max <- NULL

  # Extend genes & transcripts to match the extended UTRs
  if (length(gtf_gene) > 0) {
    gtf_gene <- extend(gtf_gene, upstream = opt$fiveutr, downstream = opt$threeutr)
  }
  if (length(gtf_transcript) > 0) {
    gtf_transcript <- extend(gtf_transcript, upstream = opt$fiveutr, downstream = opt$threeutr)
  }

  # Merge all the parts and export
  out <- c(cds, utrs, gtf_gene, gtf_transcript)
  gtf_getback <- gtf.bckp[!(gtf.bckp$rows %in% out$rows)] # Put back all the rows we have left behind
  if (sum(gtf_getback$type %in% c("five_prime_utr", "three_prime_utr")) != 0) {
    print("We found five_prime_utr or three_prime_utr where we didn't expect them. Removing.")
    gtf_getback <- gtf_getback[!gtf_getback$type %in% c("five_prime_utr", "three_prime_utr")]
  }
  out <- c(out, gtf_getback) # Put back all the rows we have left behind
#  out <- c(out, gtf.bckp[!(gtf.bckp$rows %in% out$rows)]) # put back all the rows we have left behind
  print("If TRUE the UTRs were added correctly (by number of rows)")
  if (!utr_exist) {
    length(gtf_getback) + length(utrs) == length(out) # does original annotation + new utr rows equal to the total number of rows?
  } else {
    length(gtf_getback) == length(out) # does original annotation + new utr rows equal to the total number of rows?
  }
  out$rows <- NULL

  # Fix chromosome overflow
  start(out)[start(out) < 0] <- 1

  out$chr_end <- out %>%
    as.data.frame() %>%
    left_join(gen_ind) %>%
    pull(chr_end)
  end(out)[end(out) > out$chr_end] <- out$chr_end[end(out) > out$chr_end]
  out$chr_end <- NULL

  # Sort the result
  out <- sortSeqlevels(out)
  out <- sort(out, ignore.strand = TRUE)

  rtracklayer::export(out, opt$output, format = "gtf")
}
