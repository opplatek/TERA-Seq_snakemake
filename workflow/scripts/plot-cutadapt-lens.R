#!/usr/bin/env Rscript
#
# Plot parsed Cutadapt trimmed adapter lengths per library
#
# Two arguments - input parsed Cutadapt tsv (you can use parse-cutadapt-lens.py), output pdf
#

# args <- commandArgs(trailingOnly = TRUE)

### TESTING VARIABLES ###
# args <- NULL
# args <- "data/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/log/cutadapt.len.tsv"
### TESTING VARIABLES ###

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("ggplot2"))

ifile <- snakemake@input[[1]] # args[1]
ofile <- snakemake@output[[1]] # args[2]

trimming <- read.table(ifile, sep = "\t", header = T) %>%
  select("library", "length", "count")

trimming.m <- reshape2::melt(trimming, id.vars = c("length", "library"), value.name = "count")

trimming.m <- trimming.m %>%
  group_by(variable) %>%
  mutate(freq = count / sum(count)) %>%
  group_by(library) %>%
  mutate(tot_trimmed = paste0(library, " (tot. trimmed: ", sum(count), ")")) %>%
  ungroup() %>%
  mutate(library = tot_trimmed)

max_val <- trimming.m %>%
  group_by(library) %>%
  top_n(1, freq) %>%
  select(length, library)

breaks <- c(seq(20, 80, by = 10), median(max_val$length)) #, 58) # generate break positions; for 5tera-long
labels <- as.character(breaks) # and labels

p <- ggplot(trimming.m, aes(x = length, y = freq, color = library)) +
  geom_line() +
  scale_x_continuous(limits = c(20, 80), breaks = breaks, labels = labels) +
  theme_classic() +
  theme(legend.position = "bottom") +
  xlab("Removed length") +
  ylab("Frequency") +
  ggtitle("Removed sequence length")

p <- p +
  geom_vline(xintercept = median(max_val$length), color = "orange") +
#  geom_vline(xintercept = 58, color = "grey") + # for 5tera-long
  theme(plot.title = element_text(hjust = 0.5))

pdf(ofile)
  print(p)
dev.off()
