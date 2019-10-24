################### Wilderlab eDNA sequence analysis pipeline  #######################
######################### compiled by Shaun Wilkinson ################################

## set working directory
setwd(dirname(file.choose()))


## install and load required packages
if(!("insect" %in% list.files(.libPaths()))) devtools::install_github("shaunpwilkinson/insect")
if(!("vegan" %in% list.files(.libPaths()))) install.packages("vegan")
if(!("BiocManager" %in% list.files(.libPaths()))) install.packages("BiocManager")
if(!("dada2" %in% list.files(.libPaths()))) BiocManager::install("dada2")
if(!("lulu" %in% list.files(.libPaths()))) devtools::install_github("tobiasgf/lulu")
library(insect)
library(vegan)
library(dada2)



## download example source files from
## https://www.dropbox.com/sh/0gyxeg44spprz2e/AACXq3P7BJkQP2bQMQlX5zupa?dl=1
## and unzip to working directory


## alternatively download and read in the example FASTQ files and index table using R
u1 <- "https://www.dropbox.com/s/hoxx6jr7j0qb3w7/MSRun306-FTP244_S1_L001_R1_001.fastq.gz?dl=1"
u2 <- "https://www.dropbox.com/s/3k2l53ydpuzzxcd/MSRun306-FTP244_S1_L001_R2_001.fastq.gz?dl=1"
u3 <- "https://www.dropbox.com/s/ph8dcrda4o54bnx/indices.csv?dl=1"
download.file(u1, destfile = "MSRun306-FTP244_S1_L001_R1_001.fastq.gz", mode = "wb")
download.file(u1, destfile = "MSRun306-FTP244_S1_L001_R2_001.fastq.gz", mode = "wb")
download.file(u3, destfile = "indices.csv", mode = "wb")


## read in three-column csv file with header specifying sample names, forward tags and reverse tags
indices <- read.csv("indices.csv", stringsAsFactors = FALSE, header = TRUE)
tags <- paste0(indices[, 2], ":", indices[, 3])
names(tags) <- indices[, 1]


## demultiplex and export trimmed fastq files to new folder called "demux"
demultiplex(R1 = "MSRun306-FTP244_S1_L001_R1_001.fastq.gz",
            R2 ="MSRun306-FTP244_S1_L001_R2_001.fastq.gz",
            tags = tags, up = "GARTCTTTGAACGCAAATGGC",
            down = "GCTTATTAATATGCTTAAATTCAGCG", destdir = "demux")


################ filter and merge with DADA2 pipeline #########################

path <- "demux"
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1:2)
# note original dada script produced duplicates
sample.names <- apply(sample.names, 2, paste0, collapse = "_")
plotQualityProfile(fnFs[1:16])
plotQualityProfile(fnRs[1:16])
filtFs <- file.path("filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("filtered", paste0(sample.names, "_R_filt.fastq.gz"))
# On Windows set multithread=FALSE
## fwds are 231bp, revs are 226bp
## trim to 230 & 224 = 454bp
## for Pocillopora we need 441bp so need to overlap 14 bp
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,224),
                     maxN=0, maxEE = c(2,5), truncQ = 2, rm.phix=TRUE,
                     compress = TRUE, multithread = TRUE)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool = TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool = TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
rownames(seqtab.nochim) <- sub("_R1$", "", rownames(seqtab.nochim))
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN),
               sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls:
# e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track) # gives you how many were dropped at each stage
write.csv(track, file = "trackprocess.csv", row.names = FALSE)


## remove positive control
seqtab.nochim <- seqtab.nochim[-(match("PosiC", rownames(seqtab.nochim))),]
colsums <- apply(seqtab.nochim, 2, sum)
discards <- colsums < 2
seqtab.nochim <- seqtab.nochim[, !discards]
saveRDS(seqtab.nochim, file = "seqtabnochim.rds")
write.csv(seqtab.nochim, file = "seqtabnochim.csv", row.names = FALSE)



################### taxon assignment with INSECT ############################

## extract sequences from seqtab.nochim table
x <- char2dna(colnames(seqtab.nochim))
names(x) <- paste0("ASV", seq_along(x))


## taxonomic classification with 'insect'
## full list of classifiers from https://github.com/shaunpwilkinson/insect
download.file("https://www.dropbox.com/s/f07cka6308ebk2o/classifier.rds?dl=1",
              destfile = "classifier.rds", mode = "wb") ## ITS2 classifier
classifier <- readRDS("classifier.rds")


#classifier <- readRDS("~/Dropbox/R/insect/trees/ITS2/cnidarian_ITS2_marine/version_5/classifier.rds")
longDF <- classify(x, classifier, mincount = 3, threshold = 0.8, decay = FALSE, cores = 4)
longDF <- cbind(longDF, t(seqtab.nochim))
write.csv(longDF, file = "longDF.csv", row.names = FALSE)


## aggregate table so only one row for each unique taxon
taxa <- aggregate(longDF[3:12], longDF["taxID"], head, 1)
counts <- aggregate(longDF[13:ncol(longDF)], longDF["taxID"], sum)
shortDF <- merge(taxa, counts, by = "taxID")
write.csv(shortDF, file = "shortDF.csv", row.names = FALSE)



#################### taxon-independent analysis ###########################

## install VSEARCH executable locally
## note: change 'linux' to 'macos' if working on a mac
system("wget https://github.com/torognes/vsearch/releases/download/v2.10.2/vsearch-2.10.2-linux-x86_64.tar.gz")
system("tar xzf vsearch-2.10.2-linux-x86_64.tar.gz")


## OTU clustering with VSEARCH
tf1 <- tempfile(fileext = ".fa")
tf2 <- tempfile(fileext = ".aln")
writeFASTA(x, tf1)
system(paste0("./vsearch-2.10.2-linux-x86_64/bin/vsearch --cluster_smallmem ", tf1, " --alnout ", tf2, " --id 0.97 --usersort"))
otus <- scan(tf2, what = "", sep = "\n", quiet = TRUE)
queries <- otus[grepl("^Query", otus)]
queries <- gsub("^Query >", "", queries)
targets <- otus[grepl("^Target", otus)]
targets <- gsub("^.+>", "", targets)
centers <- unique(targets)
queries <- c(queries, centers)
targets <- c(targets, centers)
singles <- names(x)[!names(x) %in% queries]
queries <- c(queries, singles)
targets <- c(targets, singles)
clusters <- insect:::.point(targets)
names(clusters) <- queries
clusters <- clusters[match(names(x), names(clusters))]
centrals <- names(clusters) %in% unique(targets)
OTUs <- t(seqtab.nochim)
rownames(OTUs) <- paste0("OTU", clusters)
OTUs <- aggregate(OTUs, by = list(rownames(OTUs)), sum)
rownames(OTUs) <- OTUs[, 1]
OTUs <- OTUs[, -1]
OTUs <- OTUs[order(as.integer(gsub("OTU", "", rownames(OTUs)))), ]
write.csv(OTUs, file = "OTUs.csv", row.names = FALSE)


## write sequences to fasta file for VSEARCH neighbor analysis
lulufa <- x[centrals]
names(lulufa) <- paste0("OTU", clusters[centrals])
lulufa <- lulufa[order(as.integer(gsub("OTU", "", names(lulufa))))]
tf1 <- tempfile(fileext = ".fa")
tf2 <- tempfile(fileext = ".txt")
writeFASTA(lulufa, tf1)
system(paste0("./vsearch-2.10.2-linux-x86_64/bin/vsearch --usearch_global ",
              tf1, " --db ", tf1, " --self --id .84 --iddef 1 --userout ", tf2,
              " -userfields query+target+id --maxaccepts 0 --query_cov .9 --maxhits 10"))
matchlist <- read.delim(tf2, header = FALSE, stringsAsFactors = FALSE)


## merge overextended OTUs with LULU
curated_result <- lulu::lulu(OTUs, matchlist)
lulutab <- curated_result$curated_table
lulutab <- lulutab[order(as.integer(gsub("OTU", "", rownames(lulutab)))), ]
saveRDS(lulutab, file = "lulutab.rds")
write.csv(lulutab, file = "lulutab.csv", row.names = FALSE)
discards <- names(lulufa) %in% curated_result$discarded_otus
lulufa <- lulufa[!discards]
writeFASTA(lulufa, file = "lulu.fa")


## heirarchical clustering of sites
Y <- t(lulutab)
Y <- Y[apply(Y, 1, sum) > 0, ] # remove sites with zero seqs, now 92 x 1886
distmat <- vegan::vegdist(Y, method = "bray", binary = TRUE)
dnd <- as.dendrogram(hclust(distmat, method = "average"))
plot(dnd)


## then do alpha, beta diversity indices, pca plots etc
H <- vegan::diversity(Y) # ...

#######################################
