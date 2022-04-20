Code used for processing raw microbiome fastq files to generate an ASV file for (Raina et al. 2022) Chemotaxis shapes the microscale organisation of the ocean's microbiome. This DADA2 pipeline is based on Callahan, B.J. et al. 2016. In this pipeline R1 and R2 were trimmed to remove low quality terminal ends (trunc(R1= 260; R2= 255)), this was found to be the trunc lengths that produced the most quality merged reads for this dataset.

#librarys for this work must be installed and loaded in R
library(R.utils);
library(dada2);
library(ShortRead); packageVersion("ShortRead") # dada2 depends on this
library(tidyverse); packageVersion("dplyr") # for manipulating data
library(Biostrings);  # for creating the final graph at the end of the pipeline
library(Hmisc); packageVersion("Hmisc") # for creating the final graph at the end of the pipeline
library(plotly); packageVersion("plotly") # enables creation of interactive graphs, especially helpful for quality plots
library(here); packageVersion("here");
here();
library(readr);
library(Biostrings);
library(DECIPHER); packageVersion("DECIPHER");

#STEP1: Set Directories and make subfolders for your files
home.dir <- ("/home/anna/jb_isca_2022/"); #BPA dataset location

setwd(home.dir);
plates<-read_csv('plates');

args = commandArgs(trailingOnly=TRUE);
p=args[1];
p=1
print(p);
plates[p,1];

home.dir <- ("/home/anna/jb_isca_2022/");

base.dir <- (paste0(home.dir, plates[p,1], "/"));
dir.create(paste0(plates[p,1], ".2021_exports_260_255"));
export_dir <-(paste0(plates[p,1], ".2021_exports_260_255/"));
dir.create(paste0(export_dir, "final"));
trimmed_dir <-(paste0(plates[p,1], ".trimmed/"));
dir.create(paste0(trimmed_dir));
trunc_dir <-(paste0(export_dir, "final"));
trimLeng <-(paste0(trimmed_dir, "final"));
dir.create(paste0(trimLeng));

files.fp.gz<-list.files(path=base.dir,  pattern=".fastq.gz");
data.fp.gz <- paste0(plates[p,1], "/");
print(data.fp.gz);

for (i in seq_along (files.fp.gz)){
  gunzip(filename=paste0(data.fp.gz,files.fp.gz[i]), overwrite=T)};

files.fp.gz<-list.files(path=base.dir,  pattern=".fastq");
data.fp.gz <- paste0(plates[p,1], "/");
print(data.fp.gz);

fnFs <- sort(list.files(data.fp.gz, pattern="_R1.fastq", full.names = TRUE));
fnRs <- sort(list.files(data.fp.gz, pattern="_R2.fastq", full.names = TRUE));

sample.names <- sapply(strsplit(basename(fnFs), "-"), `[`, 1);

fwd.plot.quals <- plotQualityProfile(fnFs[1:6])
ggsave(file = paste0(export_dir,"fwd.qualplot.pdf"), fwd.plot.quals, device="pdf");

rev.plot.quals <- plotQualityProfile(fnRs[1:6])
ggsave(file = paste0(export_dir,"rev.qualplot.pdf"), rev.plot.quals, device="pdf");

FWD <- "AGAGTTTGATCMTGGCTCAG";  ## CHANGE: Bacteria 16S: 27f V1_V3
REV <- "GWATTACCGCGGCKGCTG"; ## CHANGE: Bacteria 16S: 519r V1_V3

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD);
REV.orients <- allOrients(REV);

data.fp <- trimmed_dir;

fnFs.filtN <- file.path(data.fp, "filtN", basename(fnFs)); # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(data.fp, "filtN", basename(fnRs));

filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE, compress=F, matchIDs=TRUE);

passed.filtN <- file.exists(fnFs.filtN) # TRUE/FALSE vector of which samples passed the filter
fnFs.filtN <- fnFs.filtN[passed.filtN] # Keep only those samples that passed the filter
fnFs <- fnFs[passed.filtN] # Keep only those samples that passed the filter
fnRs.filtN <- fnRs.filtN[passed.filtN] # Keep only those samples that passed the filter
fnRs <- fnRs[passed.filtN] # Keep only those samples that passed the filter

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0));
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]));

#cutadapt was used to remove FWD and REV primers
cutadapt <- "/usr/bin/cutadapt";

data.fp <- paste0(trimmed_dir, "filtN");

#system2(cutadapt, args = "--version"); # Run shell commands from R

path.cut <- file.path(trimmed_dir, "cutadapt");
if(!dir.exists(path.cut)) dir.create(path.cut);
fnFs.cut <- file.path(path.cut, basename(fnFs));
fnRs.cut <- file.path(path.cut, basename(fnRs));

FWD.RC <- dada2:::rc(FWD);
REV.RC <- dada2:::rc(REV);
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC);
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC);

# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files Default error rate: 0.1
}

data.fp <- paste0(trimmed_dir, "cutadapt");

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0));
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]));

data.fp <- paste0(trimmed_dir, "cutadapt");

passed.fnFs.cut <- file.exists(fnFs.cut) # TRUE/FALSE vector of which samples passed the filter
fnFs.cut <- fnFs.cut[passed.fnFs.cut] # Keep only those samples that passed the filter
fnFs <- fnFs[passed.fnFs.cut] # Keep only those samples that passed the filter
fnRs.cut <- fnRs.cut[passed.fnFs.cut] # Keep only those samples that passed the filter
fnRs <- fnRs[passed.fnFs.cut] # Keep only those samples that passed the filter


sample.names <- sapply(strsplit(basename(fnFs), "-"), `[`, 1);

#STEP5: Dada2 trim step (CHANGE TRIM lengths in this step!)
#here you are renaming your new files to be sample name_R1_trim.fastq, no need to change
trimFs <- file.path(trimmed_dir, paste0(sample.names, "-R1_trim.fastq"));
trimRs <- file.path(trimmed_dir, paste0(sample.names, "-R2_trim.fastq"));

names(trimFs) <- sample.names;
names(trimRs) <- sample.names;
head(sample.names);


#In this pipeline R1 and R2 were trimmed to remove low quality terminal ends (trunc(R1= 260; R2= 255)), this was #found to be the trunc lengths that produced the most quality merged reads for this dataset.

out <- filterAndTrim(fnFs.cut, trimFs, fnRs.cut, trimRs, truncLen=c(260,255),
                     maxN=0, maxEE=c(4,6), truncQ=6, rm.phix=TRUE,
                     compress=FALSE, multithread=10, minLen = 50, matchIDs=TRUE); # On Windows set multithread=FALSE
head(out);

passed.trim <- file.exists(trimFs); # TRUE/FALSE vector of which samples passed the filter
trimFs <- trimFs[passed.trim]; # Keep only those samples that passed the filter
trimRs <- trimRs[passed.trim]; # Keep only those samples that passed the filter

data.fp <- paste0(trimLeng);

errF <- learnErrors(trimFs, nbases =1e8, verbose=TRUE, multithread=10, MAX_CONSIST=20);
errR <- learnErrors(trimRs, nbases =1e8, verbose=TRUE, multithread=10, MAX_CONSIST=20);

fwd.plot.errors <- plotErrors(errF, nominalQ=TRUE);
ggsave(file = paste0(trunc_dir,"/","fwd.errors.pdf"), fwd.plot.errors, width = 10, height = 10, device="pdf");
rev.plot.errors <- plotErrors(errR, nominalQ=TRUE);
ggsave(file = paste0(trunc_dir,"/","rev.errors.pdf"), rev.plot.errors, width = 10, height = 10, device="pdf");

#STEP7:  DEREPLICATION, DADA2 Step, MERGE, and COLLAPSE if the same! SAVE RDS
derepFs<-derepFastq(trimFs);
derepRs<-derepFastq(trimRs);

dadaFs <- dada(derepFs, err=errF, multithread=10, pool="pseudo"); #"pseudo" # or for pooled go: TRUE no ""
dadaRs <- dada(derepRs, err=errR, multithread=10, pool="pseudo"); 
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 12, maxMismatch = 1, verbose=TRUE);

# Inspect the merger data.frame from the first sample
head(mergers[[1]]);
seqtab <- makeSequenceTable(mergers);
dim(seqtab);
saveRDS(seqtab, file = paste0(trunc_dir,"/",plates[p,1],"seqtab_260_255_final_mm1.RDS"), ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)

#chimera removal

getN <- function(x) sum(getUniques(x));

table(nchar(getSequences(seqtab)));
seqtab.nochim.1 <- removeBimeraDenovo(seqtab, method="consensus", multithread=10, verbose=TRUE, minFoldParentOverAbundance=1); #consensus #pooled
dim(seqtab.nochim.1);
sum(seqtab.nochim.1)/sum(seqtab);

#collapse identical sequences
seqtab.nochim.1 <- collapseNoMismatch(seqtab.nochim.1, minOverlap = 200, verbose = T)

#save RDS files of the seqtab files of importance:
saveRDS(seqtab.nochim.1, file = paste0(trunc_dir,"/",plates[p,1],"seqtab.1_260_255_final_mm1.RDS"), ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)

#seqtab.1_260_255_final_mm1.RDS is used for assigning taxonomy with silva and for subsequent data analysis.