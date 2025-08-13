############### Load required packages ##############################

library(GenomicRanges)
library(naturalsort)
library(Biostrings)
library(data.table)
library(itertools)
library(fastmatch)
library(openxlsx)
library(foreach)
library(stringr)
library(seqinr)
library(tidyr)
library(dplyr)
library(vroom)
library(doMC)

################ Set global parameters ##############################

wd <- "/gpfs/ycga/work/brash/kk967/demo"  # demo working directory
if(getwd() != wd) setwd(wd)

hg38 <- "hg38/"
fastq <- "fastq/"
psdir <- hg38

blatopt <- "-ooc=11.ooc -fastMap -mask=lower -noHead -minMatch=8"
pslcom <- "pslPretty %s %s%s:%d-%d.fa %s -long -axt %s"
blatcom <- "blat %s%s:%d-%d.fa %s %s %s"

ct <- "iiiiiiiicciiiciiiiccc"
seqpat <- "'^[NACGT]+$'"

Sys.setenv("VROOM_SHOW_PROGRESS" = 0)
Sys.setenv("VROOM_THREADS" = 1)

if(!exists("probe.info.gr")) load("annot/probe.info.RData")
if(!exists("CA.info")) load("annot/CA.info.RData")

if(!exists("probes")) {
  load("annot/probes.RData")
  probes <- probes[['PS5']]
} 

flank = 250
mindist = 5

################ Set up parallel processing #########################

n.chunks <- 24
n.workers <- 24
registerDoMC(n.workers)

################ Run TruAlign pipeline on sample S ##################

mkTruAlign <- function(S, fq=fastq, clean=TRUE, split=TRUE, fa=TRUE, full=TRUE) {
  Sdir <- paste0(fq, S)
  if(fa) {
    mkfa(Sdir, clean, split)
  }
  if(full) {
    mkfullpanel(Sdir, clean)
  }
}

###################################

mkfa <- function(Sdir, clean, split, top=200) {
  fastqdir <- paste0(Sdir, "/Unaligned")
  chunkdir <- paste0(Sdir, "/chunks")
  if(!dir.exists(chunkdir)) dir.create(chunkdir)

  fqfls <- sprintf("%s/*_R%d_*fastq.??", fastqdir, 1:2)
  bcfls <- sprintf("%s/barcodes_R%d.txt", Sdir, 1:2)
  sqfls <- sprintf("%s/sequence_R%d.txt", Sdir, 1:2)

  if(clean) {
    print("Cleaning chunkdir..."); flush.console()
    file.remove(list.files(chunkdir, full.names = TRUE))
    if(!file.exists(sqfls[1])) {
      print("Extracting insert seqs..."); flush.console()
      foreach(I=fqfls, O=sqfls) %dopar% {
        system(sprintf("%s %s |egrep %s |cut -c10- > %s", "zcat", I, seqpat, O))
      }
    }
  }
  if(split) {
    print("Splitting inserts into chunks..."); flush.console()
    foreach(k = 1:2) %dopar% {
      prefix <- sprintf("%s/R%d_", chunkdir, k)
      system(sprintf("split -nl/%d %s %s", n.chunks, sqfls[k], prefix))
    }
  }
  if(!file.exists(bcfls[1])) {
    print("Extracting barcodes..."); flush.console()
    foreach(I=fqfls, O=bcfls) %dopar% {
      system(sprintf("%s %s |egrep %s |cut -c-8 > %s", "zcat", I, seqpat, O))
    }
  }
  if(!exists("good")) {
    print("Analyzing R1 barcodes..."); flush.console()
    bc1 <- read.delim(bcfls[1], header = FALSE, skipNul = TRUE)[,1]
    ta <- rev(sort(table(bc1)))[1:(top+1)]
    nn <- which.max(ta[1:top]/ta[2:(top+1)])
    BC1 <- names(ta)[1:nn]
    cat(nn, sprintf("%.2f%% good barcodes1\n", 100*sum(ta[1:nn])/sum(ta)))
    rm(ta)

    print("Analyzing R2 barcodes..."); flush.console()
    bc2 <- read.delim(bcfls[2], header = FALSE, skipNul = TRUE)[,1]
    ta <- rev(sort(table(bc2)))[1:(top+1)]
    nn <- which.max(ta[1:top]/ta[2:(top+1)])
    BC2 <- names(ta)[1:nn]
    cat(nn, sprintf("%.2f%% good barcodes2\n", 100*sum(ta[1:nn])/sum(ta)))
    rm(ta)

    print("Identifying good barcodes..."); flush.console()
    good <<- ((bc1 %fin% BC1) &
                (bc2 %fin% BC2) &
                (bc1 != bc2))
    Ng <- length(good)
    cat(sprintf("%.2f%% good barcode pairs in %d\n\n", 100*sum(good)/Ng, Ng))

    print("Constructing forward and reverse ditags..."); flush.console()
    FWD <<- paste0(bc1, bc2)
    REV <<- paste0(bc2, bc1)
    rm(bc1, bc2)
  }
  Isqs1 <- list.files(chunkdir, pattern = "R1_[a-z]{2}", full.names = TRUE)
  Isqs2 <- list.files(chunkdir, pattern = "R2_[a-z]{2}", full.names = TRUE)
  Isqs <- c(Isqs1, Isqs2)

  Osqs <- str_replace(Isqs, "[a-z]{2}$",
                      paste0(c(1:n.chunks, 1:n.chunks), ".fa"))

  print("Counting insert seqs in chunks..."); flush.console()
  cmd <- sprintf("wc -l %s", paste(Isqs, collapse = " "))
  Lens <- as.integer(str_extract(system(cmd, intern = TRUE), "[0-9]+"))
  Lens <- Lens[1:length(Isqs)]
  cumLens1 <- cumsum(Lens[1:n.chunks])
  cumLens2 <- cumsum(Lens[n.chunks + 1:n.chunks])
  cumLens <- c(cumLens1, cumLens2)

  print("Writing fasta files..."); flush.console()
  foreach(idx = isplitIndices(2*n.chunks, chunkSize = n.chunks/4)) %do% {
    print(idx); flush.console()
    lens <- Lens[idx]
    cumlens <- cumLens[idx]
    ofas <- Osqs[idx]
    isqs <- Isqs[idx]

    isqs <- foreach(isq = isqs) %dopar% {
      return(read.delim(isq, header = FALSE, skipNul = TRUE)[,1])
    }
    tmp <- foreach(inserts=isqs, ofa=ofas, len=lens, cumlen=cumlens) %dopar% {
      rng <- 1 + cumlen - (len:1)
      ok <- good[rng]
      ids <- as.integer(rng[ok])  # avoid float
      fanames <- paste0('sp', paste(ids, FWD[ids], REV[ids], sep = '.'))
      write.fasta(as.list(inserts[ok]), fanames, ofa, nbchar=300, as.string=TRUE)
      rm(fanames)
      return(sum(ok))
    }
    print(unlist(tmp))
  }
  file.remove(c(Isqs, bcfls, sqfls))
  rm(good, FWD, REV, inherits = TRUE)
}

mkfullpanel <- function(Sdir, clean, debug=FALSE, minsize=1e6) {
  if(clean) {
    file.remove(list.files(Sdir, pattern = ".+[.].+", full.names = TRUE))
  }
  if(!exists("psl.done")) psl.done <<- FALSE
  if(!exists("axt.done")) axt.done <<- FALSE

  for(idx in 1:length(probes)) {
    gr <- probes[idx]
    chrloc <- as.character(gr)

    xlsxfile <- sprintf("%s/%s.xlsx", Sdir, chrloc)
    if(file.exists(xlsxfile)) next               # already done here

    fafile <- sprintf("%s/%s_consensus.fa", Sdir, chrloc)
    if(file.exists(fafile)) {                    # resume from here
      mkcri(Sdir, gr)
      mkcri(Sdir, gr, mnv=FALSE)
      next
    }
    rdfile <- sprintf("%s/%s.RData", Sdir, chrloc)
    if(file.exists(rdfile)) {                    # resume from here
      if(file.size(rdfile) < minsize) next
      mkconsensus(Sdir, gr)
      mkcri(Sdir, gr)
      mkcri(Sdir, gr, mnv=FALSE)
      next
    }
    if(!psl.done) mkpsl(Sdir, gr)
    if(psl.done < 0) {
      psl.done <<- FALSE
      next
    }
    if(!axt.done) mkaxt(Sdir, gr)
    mkrdata(Sdir, gr, debug)
    if(file.size(rdfile) < minsize) next
    mkconsensus(Sdir, gr)
    mkcri(Sdir, gr)
    mkcri(Sdir, gr, mnv = FALSE)
  }
  mkallcri(Sdir)
  mkallcri(Sdir, mnv=FALSE)
  mkSNVcvg(PS, Sdir)
  mkSNVcvgsplits(Sdir)
}

mkpsl <- function(Sdir, gr) {

  axt.done <<- psl.done <<- FALSE

  S <- sub(".+/", "", Sdir)
  chrloc <- as.character(gr)
  print(paste("mkpsl", S, chrloc)); flush.console()
  chr <- as.character(seqnames(gr))
  start <- start(gr)
  end <- end(gr)

  chunkdir <- paste0(Sdir, "/chunks")

  fa1.chunks <- list.files(chunkdir, patt = "R1.+fa", full=TRUE)
  fa2.chunks <- list.files(chunkdir, patt = "R2.+fa", full=TRUE)
  fa12.chunks <- c(fa1.chunks, fa2.chunks)

  psl1.chunks <- sub(".fa$", paste0("_", chrloc, ".psl"), fa1.chunks)
  psl2.chunks <- sub(".fa$", paste0("_", chrloc, ".psl"), fa2.chunks)
  psl12.chunks <- c(psl1.chunks, psl2.chunks)

  z <- foreach(fa12=fa12.chunks, psl12=psl12.chunks) %dopar% {
    system(sprintf(blatcom,psdir,chr,start,end,fa12,blatopt,psl12))
    tmp <- vroom::vroom(psl12, delim = "\t", col_names = FALSE, col_types = ct)
    ok <- with(tmp, !(X11 + X12 - X13))
    tmp <- tmp[ok,]
    isdup <- duplicated(tmp$X10)
    if(any(isdup)) {
      ok <- with(tmp, !(X10 %in% unique(X10[isdup])))
      tmp <- tmp[ok,]
    }
    N <- nrow(tmp)
    if(!N) return(0)
    vroom::vroom_write(tmp, psl12, col_names = FALSE, quote = "none")
    return(N)
  }
  z <- unlist(z)
  if(any(!z)) {
    if(any(!!z)) file.remove(psl12.chunks[!!z])
    psl.done <<- -1
  } else psl.done <<- TRUE
}

mkaxt <- function(Sdir, gr) {
  axt.done <<- FALSE

  S <- sub(".+/", "", Sdir)
  chrloc <- as.character(gr)
  print(paste("mkaxt", S, chrloc)); flush.console()

  chr <- as.character(seqnames(gr))
  start <- start(gr)
  end <- end(gr)

  chunkdir <- paste0(Sdir, "/chunks")

  fa1.chunks <- list.files(chunkdir, patt = "R1.+fa", full=TRUE)
  fa2.chunks <- list.files(chunkdir, patt = "R2.+fa", full=TRUE)
  fa12.chunks <- c(fa1.chunks, fa2.chunks)

  psl1.chunks <- list.files(chunkdir,
                            patt = paste0("R1.+", chrloc, ".psl"), full=TRUE)
  psl2.chunks <- list.files(chunkdir,
                            patt = paste0("R2.+", chrloc, ".psl"), full=TRUE)
  psl12.chunks <- c(psl1.chunks, psl2.chunks)

  axt1.chunks <- sub("psl", "axt", psl1.chunks)
  axt2.chunks <- sub("psl", "axt", psl2.chunks)
  axt12.chunks <- c(axt1.chunks, axt2.chunks)

  foreach(psl12=psl12.chunks, fa12=fa12.chunks, axt12=axt12.chunks,
          .inorder=FALSE) %dopar% {
    system(sprintf(pslcom, psl12, psdir, chr, start, end, fa12, axt12))
  }
  file.remove(psl12.chunks)
  psl.done <<- FALSE
  axt.done <<- TRUE
}

mkoverlap <- function(df12) {
  attach(df12)

  s1 <- start
  e1 <- start + nchar(q1) -1
  e2 <- end
  s2 <- end - nchar(q2) +1

  a <- substr(q1, 1, s2-s1)
  c <- reverse(substr(reverse(q2), 1, e2-e1))
  b1 <- reverse(substr(reverse(q1), 1, 1+e1-s2))
  b2 <- substr(q2, 1, 1+e1-s2)

  detach(df12)

  pick <- which(b1 != b2)
  if(!!length(pick)) {
    bb1 <- strsplit(b1[pick], "", fixed=TRUE)
    bb2 <- strsplit(b2[pick], "", fixed=TRUE)

    b1[pick] <- foreach(x=bb1, y=bb2, .combine=c) %do% {
      x[x!=y] <- "N"
      return(paste(x, collapse = ""))
    }
    rm(b2,bb1,bb2)
  }
  return(paste0(a,b1,c))
}

mkcigar <- function(q) {
  pos <- str_locate_all(q, c("[NACGT]+","[-]+"))
  num <- lapply(pos, function(x) return(1 + x[,2] - x[,1]))
  n <- length(num[[2]])
  cig <- paste(sapply(1:n, function(k) {
    return(paste0(num[[1]][k], "M", num[[2]][k], "D"))
  }), collapse="")

  return(paste0(cig, num[[1]][n+1], "M"))
}

mkrdata <- function(Sdir, gr, debug=FALSE) {
  S <- sub(".+/", "", Sdir)
  chrloc <- as.character(gr)
  print(paste("mkrdata", S, chrloc)); flush.console()

  chunkdir <- paste0(Sdir, "/chunks")

  df12 <- sprintf("%s/%s.RData", Sdir, chrloc)

  chr <- as.character(seqnames(gr))
  start <- start(gr)
  end <- end(gr)

  axt1.chunks <- list.files(paste0(Sdir, "/chunks"),
                            patt = paste0("R1_.+_", chrloc, ".axt"), full=TRUE)
  axt2.chunks <- list.files(paste0(Sdir, "/chunks"),
                            patt = paste0("R2_.+_", chrloc, ".axt"), full=TRUE)
  axt12.chunks <- list(axt1.chunks, axt2.chunks)

  axt12 <- lapply(axt12.chunks, function(axts) {
    axt <- vroom::vroom(axts, delim = "\t", col_names = FALSE, col_types = "c")
    return(axt[, "X1", drop = TRUE])
  })
  print("Processing R1 and R2 axt files..."); flush.console()
  DT <- mclapply(axt12, function(a) {
    a_idx_info <- grep("^[1-9]", a)
    a_info <- a[a_idx_info]

    dt <- data.table::data.table(
      g = a[1 + a_idx_info],    # genomic reference sequence
      q = a[2 + a_idx_info],    # query sequence
      s = as.integer(sub("^([^ ]+ ){2}([0-9]+) .+", "\\2", a_info)),
      e = as.integer(sub("^([^ ]+ ){3}([0-9]+) .+", "\\2", a_info)),
      i = sub("^([^ ]+ ){4}([^.]+)[.].+", "\\2", a_info),
      f = sub("^([^ ]+ ){4}[^.]+[.]([^.]+)[.].+", "\\2", a_info),
      r = sub("^([^ ]+ ){4}[^.]+[.][^.]+[.]([^ ]+) .+", "\\2", a_info)
    )
    rm(a, a_idx_info, a_info)

    # remove all called insertions from q ("-" in g)
    idx <- grep("[-]", dt[,g])
    if(!!length(idx)) {
      G <- strsplit(dt[idx,g], split = "")
      Q <- strsplit(dt[idx,q], split = "")
      dt[idx, q := unlist(foreach(gg=G, qq=Q) %dopar% {
        paste(qq[gg != "-"], collapse = "")
      })]
      rm(G, Q)
    }
    return(dt)
  }, mc.cores = 2L)
  DT1 <- DT[[1]]
  DT2 <- DT[[2]]
  rm(DT, axt12)

  colnames(DT1) <- paste0(colnames(DT1), 1)
  colnames(DT2) <- paste0(colnames(DT2), 2)

  print("Identifying good read pairs..."); flush.console()

  if(debug) print("ids12...")
  ids12 <- intersect(DT1[,i1], DT2[,i2])

  if(debug) print("idx12...")
  idx12 <- data.frame(fmatch(ids12, DT1[,i1]), fmatch(ids12, DT2[,i2]))
  rm(ids12)

  # restrict analysis to paired reads
  if(debug) print("DT12...")
  DT12 <- cbind(DT1[idx12[,1]], DT2[idx12[,2]])

  rm(DT1, DT2, idx12)

  # assign mapped ranges and gaps to read pairs
  if(debug) print("start, end, gap...")
  DT12[,':='(start = pmin(s1, s2))]
  DT12[,':='(end = pmax(e1, e2))]
  DT12[,':='(gap = pmax(s1, s2) - pmin(e1, e2))]

  # determine unambiguous strand
  if(debug) print("strand...")
  DT12[,':='(strand = (sign(s2-s1) + sign(e2-e1))/2)]
  DT12 <- DT12[strand %in% c(-1, 1)]
  forward <- DT12[,strand] > 0

  # assign oriented duplex tags to read pairs
  if(debug) print("ditag...")
  DT12[,':='(mm = "X")]
  DT12[forward, mm := f1]
  DT12[!forward, mm := r1]

  # add cigar codes
  if(debug) print("cigar...")
  isdel <- grepl("[-]", DT12[,q1])
  DT12[,':='(c1 = "X")]
  DT12[!isdel, c1 := paste0(nchar(q1), "M")]
  DT12[isdel, c1 := unlist(foreach(q=q1) %dopar% mkcigar(q))]

  isdel <- grepl("[-]", DT12[,q2])
  DT12[,':='(c2 = "X")]
  DT12[!isdel, c2 := paste0(nchar(q2), "M")]
  DT12[isdel, c2:= unlist(foreach(q=q2) %dopar% mkcigar(q))]

  DF12 <- as.data.frame(DT12[, .(mm,q1,q2,c1,c2,start,end,gap,strand)])
  save(DF12, file = df12)

  rm(DT12, DF12)
  file.remove(unlist(axt12.chunks))
  psl.done <<- axt.done <<- FALSE
}

mkconsensus <- function(Sdir, gr, min.size=2, purity=3/4) {
  S <- sub(".+/", "", Sdir)
  chrloc <- as.character(gr)
  print(paste("mkconsensus", S, chrloc)); flush.console()

  outfile <- sprintf("%s/%s_consensus.fa", Sdir, chrloc)

  chr <- as.character(seqnames(gr))
  start <- start(gr)
  end <- end(gr)

  cl <- file(sprintf("%s/%s_consensus_countlog.txt",
                     Sdir, chrloc), open = "a")
  cat(date(), "\n\t", file = cl)

  infile <- sprintf("%s/%s.RData", Sdir, chrloc)
  load(infile)

  # construct dicigar codes
  f12 <- with(DF12, paste0(c1, c2))
  r12 <- with(DF12, paste0(c2, c1))
  cc <- ifelse(DF12$strand >0, f12, r12)

  # identify onTarget read families
  fam <- with(DF12, paste(mm, cc, start, end, gap, sep = ":"))
  DF.by.mlt <- split(DF12[c("q1","q2","strand")], fam)

  cat(nrow(DF12),
      "\tcorrectly structured read pairs uniquely mapped onto a CRI in",
      chr, "\n\t", file = cl)
  cat(length(DF.by.mlt),
      "\tonTarget ditag read families\n\t", file = cl)

  rm(DF12)

  # restrict families to preset size range
  DF.sizes <- sapply(DF.by.mlt, nrow)
  DF.pick <- which(DF.sizes > min.size)
  DF.by.mlt <- DF.by.mlt[DF.pick]

  cat(length(DF.pick),
      "\tonTarget ditag read families that have size >",
      min.size, "\n\t", file = cl)

  # count top and bottom strands
  DF.strand.freqs <- sapply(DF.by.mlt, function(x) table(x$strand))

  # library both strands
  DF.pick <- which(sapply(DF.strand.freqs, length) == 2)
  DF.strand.freqs <- DF.strand.freqs[DF.pick]
  DF.by.mlt <- DF.by.mlt[DF.pick]

  cat(length(DF.pick),
      "\tonTarget ditag read families that have size >",
      min.size, "and have both strands\n\t", file = cl)

  N <- 2*length(DF.by.mlt) - sum(grepl("[-]", names(DF.by.mlt)))
  cat(N, "\tonTarget ditag DCSs\n", file = cl)

  close(cl)

  # get duplex consensus sequences for Left and Right ends
  dcs <- foreach(x = DF.by.mlt) %dopar% {
    return(as.data.frame(sapply(c("L", "R"), function(lr) {
      sscs <- lapply(c(1, -1), function(s) {    # forward and reverse strands
        xs <- x[x$strand == s,]
        if(lr == "L") {
          q <- ifelse(s > 0, xs$q1, xs$q2)
          delta <- 0
        } else {
          q <- ifelse(s < 0, xs$q1, xs$q2)
          nc <- nchar(q)
          delta <- max(nc) - nc
        }
        rm(xs)
        cm <- consensusMatrix(q, shift = delta)
        cs <- colSums(cm)
        cm <- cm[, cs == length(q)]
        rm(q)
        return(consensusString(cm, ambig="N", thresh=purity))
      })
      q <- unlist(sscs)
      rm(sscs)
      delta <- 0
      if(lr == "R") {
        nc <- nchar(q)
        delta <- max(nc) - nc
      }
      cm <- consensusMatrix(q, shift = delta)
      cs <- colSums(cm)
      cm <- cm[, cs == 2]
      rm(q)
      return(consensusString(cm, ambig="N", thresh=purity))  # a single-end DCS
    }, simplify = FALSE)))
  }
  dcs <- as.data.frame(rbindlist(dcs))
  rownames(dcs) <- names(DF.by.mlt)

  # get unjoined consensus seqs
  dcs1 <- dcs[-grep("[-]", rownames(dcs)),]
  ujseq.L <- dcs1$L
  ujseq.R <- dcs1$R
  ujnames.L <- paste0(rownames(dcs1), ":L") # same root names for L & R seqs
  ujnames.R <- paste0(rownames(dcs1), ":R")

  # join overlapped consensus seqs
  dcs2 <- dcs[grep("[-]", rownames(dcs)),]
  colnames(dcs2) <- c("q1", "q2")
  info <- strsplit(rownames(dcs2), ":", fixed=TRUE)
  dcs2$start <- as.integer(sapply(info, "[", 3))
  dcs2$end <- as.integer(sapply(info, "[", 4))
  jseq.LR <- mkoverlap(dcs2)
  jnames.LR <- paste0(rownames(dcs2), ":LR")

  write.fasta(as.list(gsub("[-]", "", c(ujseq.L, ujseq.R, jseq.LR))),
              c(ujnames.L, ujnames.R, jnames.LR),
              outfile, nbchar=300, as.string=TRUE)

  rm(dcs, dcs1, dcs2, DF.by.mlt)
}

mkcri <- function(Sdir, gr, CRI=TRUE, maxlocs=5, mnv=TRUE,
                  targets=integer(0), chrloc=NA) {

  stopifnot(CRI | !!length(targets))

  S <- sub(".+/", "", Sdir)
  if(!is.na(chrloc)) {
    chr <- sub(":.+", "", chrloc)
    start <- as.integer(sub(".+:(.+)[-].+", "\\1", chrloc))
    end <- as.integer(sub(".+[-](.+)", "\\1", chrloc))
    GR <- GRanges(chr, IRanges(start, end))
  } else {
    chrloc <- as.character(gr)
    chr <- as.character(seqnames(gr))
    start <- start(gr)
    end <- end(gr)
  }
  print(paste("mkcri", S, chrloc)); flush.console()

  fa12 <- sprintf("%s/%s_consensus.fa", Sdir, chrloc)
  if(!file.exists(fa12)) {
    print("Unused probe")
    return()
  }
  axt12 <- sprintf("%s/%s.axt", Sdir, chrloc)
  psl12 <- sprintf("%s/%s.psl", Sdir, chrloc)

  system(sprintf(blatcom,psdir,chr,start,end,fa12,blatopt,psl12))

  tmp <- vroom::vroom(psl12, delim = "\t", col_names = FALSE, col_types = ct)
  ok <- with(tmp, !(X11 + X12 - X13))
  tmp <- tmp[ok,]
  isdup <- duplicated(tmp$X10)
  ok <- with(tmp, !(X10 %in% unique(X10[isdup])))
  tmp <- tmp[ok,]
  vroom::vroom_write(tmp, psl12, col_names = FALSE, quote = "none")
  rm(tmp)

  system(sprintf(pslcom, psl12, psdir, chr, start, end, fa12, axt12))
  axt <- vroom::vroom(axt12, delim = "\t", col_names = FALSE, col_types = "c")
  axt <- axt[, "X1", drop = TRUE]
  file.remove(psl12, axt12)

  axt_idx_info <- grep("^[1-9]", axt)
  axt_info <- strsplit(axt[axt_idx_info], split=" ")
  cri_refs <- strsplit(toupper(axt[1+axt_idx_info]), split="")
  cri_quer <- strsplit(toupper(axt[2+axt_idx_info]), split="")

  cri_starts <- as.integer(sapply(axt_info, "[", 3))  # relative to (start - flank)
  cri_ends <- as.integer(sapply(axt_info, "[", 4))    # relative to (start - flank)
  cri_mlts <- sapply(axt_info, "[", 5)
  cri_ir <- IRanges(cri_starts,cri_ends)

  CVG <- coverage(cri_ir)
  locN <- NULL

  cri <<- foreach(q=cri_quer, r=cri_refs, m=cri_mlts,
                  s=cri_starts, e=cri_ends) %do% {

                    # get rid of any insertions
                    ins <- r == "-"
                    if(any(ins)) {
                      q <- q[!ins]
                      r <- r[!ins]
                    }

                    if(!CRI) {
                      idx <- (1 +targets -s -start +flank)
                      if(any(idx < 1)) return(NULL)

                      rr <- paste(r[idx], collapse = "")
                      qq <- paste(q[idx], collapse = "")

                      mut <- sprintf("%s>%s", rr, qq)
                      ok <- !grepl('[-N]', mut)

                      return(mut[ok])
                    }

                    # locate Ns relative to probe range
                    if(any(q == "N")) {
                      locN <- c(locN, which(q == "N") + s-1)
                    }

                    # locate all variants
                    mm <- which(q != r)
                    if(!length(mm)) return(NULL)

                    # ignore variants within mindist of an end
                    ok <- (mm > mindist) & (mm <= length(q) - mindist)
                    if(!any(ok)) return(NULL)
                    mm <- mm[ok]

                    # ignore Ns in reference
                    ok <- r[mm] != "N"
                    if(!any(ok)) return(NULL)
                    mm <- mm[ok]

                    # ignore Ns in query
                    ok <- q[mm] != "N"
                    if(!any(ok)) return(NULL)
                    mm <- mm[ok]

                    # exclude apparently mismapped contaminants
                    ndel <- sum(q[mm] == "-")
                    if(length(mm) > maxlocs + ndel) return(NULL)

                    locs <- mm + s + start - flank - 1   # chromosomal coordinates

                    if(length(mm) > 1) {
                      if(mnv) {
                        Locs <- paste(locs, collapse = "#")
                        Muts <- paste(paste(r[mm], collapse = ""),
                                      paste(q[mm], collapse = ""),
                                      sep = ">")
                        res <- paste(Locs, ".", Muts, ".", sep = ":")
                      } else {
                        tns <- sapply(mm, function(m) paste(r[(m-1):(m+1)],
                                                            collapse = ""))
                        muts <- paste(r[mm], q[mm], sep = ">")
                        cvgs <- CVG[mm+s-1]
                        res <- paste(locs, tns, muts, cvgs, sep = ":")
                      }
                    } else {
                      tns <- paste(r[(mm-1):(mm+1)], collapse = "")
                      mut <- paste(r[mm], q[mm], sep = ">")
                      cvg <- CVG[mm+s-1]
                      res <- paste(locs, tns, mut, cvg, sep = ":" )
                    }
                    return(res)
                  }
  if(all(sapply(cri, is.null))) {
    if(!CRI) {
      print(sprintf("%s %s %s", S, chr, paste(targets, collapse=", ")))
      print("No mutations.")
      return(0)
    }
    fname <- sprintf("%s/%s%s.xlsx", Sdir, ifelse(mnv,"","snv_"), chrloc)
    DF <- data.frame(chr=NA, location=NA, triN=NA, mutation=NA,
                     count=NA, coverage=NA, pctN=NA, gene=NA, ROI=NA)
    write.xlsx(DF, fname, rowNames=F, overwrite = TRUE)
    return(0)
  }
  names(cri) <<- cri_mlts
  present <- function(x) !is.na(x)
  cri_locmut <- Filter(present, unlist(cri))
  cri_tab <<- table(cri_locmut)

  if(!CRI) {
    DF <- as.data.frame(cri_tab)
    DF[,1] <- as.character(DF[,1])
    print(sprintf("%s %s %s", S, chr, paste(targets, collapse=" ")))
    print(DF, row.names = FALSE)
  } else {
    cnts <- as.integer(cri_tab)
    info <- strsplit(names(cri_tab), split = ":")
    locs <- sapply(info, "[", 1)
    tris <- sapply(info, "[", 2)
    muts <- sapply(info, "[", 3)
    cvgs <- sapply(info, "[", 4)

    ok <- !grepl("#",names(cri_tab))
    Locs <- as.integer(locs[ok])
    Cvgs <- as.integer(cvgs[ok])

    cvgs[ok] <- Cvgs
    is.na(cvgs[!ok]) <- TRUE

    full_width <- end - start + 2*flank + 4  # extra width for safety

    if(is.null(locN)) {
      freqN <- Rle(0, full_width)
    } else freqN <- Rle(sort(locN))

    countN <- Rle(0, full_width)
    countN[runValue(freqN)] <- runLength(freqN)

    pctN <- rep(NA, length(ok))
    pctN[ok] <- round(100*as.integer(countN[Locs - start + flank])/Cvgs, 3)

    DF <- data.frame(chr=chr, location=locs, triN=tris, mutation=muts,
                     count=cnts, coverage=cvgs, pctN=pctN, gene=NA, ROI=NA)

    if(mnv) {
      LOCS <- lapply(strsplit(locs, split="#"), as.integer)
    } else {
      LOCS <- lapply(locs, as.integer)
    }
    LOCS.gr <- lapply(LOCS, function(x) GRanges(chr, IRanges(x, width=1)))
    if(!!length(gr$info)) {
      idx <- which(!!sapply(LOCS.gr, function(x) sum(countOverlaps(x, GR))))
    } else {
      idx <- which(!!sapply(LOCS.gr, function(x) sum(countOverlaps(x, gr))))
    }
    if(!!length(idx)) {
      if(!!length(gr$info)) {
        DF$gene[idx] <- paste(unique(sub("_ROI.+", "", gr$info)), collapse="//")
        ovl <- lapply(LOCS.gr[idx], function(x) findOverlaps(x, gr))
        if(any(!!lengths(ovl))) {
          id1 <- which(!!lengths(ovl))
          id2 <- lapply(ovl[id1], function(x) unique(to(x)))
          DF$ROI[idx[id1]] <- sapply(id2, function(ids) {
            paste(gr$info[ids], collapse = ", ")
          })
        }
      } else {
        DF$gene[idx] <- gr$gene
        ovl <- lapply(LOCS.gr[idx], function(x) findOverlaps(x, probe.info.gr,
															 maxgap = -1))
        if(any(!!lengths(ovl))) {
          id1 <- which(!!lengths(ovl))
          id2 <- lapply(ovl[id1], function(x) unique(to(x)))
          DF$ROI[idx[id1]] <- sapply(id2, function(ids) {
            paste(probe.info.gr$ROI[ids], collapse = ", ")
          })
        }
      }
    }
    fname <- sprintf("%s/%s%s.xlsx", Sdir, ifelse(mnv,"","snv_"), chrloc)
    write.xlsx(DF, fname, rowNames=F)
  }
}

mkallcri <- function(Sdir, mnv=TRUE, group="demo", maxlocs=4, WB=TRUE,
                     BG=FALSE, normalize=FALSE, redo=FALSE) {
  if(normalize) WB <- FALSE

  S <- sub(".+/", "", Sdir)

  if(WB) {
    wbfile <- paste0(fastq, group, ifelse(mnv,"_m","_s"), "nv_panels.xlsx")
    if(file.exists(wbfile)) {
      wb <- loadWorkbook(wbfile)
      if(S %in% names(wb)) return()
    } else {
      wb <- createWorkbook()
    }
  }
  if(redo) {
    xlsxfiles <- list.files(Sdir, patt="^chr.+xlsx", full.names = TRUE)
    chrlocs <- sub(".+(chr.+).xlsx", "\\1", xlsxfiles)
    chr <- sub(":.+","",chrlocs)
    start <- as.integer(sub(".+:(.+)[-].+", "\\1", chrlocs))
    end <- as.integer(sub(".+[-](.+)", "\\1", chrlocs))
    GR <- GRanges(chr, IRanges(start, end))

    ps.gr <- PS[[group]]
    invisible(foreach(k=1:length(GR), .inorder=FALSE) %dopar% {
      ovl <- as.matrix(findOverlaps(ps.gr, GR[k]))
      mkcri(Sdir=Sdir, gr=ps.gr[ovl[,1]], mnv=mnv, chrloc=chrlocs[k])
    })
  }
  pat <- "chr.+xlsx"
  xlsxfiles <- list.files(Sdir, patt=pat, full.names=TRUE)
  pat <- paste0(ifelse(mnv, "/", "/snv_"), pat)
  xlsxfiles <- xlsxfiles[grep(pat, xlsxfiles)]

  chrlocs <- sub(".+(chr.+).xlsx", "\\1", xlsxfiles)
  chr <- sub(":.+","",chrlocs)
  start <- as.integer(sub(".+:(.+)[-].+", "\\1", chrlocs))
  end <- as.integer(sub(".+[-](.+)", "\\1", chrlocs))

  ord <- naturalorder(chrlocs)

  DF <- data.frame()
  for(xlsxfile in xlsxfiles[ord]) {
    df <- read.xlsx(xlsxfile)
    if(!nrow(df)) next
    idx <- which(!is.na(df$gene)) # no flanks
    if(!length(idx)) next
    if(!mnv) {
      cvg <- as.numeric(df$coverage[idx])
      cnt <- as.numeric(df$count[idx])
      idx <- idx[3*cnt < cvg]     # no SNPs
      if(!length(idx)) next
    }
    DF <- rbind(DF, df[idx,])
  }
  if(!mnv) {
    DF$location <- as.integer(DF$location)
    DF$coverage <- as.integer(DF$coverage)
    if(normalize) DF$count <- round(DF$count * 50000/max(DF$coverage), 2)
  }
  if(WB) {
    addWorksheet(wb, S)
    writeData(wb, S, DF)
    setHeaderFooter(wb, S)
    setColWidths(wb, S, 1:ncol(DF), widths = "auto")
    freezePane(wb, S, firstActiveRow = 2, firstActiveCol = 1)
    saveWorkbook(wb, wbfile, TRUE)
  }
  if(!BG) return()

  if(mnv)  DF <- DF[str_count(DF$location, fixed("#")) < maxlocs,]

  if(mnv) {
    bg <- tidyr::separate_rows(DF[c("chr", "location", "count")],
                               location, sep = "#", convert = TRUE)
    bg$end <- as.integer(bg$location)
    bg$start <- bg$end -1L
    bg <- bg[c("chr", "start", "end", "count")]
    bg <- bg %>% group_by(chr, start, end) %>% summarize(count = sum(count))
    if(normalize) {
      bgfile <- sprintf("%s/normalized_%s%s", Sdir, S, ".bedgraph")
    } else bgfile <- sprintf("%s/%s%s", Sdir, S, ".bedgraph")

    if(file.exists(bgfile)) file.remove(bgfile)
    bgf <- file(bgfile, open = "a")
    cat("track type=bedGraph\n", file = bgf)
    write.table(bg, bgf, TRUE, col.names=F, row.names=F, quote=F, sep="\t")
    close(bgf)
  } else {
    bg <- DF[c("chr", "location", "count")]
    bg$end <- as.integer(bg$location)
    bg$start <- bg$end -1L
    bg <- bg[c("chr", "start", "end", "count")]

    is.CT <- grepl("C>T|G>A", DF$mutation)
    is.diPy <- !(grepl("[GA]C[GA]",DF$triN) | grepl("[CT]G[CT]",DF$triN))

    pick <- is.CT & is.diPy
    if(any(pick)) {
      bg1 <- bg[pick,]
      if(normalize) {
        bgfile1 <- sprintf("%s/normalized_snv1_%s%s", Sdir, S, ".bedgraph")
      } else bgfile1 <- sprintf("%s/snv1_%s%s", Sdir, S, ".bedgraph")

      if(file.exists(bgfile1)) file.remove(bgfile1)
      bgf1 <- file(bgfile1, open = "a")
      cat("track type=bedGraph\n", file = bgf1)
      write.table(bg1, bgf1, TRUE, col.names=F, row.names=F, quote=F, sep="\t")
      close(bgf1)
    }
    pick <- is.CT & !is.diPy
    if(any(pick)) {
      bg2 <- bg[pick,]
      if(normalize) {
        bgfile2 <- sprintf("%s/normalized_snv2_%s%s", Sdir, S, ".bedgraph")
      } else bgfile2 <- sprintf("%s/snv2_%s%s", Sdir, S, ".bedgraph")

      if(file.exists(bgfile2)) file.remove(bgfile2)
      bgf2 <- file(bgfile2, open = "a")
      cat("track type=bedGraph\n", file = bgf2)
      write.table(bg2, bgf2, TRUE, col.names=F, row.names=F, quote=F, sep="\t")
      close(bgf2)
    }
    pick <- !is.CT
    if(any(pick)) {
      bg3 <- bg[pick,]
      if(normalize) {
        bgfile3 <- sprintf("%s/normalized_snv3_%s%s", Sdir, S, ".bedgraph")
      } else bgfile3 <- sprintf("%s/snv3_%s%s", Sdir, S, ".bedgraph")

      if(file.exists(bgfile3)) file.remove(bgfile3)
      bgf3 <- file(bgfile3, open = "a")
      cat("track type=bedGraph\n", file = bgf3)
      write.table(bg3, bgf3, TRUE, col.names=F, row.names=F, quote=F, sep="\t")
      close(bgf3)
    }
  }
}
