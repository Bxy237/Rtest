#make a package: vizAPA
#

# Whether A in B
# AinB(c('start','end'),c('start','x','y'),all=T)
# @param A a vector
# @param B a vector
# @param all whether need all A in B
# @return TRUE/FALSE
# @examples
# AinB(c('start','end'),c('start','x','y'),all=T)
AinB<-function(A, B, all=T){
  if (is.factor(A))   A=as.character(A)
  if (is.factor(B))   B=as.character(B)
  if (!(is.vector(A) & is.vector(B))) {
    return(F)
  }
  x=sum(A %in% B)
  if ((all & x==length(A)) | (!all & x>0)) {
    return(T)
  } else {
    return(F)
  }
}




# -------- class PACdataset -------------


#' PACdataset object
#'
#' PACdataset is an S4 class to represent a poly(A) site list.
#' PACs are stored as PACdataset for filtering, statistics, and other operations.
#'
#' @slot counts a data frame denoting the counts of PAs. Each row is one PA, each column is one experiment.
#' @slot colData sample defination. Each row is one experiment, and the column is the sample group.
#' The number of rows is the same as the column number of the counts slot.
#' @slot anno detailed annotation information for PAs, with the same row number as the counts slot.
#' @slot supp a list storing additional data, like stopCodon.
#' @family PACdataset functions
PACdataset <- setClass("PACdataset", slots=list(counts="data.frame", colData="data.frame", anno="data.frame", supp="list"))

setMethod("show",
          "PACdataset",
          function(object) {
            cat('PAC#',nrow(object@counts),'\n')

            if ('ftr' %in% colnames(object@anno)) {
              cat('gene#',length(unique(object@anno$gene[getNonItgFtrId(object@anno$ftr)])),'\n')
              pn=table(object@anno$ftr)
              pn=cbind(nPAC=pn)
              print(pn)

              if ('3UTR' %in% unique(object@anno$ftr)) {
                avgutr=NA
                if ('toStop' %in% colnames(object@anno)) {
                  avgutr=mean(object@anno$toStop[object@anno$ftr=='3UTR'])
                } else if ('three_UTR_length' %in% colnames(object@anno)) {
                  avgutr=mean(object@anno$three_UTR_length[object@anno$ftr=='3UTR'])
                }
                if (!is.na(avgutr)) cat('Mean 3UTR length of PACs (bp):',floor(avgutr),'\n')
              }
            }

            cat('sample#',ncol(object@counts),'\n')
            cat(colnames(object@counts)[1:(min(ncol(object@counts),5))],'...\n')
            cat('groups:\n')
            cat(sprintf('%s...[%d x %d]%s','@colData',nrow(object@colData), ncol(object@colData),'\n'))
            if (nrow(object@colData)>0)
              print(head(object@colData[, 1:min(10, ncol(object@colData)), drop=F], 2))
            cat(sprintf('%s...[%d x %d]%s','@counts',nrow(object@counts), ncol(object@counts),'\n'))
            if (nrow(object@counts)>0)
              print(head(object@counts[, 1:min(10, ncol(object@counts)), drop=F],2))
            cat(sprintf('%s...[%d x %d]%s','@colData',nrow(object@colData), ncol(object@colData),'\n'))
            print(head(object@colData,2))
            cat(sprintf('%s...[%d x %d]%s','@anno',nrow(object@anno), ncol(object@anno),'\n'))
            print(head(object@anno,2))

            if (length(object@supp)>0) {
              cat(sprintf('%s...[%d]%s','@supp',length(object@supp),'\n'))
              cat(names(object@supp),'\n')
            }

          }
)

#nrows of PACdataset
setMethod("length",
          "PACdataset",
          function(x) {
            return(nrow(x@counts))
          }
)

#overloading of subscript operator
#pacds[1:100, 2:5]
#pacds[1:100, ]
#pacds[, 1:5]
#pacds[, c(T,T,F)]
#pacds[, c('I_18','I_19')]
setMethod("[", signature(x = "PACdataset"), def=function(x, i, j, ..., drop=FALSE) {
  if (!missing(i)) {
    x@counts=x@counts[i, , drop=drop]
    x@anno=x@anno[i,, drop=drop ]
  }
  if (!missing(j)) {
    x@counts=x@counts[, j, drop=drop ]
    cnames=colnames(x@counts)
    x@colData=x@colData[cnames, , drop=F]
  }
  for (i in 1:ncol(x@colData)) {
    x@colData[,i]=factor(x@colData[,i])
  }
  return(x)
}
)

## rbind multiple PACds
rbind.PACds <- function(...) {
  pars=list(...)
  if (length(pars)==0) return(pars)
  if (length(pars)==1) return(pars[[1]])
  if (sum(lapply(pars, class)=='PACdataset')!=length(pars)) {
    stop("rbind: not all vars are PACdataset!\n")
  }
  p=pars[[1]]
  for (i in 2:length(pars)) {
    p2=pars[[i]]
    if (ncol(p@counts)!=ncol(p2@counts)) stop("rbind: not same column number in @counts!")
    if (ncol(p@anno)!=ncol(p2@anno)) stop("rbind: not same column number in @anno!")
    if (!identical(p@colData,p2@colData)) stop("rbind: not the same @colData!")
    p@counts=rbind(p@counts, p2@counts)
    p@anno=rbind(p@anno, p2@anno)
    if (!identical(p@supp, p2@supp)) p@supp=list(p@supp, p2@supp)
  }
  return(p)
}


#' Combine multiple PACdatasets by rows
#'
#' @usage rbind(...)
#' @param ... one or more PACdataset object.
#' @return If params are PACdatasets with the same column definations, then will combine rows of PACds. See ?base::cbind for the value returned by the default methods.
#' @name rbind
#' @family PACdataset functions
#' @export
rbind <- function (...) {
  if (length(attr(list(...)[[1]], "class"))>0) {
    if (attr(list(...)[[1]], "class") == "PACdataset") return(rbind.PACds(...))
  }
  return(base::rbind(...))
}


# -------------- *** readPACdataset() **** ---------------


#' Read a PACdataset
#'
#' readPACdataset reads PAC counts and sample annotation into a PACdataset.
#'
#' @usage readPACdataset(pacFile, colDataFile, noIntergenic=TRUE, PAname='PA')
#' @param pacFile a file name or a data frame. If it is a file, it should have header, with at least (chr, strand, coord) columns.
#' This file could have other columns, including gff cols (gene/gene_type/ftr/ftr_start/ftr_end/UPA_start/UPA_end) and user-defined sample columns.
#' If there are at least one non-numeric columns other than above gff cols, then all remaining columns are considered as annotation columns.
#' If all remaining columns are numeric, then they are all treated as sample columns.
#' Use annotatePAC() first if need genome annotation of coordinates.
#' @param colDataFile a file name or a data frame. If it is a file, then it is an annotation file of samples with header,
#' rownames are samples (must be all in pacFile), columns names are sample groups.
#' There could be single or multiple columns to define the groups of samples.
#' When colDataFile=NULL, then readPACds will automately retreive sample columns and gff columns (if any) from pacFile.
#' If there is no sample columns, then will set colData as a data frame with 1 column (=group) and 1 row (=tag), and element=group1.
#' If pacfile or colDataFile is a character, then it is a file name, so readPACds will read data from file.
#' @param noIntergenic TRUE/FALSE. If TRUE, then will remove PACs in intergenic (ftr='^inter')
#' @param PAname specify how to set the name (rowname) of PACs.
#' PAname=PA, the PA name is set as 'gene:PAN'; PAname=coord, then 'gene:coord'.
#' @return A PACdataset object, with @anno being a data frame with at least three columns chr/strand/coord. If there is no sample column, then will add one sample named tag in group1.
#' @examples
#' data(PACds)
#' ## read simple PACfile that only has columns chr/strand/coord
#' pacFile=PACds@anno[,c('chr','strand','coord')]
#' colDataFile=NULL
#' p=readPACdataset(pacFile, colDataFile)

#' ## read PACfile that has columns chr/strand/coord and sample columns
#' pacFile=PACds@anno[,c('chr','strand','coord')]
#' pacFile=cbind(pacFile, PACds@counts[,c('anther1','embryo1','anther2')])
#' colDataFile=NULL
#' p=readPACdataset(pacFile, colDataFile)
#' p@colData; head(p@counts)

#' ## read PACfile that has columns chr/strand/coord, sample columns, and gff cols like gene/gene_type/ftr/ftr_start/ftr_end/UPA_start/UPA_end
#' pacFile=PACds@anno
#' pacFile=cbind(pacFile, PACds@counts[,c('anther1','embryo1','anther2')])
#' colDataFile=NULL
#' p=readPACdataset(pacFile, colDataFile)
#' p@colData; head(p@counts); head(p@anno)

## read from data frame of PACfile and colDataFile
#' pacFile=PACds@anno
#' smps=c('anther1','embryo1','anther2')
#' pacFile=cbind(pacFile, PACds@counts[,smps])
#' colDataFile=as.data.frame(matrix(c('group1','group2','group1'), ncol=1, dimnames=list(smps, 'group')))
#' p=readPACdataset(pacFile, colDataFile)
#' p@colData; head(p@counts); head(p@anno)

## read from file names of PACfile and colDataFile
#' write.table(pacFile, file='pacFile', row.names=FALSE)
#' write.table(colDataFile, file='colDataFile', row.names=TRUE)
#' p=readPACdataset(pacFile='pacFile',
#'             colDataFile='colDataFile', noIntergenic=TRUE, PAname='PA')
#'
#' @name readPACdataset
#' @export
readPACdataset<-function(pacFile, colDataFile=NULL, noIntergenic=TRUE, PAname='PA') {

  if (!(PAname %in% c('PA','coord'))) {
    stop("PAname must be PA or coord")
  }

  if (is.character(pacFile)) {
    d=read.table(pacFile, header=T)
  } else {
    d=pacFile
  }

  gffcols=c('gene','gene_type','ftr','ftr_start','ftr_end','UPA_start','UPA_end')

  if (sum(!(c('chr','strand','coord') %in% colnames(d)))!=0) {
    stop("chr,strand,coord not all in header of pacfile")
  }

  cat(nrow(d),'PACs\n')

  if ('ftr' %in% colnames(d)) {
    if (noIntergenic) {
      d=d[getNonItgFtrId(d$ftr),]
      cat(nrow(d),'No intergenic PACs\n')
    }
    #WB's data 20190620
    d$ftr[d$ftr=='three_prime_UTR']='3UTR'
    d$ftr[d$ftr=='five_prime_UTR']='5UTR'
  }

  #d[d=='unkown']=NA

  if ('ftr_start' %in% colnames(d)) {
    if (!is.numeric(d$ftr_start)) {d$ftr_start=as.numeric(d$ftr_start)}
    if (!is.numeric(d$ftr_end)) {d$ftr_end=as.numeric(d$ftr_end)}
  }

  if ('gene' %in% colnames(d)) {
    idx=which(d$gene=='NULL')
    if (length(idx)>0) {
      cat(length(idx),'gene name is NULL, change to chrStrand')
      d$gene[idx]=paste0(d$chr[idx],d$strand[idx])
    }

    #order 5' to 3'
    d1=d[d$strand=='+',]
    d1=d1[order(d1$gene,d1$coord,decreasing = FALSE),]
    d2=d[d$strand=='-',]
    d2=d2[order(d2$gene,d2$coord,decreasing = TRUE),]
    d=rbind(d1,d2)

    if (PAname=='coord') {
      paid=paste0(d$gene,':',d$coord)
    } else if (PAname=='PA') {
      rg=rle(d$gene)
      paid=paste0(d$gene,':PA',unlist(lapply(rg$lengths,seq)))
    }
  } else {
    paid=paste0('PA',1:nrow(d))
  }

  rownames(d)=paid

  allcols=colnames(d)


  # group defination of sample columns
  if (!is.null(colDataFile)) {
    if (is.vector(colDataFile)) {
      colData=read.table(colDataFile, colClasses="character")
    } else {
      colData=colDataFile
    }

    if (sum(rownames(colData) %in% colnames(d))!=nrow(colData)) {
      stop("rownames of annofile not all in columns of pacfile")
    }

  } else {
    #remove gffcolid and chr/strand/coord, other columns are sample columns, and the sample group is 'group', groupname is 'group1'
    smpcols=which(allcols %in% c(gffcols,'chr','strand','coord'))
    smpcols=allcols[-smpcols]
    #no tag columns, then add tag=1 columns
    if (length(smpcols)==0) {
      smpcols='tag'
      d$tag=1
    } else { #one columns is chr, then they are all annotations
      for (i in smpcols) {
        if (!(is.numeric(d[,i]))) {
          smpcols='tag'
          d$tag=1
          break
        }
      }
    }
    colData=as.data.frame(matrix( rep('group1',length(smpcols)), ncol=1, dimnames =list(smpcols,'group') ))
  }

  for (i in 1:ncol(colData)) {
    colData[,i]=factor(colData[,i])
  }
  colData=colData[order(rownames(colData)), , drop=F]
  anno=d[,-which(colnames(d) %in% rownames(colData))]
  anno[anno=='unkown']=NA

  counts=d[, rownames(colData), drop=F]
  PACds=new("PACdataset",counts=counts, colData=colData, anno=anno)
  return(PACds)
}





# -------------- *** useGff() **** ---------------

#' Parse genome annotation file
#'
#'  parse genome annotation file of  gtf/gff3/gff format or genome annotation object of TxDb.
#'  If specify species, then specify the  source of gff file from BioMart.
#' @usage useGff(gff, txdb,species)
#' @param gff annotation gtf/gff3/gff fotmat  file
#' @param txdb a TxDb object
#' @param species  specify species, then specify the  source of gff file from BioMart.
#' @return If input gtf.file,then return parsed genome annotation object, which is a data.frame
#' If input txdb ,then return parsed genome annotation object, which is a list of three elements (anno.rna, anno.need, anno.frame).
#' If input a species, then return gff file from BioMart.
#' @examples
#' ## Way1: Based on an annotation file in gtf format
#' gff.path <- "Arabidopsis_thaliana.TAIR10.49.gtf"
#' gtf.data <- useGff(gff =gff.path)
#'
#' ##Way2: Based on a TxDb object generated from BioMart.
#' # Parse Arabidopsis Txdb
#' library(TxDb.Athaliana.BioMart.plantsmart28)
#' txdb <- useGff(TxDb.Athaliana.BioMart.plantsmart28)
#' # Parse mm10 Txdb
#' BiocManager::install("TxDb.Mmusculus.UCSC.mm10.ensGene")
#' library(TxDb.Mmusculus.UCSC.mm10.ensGene)
#' txdb <- useGff(TxDb.Mmusculus.UCSC.mm10.ensGene)
#' @name useGff
#' @import biomaRt
#' @export
useGff<- function (gff=NULL,txdb=NULL){
  if (is.null(gff) & is.null(txdb)) {
    return("must provide gff file or Txdb objict")
  }
  if (!is.null(gff) & !is.null(txdb)){
    return("can't provide gff file and Txdb object at the same time, choose one to provide")
  }
  if (!is.null(gff)) {
    if (grepl('\\.gff3|\\.gtf',tolower(gff))) {
      rt=parseGff(gff)   #当输入gff3/gff/gtf格式的文件时，使用底层函数parseGff加载文件

    }
    invisible(gc())
    return(rt)}
  #~chr
  if(!is.null(txdb)){
    if (class(txdb)=='TxDb'){rt=parseTxdb(txdb)}#当输入"Txdb"对象时，使用底层函数parseTxdb
    invisible(gc())
    return(rt)
  }
}

##function:parseTxdb
##parses genome annotation object of TxDb
parseTxdb <- function (an.txdb) {

  if(class(an.txdb)!='TxDb') stop("an.txdb not of class TxDb!")

  genes <- genes(an.txdb,columns=c("tx_type","gene_id"))
  genes <- as.data.frame(genes)

  genes <- data.frame(seqnames=as.character(genes$seqnames) ,start=as.integer(genes$start),
                      end=as.integer(genes$end),width=as.integer(genes$width),
                      strand=as.character(genes$strand),type="gene",
                      ID =as.character(genes$gene_id),biotype=as.character(genes$tx_type),
                      gene_id =as.character(genes$gene_id),Parent=NA,transcript_id=NA)
  #setdiff(colnames(genes),colnames(tari))
  rnas <- transcripts(an.txdb,columns=c("tx_name","tx_type","gene_id"))
  rnas<- as.data.frame(rnas)

  #test <- strsplit(as.character(rnas$gene_id) ,"\\s+")
  #temp3 <- paste("",lapply(test,"[[",1),sep="");
  #head(temp3)
  rnas <- data.frame(seqnames=as.character(rnas$seqnames) ,start=as.integer(rnas$start),
                     end=as.integer(rnas$end),width=as.integer(rnas$width),
                     strand=as.character(rnas$strand),type="RNA",
                     ID =as.character(rnas$tx_name),biotype=as.character(rnas$tx_type),
                     gene_id =as.character(rnas$gene_id),Parent=as.character(rnas$gene_id),
                     transcript_id=as.character(rnas$tx_name))
  # exons <- exons(an.txdb,columns=c("exon_name","tx_name","tx_type","gene_id"))
  # exons <- as.data.frame(exons)
  # head(exons)

  exons <- exonsBy(an.txdb,by=c("tx"),use.names=TRUE)
  exons <- as.data.frame(exons)
  exons <- data.frame(seqnames=as.character(exons$seqnames) ,start=as.integer(exons$start),
                      end=as.integer(exons$end),width=as.integer(exons$width),
                      strand=as.character(exons$strand),type="exon",
                      ID =as.character(exons$exon_name),biotype=NA,
                      gene_id =NA,Parent=as.character(exons$group_name),
                      transcript_id=as.character(exons$group_name))
  index <- match(exons$Parent,rnas$transcript_id)
  #which(is.na(index))
  exons$gene_id <- rnas$Parent[index]
  exons$biotype <- rnas$biotype[index]

  #==================================
  #CDS
  cdss <- cdsBy(an.txdb,by=c("tx"),use.names=TRUE)
  cdss <- as.data.frame(cdss)
  cdss <- data.frame(seqnames=as.character(cdss$seqnames) ,start=as.integer(cdss$start),
                     end=as.integer(cdss$end),width=as.integer(cdss$width),
                     strand=as.character(cdss$strand),type="CDS",
                     ID =as.character(cdss$cds_name),biotype=NA,
                     gene_id =NA,Parent=as.character(cdss$group_name),
                     transcript_id=as.character(cdss$group_name))
  index <- match(cdss$Parent,rnas$transcript_id)
  #which(is.na(index))
  cdss$gene_id <- rnas$Parent[index]
  cdss$biotype <- rnas$biotype[index]
  #head(cdss)
  #cdss <- cds(an.txdb,columns=c("cds_name","tx_name","tx_type","gene_id"))

  #==================================
  #introns
  introns <- intronsByTranscript(an.txdb,use.names=TRUE)
  introns <- as.data.frame(introns)
  introns <- data.frame(seqnames=as.character(introns$seqnames) ,start=as.integer(introns$start),
                        end=as.integer(introns$end),width=as.integer(introns$width),
                        strand=as.character(introns$strand),type="intron",
                        ID =NA,biotype=NA,
                        gene_id =NA,Parent=as.character(introns$group_name),
                        transcript_id=as.character(introns$group_name))
  index <- match(introns$Parent,rnas$transcript_id)
  #which(is.na(index))
  introns$gene_id <- rnas$Parent[index]
  introns$biotype <- rnas$biotype[index]
  #head(introns)

  #===================================================
  #five UTR
  fiveUTRs <- fiveUTRsByTranscript(an.txdb,use.names=TRUE)
  fiveUTRs <- as.data.frame(fiveUTRs)
  fiveUTRs <- data.frame(seqnames=as.character(fiveUTRs$seqnames) ,start=as.integer(fiveUTRs$start),
                         end=as.integer(fiveUTRs$end),width=as.integer(fiveUTRs$width),
                         strand=as.character(fiveUTRs$strand),type="five_prime_UTR",
                         ID =NA,biotype=NA,
                         gene_id =NA,Parent=as.character(fiveUTRs$group_name),
                         transcript_id=as.character(fiveUTRs$group_name))
  index <- match(fiveUTRs$Parent,rnas$transcript_id)
  #which(is.na(index))
  fiveUTRs$gene_id <- rnas$Parent[index]
  fiveUTRs$biotype <- rnas$biotype[index]
  #head(fiveUTRs)


  #===========================================
  #three UTR
  threeUTRs <- threeUTRsByTranscript(an.txdb,use.names=TRUE)
  threeUTRs <- as.data.frame(threeUTRs)
  threeUTRs <- data.frame(seqnames=as.character(threeUTRs$seqnames) ,start=as.integer(threeUTRs$start),
                          end=as.integer(threeUTRs$end),width=as.integer(threeUTRs$width),
                          strand=as.character(threeUTRs$strand),type="three_prime_UTR",
                          ID =NA,biotype=NA,
                          gene_id =NA,Parent=as.character(threeUTRs$group_name),
                          transcript_id=as.character(threeUTRs$group_name))
  index <- match(threeUTRs$Parent,rnas$transcript_id)
  #which(is.na(index))
  threeUTRs$gene_id <- rnas$Parent[index]
  threeUTRs$biotype <- rnas$biotype[index]

  anno.frame <- rbind(genes,rnas,exons,cdss,introns,fiveUTRs,threeUTRs)
  anno.frame$type <- factor(anno.frame$type,levels=c("gene","RNA","five_prime_UTR","exon","CDS","intron",
                                                     "three_prime_UTR"))
  #anno.frame <- anno.frame[order(anno.frame$transcript_id,anno.frame$gene_id,
  #                               anno.frame$start,anno.frame$strand,anno.frame$type),]
  anno.need <- rbind(exons,cdss,introns,fiveUTRs,threeUTRs)
  anno.rna <- rnas

  return(list(anno.need=anno.need, anno.rna=anno.rna, anno.frame=anno.frame))
}


##funtion:parseGff
##parses genome annotation file of gff3/gtf format
parseGff <- function(aGFF) {

  if (!is.character(aGFF)) stop("aGFF not a character string!")
  if (!grepl('\\.gff3|\\.gtf', tolower(aGFF))) stop('aGFF not .gff3/.gff/.gtf!')

  if (grepl('\\.gff3|\\.gff', tolower(aGFF))) {
    #------------------------------------------------------
    #Loading annotation (gff3 format)
    #-------------------------------------------------------
    gff.path=aGFF
    anno <- import.gff3(gff.path)
    anno.frame <- as.data.frame(anno,stringsAsFactors =FALSE)
    anno.frame$seqnames <- as.character(anno.frame$seqnames)
    anno.frame$strand <- as.character(anno.frame$strand)
    anno.frame$type <- as.character(anno.frame$type)
    #print("###annotation file type information")
    #print(table(anno.frame$type))
    #delete chromosome information
    anno.frame$Parent <- sub(pattern="\\S+\\:",replacement = "",anno.frame$Parent)
    anno.frame$ID <- sub(pattern="\\S+\\:",replacement = "",anno.frame$ID)
    if(length(which(anno.frame$type=="chromosome"))){
      anno.frame <- anno.frame[-which(anno.frame$type=="chromosome"),]
    }
    #instead transcript to RNA
    anno.frame$type[which(anno.frame$type == "transcript")] <-"RNA"
    #getting RNA row
    rna.id <- grep("RNA$",anno.frame$type,ignore.case = FALSE)
    anno.rna <- anno.frame[rna.id,]

  } else if (grepl('\\.gtf', tolower(aGFF))) {
    gtf.path=aGFF
    anno <- import(gtf.path,format="gtf")
    anno.frame <- as.data.frame(anno,stringsAsFactors =FALSE)
    anno.frame$seqnames <- as.character(anno.frame$seqnames)
    anno.frame$strand <- as.character(anno.frame$strand)
    anno.frame$type <- as.character(anno.frame$type)
    anno.frame$Parent <- as.character(anno.frame$transcript_id)
    anno.frame$type[which(anno.frame$type == "transcript")] <-"RNA"
    #getting RNA row
    trans.id <- grep("transcript",anno.frame$type,ignore.case = FALSE)
    if(length(trans.id)){
      anno.frame$type[which(anno.frame$type == "transcript")] <-"RNA"
      rna.id <- grep("RNA$",anno.frame$type,ignore.case = FALSE)
    }else{
      rna.id <- grep("RNA$",anno.frame$type,ignore.case = FALSE)
    }
    if(length(rna.id)==0){
      anno.frame <- add_rna(anno.frame)
      rna.id <- grep("RNA$",anno.frame$type,ignore.case = FALSE)
    }
    if(length(which(anno.frame$type=="gene"))==0){
      anno.gene <- anno.frame[rna.id,]
      anno.gene$type<- "gene"
      anno.gene$Parent <- ""
      anno.frame <- rbind(anno.frame,anno.gene)
    }
    anno.frame$ID <- anno.frame$Parent
    if(length(which(anno.frame$type=="chromosome"))){
      anno.frame <- anno.frame[-which(anno.frame$type=="chromosome"),]
    }
    anno.rna <- anno.frame[rna.id,]
  } #~gtf


  #If the comment is incomplete, losing transcript_id
  if(!length(which(colnames(anno.rna) == "transcript_id"))){
    anno.rna$transcript_id<-anno.rna$ID
  }
  #table(anno.rna$type)

  #ID=transcript:AT1G01010.1;Parent=gene:AT1G01010;biotype=protein_coding;transcript_id=AT1G01010.1
  #1	araport11	five_prime_UTR	3631	3759	.	+	.	Parent=transcript:AT1G01010.1
  #confirm that the names of transcript is consistent with exon/cds/utr
  if(length(setdiff(anno.rna$transcript_id,anno.frame$Parent))){
    stop("Not consistent between transcript id in rna and exon/cds/utr")
  }
  #anno.frame$Parent[which(anno.frame$type=="gene")]<- ""
  # anno.frame.raw <- anno.frame
  if(is.na(match("three_prime_UTR",unique(anno.frame$type)))){
    if(is.na(match("CDS",unique(anno.frame$type)))){
      warning("This annotation without CDS, we can't identify UTR region")
    }else{
      print("Extracting UTR region")
      anno.frame <- add_utr(anno.frame)
    }
  }
  #=========================================================================
  #anno.need store cds/exon/5utr/3utr information
  anno.need <- anno.frame[which(anno.frame$Parent %in% anno.rna$transcript_id),]

  need.rna.id <- grep("RNA$",anno.need$type,ignore.case = FALSE)
  if(length(need.rna.id)){
    anno.need<-anno.need[-need.rna.id,]
  }

  index <- match(anno.need$Parent,anno.rna$transcript_id)

  if(length(which(is.na(index)))){
    stop("error can't find exon/cds/5utr/3utr 's parent")
  }
  anno.need$gene_id <- anno.rna$Parent[index]

  if(is.na(match("biotype",colnames(anno.rna)))){
    anno.rna$biotype <- NA
  }
  anno.need$biotype <- anno.rna$biotype[index]
  #====================================================================
  #ann.intron stores intron information
  exon.id <- grep("exon",anno.need$type,ignore.case = FALSE)
  ann.exon <- anno.need[exon.id,]
  if(length(which(is.na(ann.exon$Parent)))){
    print("exist some exon can't find parent id ")
  }
  ann.exon <- ann.exon[order(ann.exon$Parent,ann.exon$start,ann.exon$strand),]
  ann.exon.1 <- ann.exon[seq(1,nrow(ann.exon),2),]
  ann.exon.2 <- ann.exon[seq(2,nrow(ann.exon),2),]
  ann.exon.3 <- ann.exon[seq(3,nrow(ann.exon),2),]

  keep.num1 <- min(nrow(ann.exon.1),nrow(ann.exon.2))
  ann.exon.k1<-ann.exon.1[1:keep.num1,]
  ann.exon.k2<-ann.exon.2[1:keep.num1,]
  index <- which(ann.exon.k1$Parent == ann.exon.k2$Parent)
  if(!identical(ann.exon.k1$Parent[index],ann.exon.k2$Parent[index])){
    stop("something error with extart intron region")
  }
  ann.intron1 <- ann.exon.k1[index,]
  ann.intron1$type <- "intron"
  ann.intron1$start <- ann.exon.k1$end[index]+1
  ann.intron1$end <- ann.exon.k2$start[index]-1


  keep.num2 <- min(nrow(ann.exon.2),nrow(ann.exon.3))
  ann.exon.kk2<-ann.exon.2[1:keep.num2,]
  ann.exon.k3<-ann.exon.3[1:keep.num2,]
  index <- which(ann.exon.kk2$Parent == ann.exon.k3$Parent)
  if(!identical(ann.exon.kk2$Parent[index],ann.exon.k3$Parent[index])){
    stop("something error with extart intron region")
  }
  ann.intron2 <- ann.exon.kk2[index,]
  ann.intron2$type <- "intron"
  ann.intron2$start <- ann.exon.kk2$end[index]+1
  ann.intron2$end <- ann.exon.k3$start[index]-1
  ann.intron <- rbind(ann.intron1,ann.intron2)
  ann.intron <- ann.intron[order(ann.intron$Parent,ann.intron$start,ann.intron$strand),]
  anno.need <- rbind(anno.need,ann.intron)


  #table(anno.need$type)
  rna.error <- grep("RNA$",anno.need$type,ignore.case = FALSE)
  if(length(rna.error)){
    anno.need <- anno.need[-rna.error,]
  }
  return(list(anno.need=anno.need, anno.rna=anno.rna, anno.frame=anno.frame))
}



#=========================================================
#------------------------------------------------------
#function:add_utr()
#Adding 3UTR and 5UTR region
#--------------------------------------------------------
#======================================================
add_utr <- function(anno.frame=NULL){
  anno.cds <- anno.frame[which(anno.frame$type=="CDS"),]
  anno.exon <- anno.frame[which(anno.frame$type=="exon"),]
  rna.id <- grep("RNA$",anno.frame$type,ignore.case = FALSE)
  anno.rna <- anno.frame[rna.id,]
  if(!length(which(colnames(anno.rna) == "transcript_id"))){
    anno.rna$transcript_id<-anno.rna$ID
  }
  anno.cds.frist <- anno.cds[order(anno.cds$Parent,anno.cds$start,anno.cds$strand,decreasing = FALSE),]
  anno.cds.last <- anno.cds[order(anno.cds$Parent,anno.cds$start,anno.cds$strand,decreasing = TRUE),]
  anno.cds.frist <- anno.cds.frist[!duplicated(anno.cds.frist$Parent),]
  anno.cds.last <- anno.cds.last[!duplicated(anno.cds.last$Parent),]
  index.frist <-match(anno.cds.frist$Parent,anno.rna$transcript_id)
  index.last <-match(anno.cds.last$Parent,anno.rna$transcript_id)
  if(length(which(is.na(c(index.frist,index.last))))){
    stop("Can't find cds parent based on input annotation file ")
  }
  anno.cds.frist$utr.start <- anno.rna$start[index.frist]
  anno.cds.frist$utr.end <- anno.cds.frist$start -1
  anno.cds.frist <- anno.cds.frist[which( (anno.cds.frist$utr.end- anno.cds.frist$utr.start) >=0),]

  anno.cds.last$utr.start <- anno.cds.last$end +1
  anno.cds.last$utr.end <- anno.rna$end[index.last]
  anno.cds.last <- anno.cds.last[which((anno.cds.last$utr.end- anno.cds.last$utr.start) >=0),]


  gr.first <- GRanges(seqnames =as.character(anno.cds.frist$Parent) ,
                      ranges =IRanges(start=as.integer(anno.cds.frist$utr.start) ,
                                      end=as.integer(anno.cds.frist$utr.end)),
                      strand =as.character(anno.cds.frist$strand))

  gr.last <- GRanges(seqnames =as.character(anno.cds.last$Parent) ,
                     ranges =IRanges(start=as.integer(anno.cds.last$utr.start) ,
                                     end=as.integer(anno.cds.last$utr.end)),
                     strand =as.character(anno.cds.last$strand))

  gr.exon <- GRanges(seqnames =as.character(anno.exon$Parent) ,
                     ranges =IRanges(start=as.integer(anno.exon$start) ,
                                     end=as.integer(anno.exon$end)),
                     strand =as.character(anno.exon$strand))

  ov.first <- findOverlaps(gr.first,gr.exon)
  ov.last <- findOverlaps(gr.last,gr.exon)
  ov.first <- as.data.frame(ov.first)
  ov.last <- as.data.frame(ov.last)
  colnames(ov.first)<-c("cdsID","exonID")
  colnames(ov.last) <- c("cdsID","exonID")


  ov.first$utr.start <- as.integer(anno.cds.frist$utr.start[ov.first$cdsID])
  ov.first$utr.end <- as.integer(anno.cds.frist$utr.end[ov.first$cdsID])
  ov.first$exon.start <- as.integer(anno.exon$start[ov.first$exonID])
  ov.first$exon.end <- as.integer(anno.exon$end[ov.first$exonID])
  ov.first$utr.start.r <- ov.first$exon.start
  ov.first$utr.end.r <- apply(ov.first[,c("utr.end","exon.end")],1,min)
  five.utr <- anno.exon[ov.first$exonID,]
  five.utr$start <- ov.first$utr.start.r
  five.utr$end <- ov.first$utr.end.r
  if(nrow(five.utr)){
    five.utr$type <- "five_prime_UTR"
    five.utr$type[which(five.utr$strand=="-")] <- "three_prime_UTR"
  }



  ov.last$utr.start <- as.integer(anno.cds.last$utr.start[ov.last$cdsID])
  ov.last$utr.end <- as.integer(anno.cds.last$utr.end[ov.last$cdsID])
  ov.last$exon.start <- as.integer(anno.exon$start[ov.last$exonID])
  ov.last$exon.end <- as.integer(anno.exon$end[ov.last$exonID])
  ov.last$utr.start.r <- apply(ov.last[,c("utr.start","exon.start")],1,max)
  ov.last$utr.end.r <- ov.last$exon.end
  three.utr <- anno.exon[ov.last$exonID,]
  three.utr$start <- ov.last$utr.start.r
  three.utr$end <- ov.last$utr.end.r
  if(nrow(three.utr)){
    three.utr$type <- "three_prime_UTR"
    three.utr$type[which(three.utr$strand=="-")] <- "five_prime_UTR"
  }
  utr <- rbind(three.utr,five.utr)
  utr <- utr[order(utr$Parent,utr$type,utr$start),]
  utr$width <- as.integer(utr$end-utr$start+1)

  #-------------------------------------
  #check result
  #  really.utr <- anno.frame[which(anno.frame$type %in% c("three_prime_UTR","five_prime_UTR")),]
  #  really.utr <- really.utr[order(really.utr$Parent,really.utr$type,really.utr$start),]
  # length(unique(really.utr$Parent))
  # length(unique(utr$Parent))
  # identical(utr$start,really.utr$start)
  # identical(utr$end,really.utr$end)
  # identical(utr$strand,really.utr$strand)
  #  write.table(really.utr,file="really_utr.txt",col.names = TRUE,row.names = FALSE,sep="\t",
  #              quote=FALSE)
  #  write.table(utr,file="build_utr.txt",col.names = TRUE,row.names = FALSE,sep="\t",
  #              quote=FALSE)
  anno.frame <-rbind(anno.frame,utr)
  return(anno.frame)
}

#=========================================================
#------------------------------------------------------
#function:add_rna()
#Adding RNA region
#--------------------------------------------------------
#======================================================
add_rna <- function(anno.frame=NULL){
  anno.exon <- anno.frame[which(anno.frame$type=="exon"),]
  anno.exon.order <- anno.exon[order(anno.exon$gene_id,anno.exon$transcript_id,
                                     anno.exon$strand,anno.exon$start,decreasing = FALSE),]
  anno.exon.rev <- anno.exon.order[nrow(anno.exon.order):1,]

  anno.exon.order.unique <- anno.exon.order[!duplicated(anno.exon.order$transcript_id),]
  anno.exon.rev.order <- anno.exon.rev[!duplicated(anno.exon.rev$transcript_id),]
  anno.rna <- anno.exon.order.unique
  index <- match(anno.rna$transcript_id,anno.exon.rev.order$transcript_id)
  anno.rna$end <- anno.exon.rev.order$end[index]
  anno.rna$Parent <- anno.rna$gene_id
  anno.rna$type <- "mRNA"
  anno.frame <-rbind(anno.frame,anno.rna)
  return(anno.frame)
}





# -------------- *** readBAM() **** ---------------
#' Read the BAM file
#'
#' @usage readBAM(bam.path=NULL,bams=NULL,set.n=NULL,bam.keyword="*sorted.bam$",file=NULL,label=cell.type.name,group=rep(1,length(cell.type.name)))
#' @param bam.path path of input of bam
#' @param bams bam files
#' @param set.n select the number of  bam files
#' @param bam.keyword selected bam files
#' @return bam file
#' @name readBAM
#' @export
readBAM<-function(bam.path=NULL,bams=NULL,set.n=NULL,bam.keyword="*sorted.bam$",file=NULL,label=cell.type.name,group=rep(1,length(cell.type.name))){
  if(is.null(bams)==FALSE){

    bam.file<-list(fileName=bams[[1]][set.n],group=bams[[2]][set.n],label=bams[[3]][set.n],bam.path=bams[[4]])
  }
  if(!is.null(file)){
    bam<-as.list(file)
    bam.file<-list(fileName=bam$fileName,group=bam$group,label=bam$label,bam.path=bam.path)
  }else if
  (!is.null(bam.path)& is.null(file)){
    if(length(grep("\\/$",bam.path))==0){
      bam.path <- paste0(bam.path,"/")
    }
    setwd(bam.path)
    bam.files <- list.files("./",bam.keyword)

    cell.type.name <- as.character(sapply(strsplit(bam.files,"\\_"),
                                          "[[",3))
    cell.type.name <- gsub("\\.sorted.bam","",cell.type.name)
    n=length(cell.type.name)
    bam.file<-list(fileName=bam.files,group=group,label=label,bam.path=bam.path)
    if(is.null(set.n)==FALSE){ bam.file<-list(fileName=bam.file[[1]][set.n],group=bam.file[[2]][set.n],label=bam.file[[3]][set.n],bam.path=bam.file[[4]])
    }
  }
  return(bam.file)
}
# -------------- *** isChrConsistent() **** ---------------

#' Verify chromosome names
#'
#' @usage isChrConsistent(pacds, obj, col=NULL, allin=FALSE)
#' @param pacds a PACdataset with @anno$chr, or df/matrix with chr column, or vector
#' @param obj BSgenome, or FaFile, or data.frame/matrix or vector
#' @param col the colname of chr in obj[if obj is df or matrix]
#' @param allin pacds$chr should all in obj's chr, otherwise intersect is allowed
#' @return TRUE (consistent or intersect) or FALSE
#' @examples
#' isChrConsistent(PACds, obj=bsgenome, allin=TRUE)
#' isChrConsistent(PACds, obj=FaFile(aFastaFile), allin=TRUE)
#' isChrConsistent(PACds, obj=unique(PACds@anno$chr), allin=TRUE)
#' isChrConsistent(PACds, obj=c('x','a'), allin=TRUE)
#' isChrConsistent(PACds, obj=d, allin=TRUE)
#' isChrConsistent(PACds, obj=d, col='strand',allin=TRUE)
#' isChrConsistent(c('1','2'), obj=c('chr1','chr2'),allin=TRUE)
#' @name isChrConsistent
#' @export

isChrConsistent<-function(pacds, obj, col=NULL, allin=FALSE) {
  if (is.vector(pacds)) {
    pacchr=pacds
  } else if  (class(pacds)=='data.frame' | class(obj)=='matrix') {
    pacchr=unique(as.character(pacds$chr))
  } else if (class(pacds)=='PACdataset' ) {
    pacchr=unique(pacds@anno$chr)
  }

  if (class(obj)=='BSgenome') {
    objchr=seqnames(bsgenome)
  } else if (class(obj)=='FaFile') {
    if (!file.exists(obj$index)) {
      cat("Indexing fa...\n")
      indexFa(obj$path)
    }
    objchr=names(seqlengths(obj))
  } else if (class(obj)=='data.frame' | class(obj)=='matrix') {
    if (is.null(col)) {
      if ('seqnames' %in% colnames(obj)) {
        col='seqnames'
      } else if ('chr' %in% colnames(obj)) {
        col='chr'
      }
    }
    if (!(col %in% colnames(obj))) {
      stop(cat('chrCol',col, 'not in obj\n'))
    }
    objchr=unique(obj[,col])
  } else {
    objchr=obj
  }
  if (allin) {
    return(AinB(pacchr, objchr, all=TRUE))
  } else {
    return(length(intersect(pacchr,objchr))!=0)
  }
}


# -------------- *** vizTracks() **** ---------------

library(Gviz)
library(GenomicRanges)

library(biomaRt)
library(GenomicFeatures)
library(rtracklayer)
library(Rsamtools)
library(bamsignals)
library(stringr)

#' track plot
#'
#' show the gene model of a gene or a chromosomal region,the PAS model and the coverage  model of a specified PAS across different cell types
#'
#' @usage vizTracks(gtf.data,bams,pas.postition,region,gene,PA, space5=300,space3=300,
#'                  showGene=TRUE, showallPA=FALSE,bm=bm,Group=FALSE,genome.version=NULL,
#'                  cov.Merge="No",showCov="Single",covStyle="cov",
#'                  bams.col="common",bams.fill="common",bams.background.col="common",
#'                  PAcol="red",PAfill="red",
#'                  cov.Merge.background.col="red",cov.Merge.col="red",cov.Merge.fill="red")
#' @param gtf.data annotation gtf fotmat  file
#' @param bams bam file,it is used to extract the expression level of PAS to plot the PAS coverage information
#' @param pas.position poly(A) site information
#' @param region specify the chromosomal region
#' @param gene specify the gene
#' @param pA specify the poly (A) site
#' @param space5 extend 5' end,default is 300 bp
#' @param space3 extend 3' end,default is 300 bp
#' @param showGene If TRUE,tnen the area where the map is drawn becomes the gene area where the PA is located
#' If FALSE,then only show that part of the area where the PA is located
#' @param showallPA If TRUE,then  show all PAs.Otherwise Only show the PA specified under this gene
#' @param bm BioMart database
#' @param Group default is FALSE.If TRUE,then plot coverage models on one track.If FALSE, then each condition is one track.
#' @param genome.version genome.version
#' @param cov.Merge Avg/Sum/No.If cov.Merge=Avg,then  average the coverage models on one track .
#' If cov.Merge=Sum,then sum the coverage models on one track.If cov.Merge=NO,then make a separate track.
#' @param showCov Both/Single/Merge.default is Single.If showCov=Both and cov.Merge=Avg/Sum,then both the merged  and the separate coverage model are showed.
#' If showCov=Single,then only show the separate coverage model.
#' If showCov=Merge and cov.Merge=Avg/Sum,then only show the merged coverage model.
#' @param covStyle cov/curve.If covStyle=cov,then  the coverage model draws in the form of an overlay.
#' If covStyle=curve ,the coverage model draws in the form of a curve.
#' @return vizTracks() reutrns Gviz's Tracks.If gtf.data is provided,then plot the gene model.
#' If pas.position is provided,then piot the pA model.
#' If bam file is provided,then plot the coverage model across different cell types.
#' @examples
#' #way1.specify the chromosomal region
#' vizTracks(gtf.data=gtf.data,bams=bams,pas.postition=pas.postition,region="chr1:+:3000:9000", bm=bm，cov.Merge = "Avg",showCov = "Both")
#' #way2.specify the gene
#' vizTracks(gtf.data = gtf.data,pas.position = pas.postition,bams=bams,gene="AT1G01040",bm=bm,cov.Merge = "Avg",showCov = "Both")
#' #way3.specify the poly (A) site
#' vizTracks(gtf.data=gtf.data,bams=bams,pas.position = pas.position,PA="chr1:+:10000447",bm=bm,
#'           genome.version = genome.version,showGene = FALSE,showallPA = TRUE,cov.Merge = "Avg",showCov = "Both")
#' @name vizTracks
#' @export
vizTracks<-function(gtf.data,bams,pas.position,region=NULL,gene=NULL,PA=NULL, space5=300,space3=300,
                    showGene=TRUE, showallPA=TRUE,Group=FALSE,bm,genome.version=NULL,
                    cov.Merge="No",showCov="Single",covStyle="cov",
                    bams.col="common",bams.fill="common",bams.background.col="common",
                    PAcol="red",PAfill="red",
                    cov.Merge.background.col="red",cov.Merge.col="red",cov.Merge.fill="red"){
  options(stringsAsFactors = FALSE)
  options(ucscChromosomeNames=FALSE)
  gtf.data <- as.data.frame(gtf.data)
  gtf.data <- subset(gtf.data,type=="gene")
  gtf.d<-gtf.data
  bs<-bams
  pas.p<-pas.position
  s5<-space5
  s3<-space3
  se<-showGene

  Gp<-Group
  genome.v<-genome.version
  cov.M<-cov.Merge
  b<-bm
  showC<-showCov
  covS<-covStyle
  bams.c<-bams.col
  bams.f<-bams.fill
  bams.b<-bams.background.col
  PAc<- PAcol
  PAf<-PAfill
  cov.Merge.background.c<-cov.Merge.background.col
  cov.Merge.c<-cov.Merge.col
  cov.Merge.f<-cov.Merge.fill

  #####when provide the region

  if(!is.null(region)){
    region<-strsplit(region,":",fixed=T)
    chromosome<-strsplit(region[[1]][1],"chr")[[1]][2]
    strand<-region[[1]][2]
    region.from<-as.numeric(region[[1]][3])
    region.to<-as.numeric(region[[1]][4])
    if(region.to-region.from>10000){
      print("the region is so big,please set again")
    }
    vizTracksOfRegion(gtf.data=gtf.d,bams=bs,pas.position=pas.p,chromosome=chromosome,strand=strand,
                      start=region.from,end=region.to,PAcoord=NULL,space5=s5,space3=s3,
                      showGene=se, showallPA=TRUE,Group = Gp,bm=b,genome.version=genome.v,
                      cov.Merge=cov.M,showCov=showC,covStyle=covS,
                      bams.col=bams.c,bams.fill=bams.f,bams.background.col=bams.b,
                      PAcol= PAc,PAfill=PAf,
                      cov.Merge.background.col=cov.Merge.background.c,cov.Merge.col=cov.Merge.c,cov.Merge.fill=cov.Merge.f)
  }

  #####when provide the gene
  if(!is.null(gene)){

    if(length(gene%in%gtf.data$gene_id)==1){
      n<-which(gtf.data$gene_id==gene)
    }else if(length(gene%in%gtf.data$gene_name)==1){
      n<-which(gtf.data$gene_name==gene_name)
    }else{
      print("warnings:there is no gene")
    }
    chromosome <-as.character(gtf.data$seqnames[n])
    strand <-as.character(gtf.data$strand[n])
    region.from<-as.numeric(gtf.data$start[n])
    region.to<-as.numeric(gtf.data$end[n])

    vizTracksOfRegion(gtf.data=gtf.d,bams=bs,pas.position=pas.p,chromosome=chromosome,strand=strand,
                      start=region.from,end=region.to,PAcoord=NULL,space5=s5,space3=s3,
                      showGene=se, showallPA=TRUE,Group = Gp,bm=b,genome.version=genome.v,
                      cov.Merge=cov.M,showCov=showC,covStyle=covS,
                      bams.col=bams.c,bams.fill=bams.f,bams.background.col=bams.b,
                      PAcol= PAc,PAfill=PAf,
                      cov.Merge.background.col=cov.Merge.background.c,cov.Merge.col=cov.Merge.c,cov.Merge.fill=cov.Merge.f)
  }

  ######when provide the PA
  if(!is.null(PA)){
    if(PA%in%rownames(pas.position==TRUE)){
      pas.n<-which(row.names(pas.position)==PA)
      chromosome<-as.character(pas.position[pas.n,]$chr)
      strand<-as.character(pas.position[pas.n,]$strand)
      coord<-pas.position[pas.n,]$coord
    }else{
      PA<-strsplit(PA,":",fixed=T)
      chromosome<-strsplit(PA[[1]][1],"chr")[[1]][2]
      strand<-PA[[1]][2]
      coord<-as.numeric(PA[[1]][3])
      pas.n1=which(pas.position$chr==chromosome)
      pas.n2=which(pas.position$strand==strand)
      pas.n3=which(pas.position$coord==coord)
      pas.n=intersect(x=pas.n3,y=intersect(x=pas.n1,y=pas.n2))

    }

    sA<-showallPA

    pas.position.sub<-pas.position[pas.n,]
    gene<-pas.position.sub$gene
    n<-which(gtf.data$gene_id==gene)
    region.from<-gtf.data$start[n]
    region.to<-gtf.data$end[n]

    vizTracksOfRegion(gtf.data=gtf.d,bams=bs,pas.position=pas.p,chromosome=chromosome,strand=strand,
                      start=region.from,end=region.to,PAcoord=coord,space5=s5,space3=s3,
                      showGene=se, showallPA=sA,Group=Gp,bm=b,genome.version=genome.v,
                      cov.Merge=cov.M,showCov=showC,covStyle=covS,
                      bams.col=bams.c,bams.fill=bams.f,bams.background.col=bams.b,
                      PAcol= PAc,PAfill=PAf,
                      cov.Merge.background.col=cov.Merge.background.c,cov.Merge.col=cov.Merge.c,cov.Merge.fill=cov.Merge.f)
  }
}


##################################vizTracksofRegion()
vizTracksOfRegion<-function(gtf.data,bams,pas.position,chromosome,strand,start, end,PAcoord=NULL,space5=300,space3=300,
                            showGene=TRUE, showallPA=TRUE,Group=FALSE,bm,genome.version=NULL,
                            cov.Merge="No",showCov="Single",covStyle="cov",
                            bams.col="common",bams.fill="common",bams.background.col="common",
                            PAcol="red",PAfill="red",
                            cov.Merge.background.col="red",cov.Merge.col="red",cov.Merge.fill="red"){
  options(stringsAsFactors = FALSE)
  options(ucscChromosomeNames=FALSE)
  colors.pattern <- c("red","orange","yellow","green","cyan", "blue", "purple","black",
                      "aqua","blueviolet","brown","chocolate","coral","darkgray","darkgreen",
                      "forestgreen","gold","gray","lightgrey","pink")
  start<-as.numeric(start)
  end<-as.numeric(end)

  ifelse(strand=="+",c( region.from<-start-space5,region.to<-end+space3),
         c(region.from<-start-space3, region.to<-end+space5))


  ########Building Gene Model

  biomTrack <- BiomartGeneRegionTrack(genome =genome.version,
                                      name = "Gene model",
                                      chromosome = chromosome,
                                      start=region.from,
                                      end=region.to,
                                      biomart = bm
  )
  displayPars(biomTrack) <- list(col ="#595959",fill="#595959",
                                 background.title = "brown",
                                 fontfamily="Arial",
                                 fontfamily.title="Arial",
                                 col.grid="black",
                                 fontfamily.group="Arial",
                                 utr5="red",
                                 utr3="black",
                                 protein_coding="#595959",
                                 col.line="black")

  ########Adding genome Axis
  gtrack <- GenomeAxisTrack()
  ########Adding PAS sites
  gene_region <- GRanges(
    seqnames=chromosome,
    ranges = IRanges(start = region.from,end = region.to),
    strand = strand
  )
  pas.region <- makeGRangesFromDataFrame(pas.position,keep.extra.columns = TRUE)
  if(showallPA==TRUE){
    ov.region <- findOverlaps(pas.region,gene_region)
    pas.position.sub <- pas.position[ov.region@from,]

    if(nrow(pas.position.sub)==0){
      print("No poly(A) sites in this gene")
      next;
    }

    if(strand=="-"){
      start<-pas.position.sub$start-space3
      end<-pas.position.sub$end+space5
    }else{
      start<-pas.position.sub$start-space5
      end<- pas.position.sub$end+space3
    }
    aTrack <- AnnotationTrack(start=start,
                              end=end,
                              chromosome=pas.position.sub$chr,
                              strand=pas.position.sub$strand,
                              id=as.character(pas.position.sub$peakID),
                              genome=biomTrack@genome, name="PAS",
                              background.title = "brown",
                              fontfamily="Arial",
                              fontfamily.title="Arial",
                              shape="box",
                              col =PAcol,fill=PAfill,fontsize.legend=14
    )
  }else{pas.n1<-which(pas.position$chr==chromosome)
  pas.n2<-which(pas.position$strand==strand)
  pas.n3<-which(pas.position$coord==PAcoord)
  pas.n<-intersect(x=pas.n3,y=intersect(x=pas.n1,y=pas.n2))

  pas.position.sub<-pas.position[pas.n,]


  if(strand=="-"){
    start<-pas.position.sub$start-space3
    end<-pas.position.sub$end+space5
  }else{start<-pas.position.sub$start-space5
  end<-pas.position.sub$end+space3}

  aTrack <- AnnotationTrack(start=start,
                            end=end,
                            chromosome=pas.position.sub$chr,
                            strand=pas.position.sub$strand,
                            id=as.character(pas.position.sub$peakID),
                            genome=biomTrack@genome, name="PAS",
                            background.title = "brown",
                            fontfamily="Arial",
                            fontfamily.title="Arial",
                            shape="box",
                            col =PAcol,fill=PAfill,fontsize.legend=14 )
  }# if showallPA=F



  ht.aTrack <- HighlightTrack(trackList = list(aTrack ),
                              start=pas.position.sub$coord,width=3,
                              chromosome =chromosome)

  if(showGene==TRUE){
    gviz.list <- list(gtrack,biomTrack,aTrack,ht.aTrack)
    names(gviz.list) <- c("genome_scale","gene_model",
                          "pas","pas_coord")
    gviz.n<-length(gviz.list)-1
    plot.n<-gviz.n+2
  }

  if(showGene==FALSE){
    gviz.list <- list(gtrack,aTrack,ht.aTrack)
    names(gviz.list) <- c("genome_scale",
                          "pas","pas_coord")
    gviz.n=length(gviz.list)-1
    plot.n=gviz.n+2
  }

  if(strand=="+"){
    bam.n<-which(str_detect(bams$fileName,"forward.sorted_"))
    bams$fileName<- bams$fileName[grep("forward.sorted_",bams$fileName)]
    bams<-list(fileName=bams$fileName,group=bams$group[bam.n],label=bams$label[bam.n],bam.path=bams$bam.path)
  }else if(strand=="-"){
    bam.n<-which(str_detect(bams$fileName,"reverse.sorted_"))
    bams$fileName<- bams$fileName[grep("reverse.sorted_",bams$fileName)]
    bams<-list(fileName=bams$fileName,group=bams$group[bam.n],label=bams$label[bam.n],bam.path=bams$bam.path)
  }else {
    print("warnings: don't find any bam file with forward_sorted or reverse.sorted")
  }
  if(Group==TRUE){
    genes <- gene_region
    coverage.data<- data.frame(start=c(start(genes):end(genes)))
    coverage.data$seqnames <-as.character(genes@seqnames)
    coverage.data$strand <- as.character(genes@strand)
    coverage.data$end <- coverage.data$start

    for(i in 1:length(bams$fileName)){
      bampath<-paste0(bams$bam.path,bams$fileName[i])
      bf<-Rsamtools::BamFile(bampath)
      covSigs <- bamCoverage(bampath, genes, verbose=FALSE)
      middle<- data.frame(covSigs[1])
      colnames(middle) <- bams$label[i]
      coverage.data <- cbind(coverage.data,middle)
    }
    if(strand=="-"){
      coverage.data$start<-rev(coverage.data$start)
      coverage.data$end<-rev(coverage.data$end)
    }
    coverage.region<-makeGRangesFromDataFrame(coverage.data,keep.extra.columns = TRUE)
    coverage.track<-DataTrack(coverage.region,name="Coverage")
    displayPars(coverage.track) <- list(
      background.title="#2e4057",groups=bams$label,
      col=colors.pattern[1:length(bams$label)],
      legend=TRUE,cex.legend=0.8,fontsize.legend=12
    )
    if(showGene==TRUE){
      plotTracks(list(gviz.list[[1]],
                      gviz.list[[2]],
                      gviz.list[[3]],coverage.track),

                 type=c("l"),window=-1,
                 from=region.from,
                 to=region.to)}
    if(showGene==FALSE){
      plotTracks(list(gviz.list[[1]],
                      gviz.list[[2]],
                      coverage.track),

                 type=c("l"),window=-1,
                 from=region.from,
                 to=region.to)}


  }else{
    for(i in 1:length(bams$fileName)){
      if(covStyle=="cov"){
        type="h"
      }else if(covStyle=="curve"){
        type="l"
      }
      data.track.test <- DataTrack(range=paste0(bams$bam.path,bams$fileName[i]),
                                   genome=genome.version,
                                   type=type,
                                   name=bams$label[i],
                                   window=-1,
                                   chromosome=chromosome,
                                   from=region.from,
                                   to=region.to
      )
      if(bams.col%in%colors.pattern==FALSE){
        bam.col=colors.pattern[i]
      }
      if(bams.fill%in%colors.pattern==FALSE){
        bam.fill=colors.pattern[i]
      }
      if(bams.background.col%in%colors.pattern==FALSE){
        bam.background.col=colors.pattern[i]
      }
      displayPars(data.track.test)<-list(
        col=bam.col,fill=bam.fill,
        background.title =bam.background.col
      )
      gviz.list[[bams$label[i]]] <- data.track.test
    }

    if(cov.Merge=="Avg"){
      genes <- gene_region
      coverage.data<- data.frame(start=c(start(genes):end(genes)))
      coverage.data$seqnames <-as.character(genes@seqnames)
      coverage.data$strand <- as.character(genes@strand)
      coverage.data$end <- coverage.data$start

      for(i in 1:length(bams$fileName)){
        bampath<-paste0(bams$bam.path,bams$fileName[i])
        bf<-Rsamtools::BamFile(bampath)
        covSigs <- bamCoverage(bampath, genes, verbose=FALSE)
        middle<- data.frame(value=covSigs[1])
        colnames(middle) <- bams$label[i]
        coverage.data <- cbind(coverage.data,middle)
      }
      if(strand=="-"){
        coverage.data$start<-rev(coverage.data$start)
        coverage.data$end<-rev(coverage.data$end)
      }
      coverage.region <- makeGRangesFromDataFrame(coverage.data,
                                                  keep.extra.columns = TRUE)
      coverage.track <- DataTrack(coverage.region,name="Coverage")
      if(covStyle=="cov"){
        type=c("a","h")
        displayPars(coverage.track)<-list(
          background.title=cov.Merge.background.col,
          col=cov.Merge.col,fill=cov.Merge.fill,type=c("a","h"),
          legend=TRUE,cex.legend=0.8,fontsize.legend=12
        )
      }else if(covStyle=="curve"){

        displayPars(coverage.track)<-list(
          background.title=cov.Merge.background.col,
          col=cov.Merge.col,type="l",fill=cov.Merge.fill,
          legend=TRUE,cex.legend=0.8,fontsize.legend=12
        )
      }
      if(showCov=="Merge"){
        if(showGene==TRUE){
          plotTracks(list(gviz.list[[1]],
                          gviz.list[[2]],
                          gviz.list[[3]],
                          coverage.track),window=-1,
                     transcriptAnnotation="symbol",
                     from=region.from,
                     to=region.to)}
        if(showGene==FALSE){
          plotTracks(list(gviz.list[[1]],
                          gviz.list[[2]],
                          coverage.track),window=-1,
                     transcriptAnnotation="symbol",
                     from=region.from,
                     to=region.to)}

      }else if(showCov=="Both"){
        gviz.list[["coverage"]]<-coverage.track
        plotTracks((gviz.list[c(1:gviz.n,plot.n:length(gviz.list))]),window=-1,
                   transcriptAnnotation="symbol",
                   from=region.from,
                   to=region.to)
      }
    }else if(cov.Merge=="Sum"){
      genes <- gene_region
      coverage.data<- data.frame(start=c(start(genes):end(genes)))
      coverage.data$seqnames <-as.character(genes@seqnames)
      coverage.data$strand <- as.character(genes@strand)
      coverage.data$end <- coverage.data$start

      for(i in 1:length(bams$fileName)){
        bampath<-paste0(bams$bam.path,bams$fileName[i])
        bf<-Rsamtools::BamFile(bampath)
        covSigs <- bamCoverage(bampath, genes, verbose=FALSE)
        middle<- data.frame(value=covSigs[1])
        colnames(middle) <- bams$label[i]
        coverage.data <- cbind(coverage.data,middle)
      }
      if(strand=="-"){
        coverage.data$start<-rev(coverage.data$start)
        coverage.data$end<-rev(coverage.data$end)
      }
      coverage.region<-makeGRangesFromDataFrame(coverage.data,keep.extra.columns = TRUE)
      coverage.track<-DataTrack(coverage.region,name="Coverage")
      if(covStyle=="cov"){

        displayPars(coverage.track)<-list(
          background.title=cov.Merge.background.col,type="h",
          col=cov.Merge.col,fill=cov.Merge.fill,
          legend=TRUE,cex.legend=0.8,fontsize.legend=12
        )
      }else if(covStyle=="curve"){

        displayPars(coverage.track)<-list(
          background.title=cov.Merge.background.col,type="l",
          col=cov.Merge.col,fill=cov.Merge.fill,
          legend=TRUE,cex.legend=0.8,fontsize.legend=12
        )
      }
      if(showCov=="Merge"){
        if(showGene==TRUE){
          plotTracks(list(gviz.list[[1]],
                          gviz.list[[2]],
                          gviz.list[[3]],
                          coverage.track),window=-1,
                     transcriptAnnotation="symbol",
                     from=region.from,
                     to=region.to)}
        if(showGene==FALSE){
          plotTracks(list(gviz.list[[1]],
                          gviz.list[[2]],
                          coverage.track),window=-1,
                     transcriptAnnotation="symbol",
                     from=region.from,
                     to=region.to)}
      }else if(showCov=="Both"){
        gviz.list[["coverage"]]<-coverage.track
        plotTracks((gviz.list[c(1:gviz.n,plot.n:length(gviz.list))]),window=-1,
                   transcriptAnnotation="symbol",
                   from=region.from,
                   to=region.to)
      }
    }else if(cov.Merge=="No"){
      if(showCov=="Single"){
        plotTracks((gviz.list[c(1:gviz.n,plot.n:length(gviz.list))]),
                   transcriptAnnotation="symbol",from = region.from,to=region.to,
                   window=-1)
      }
    }
  }
}

########################################vizPAplot

# -------------- *** vizPAplot() **** ---------------
#' box/violin/point plot
#'
#'Draws a box/violin/point plot to show all PA expression data(PA count/PA ratio) in a chromosomal region or gene across different cell types
#'
#' @usage vizPAplot(scPACds,datatype,group,object="chr",type="point",log=TRUE,cols=NULL)
#' @param scPACds A dataset of the PACdataset class
#' @param datatype count/ratio.when datatype=count,the PA data provided is PA count;when datatype=ratio,the PA data provided is PA ratio.
#' @param group the column name containing cell classification information, stored in colData
#' @param object chr/gene.Specify whether the drawing object is chromosome or genetic PA sites
#' @param type box/violin/point:specified drawing type .If type=point and object=chr,show it in the form of box plot, violin plot plus point plot;
#' if  type=point and object=gene, show it with point plot
#' @param log Plot the y axis value on log scale(calculate the natural logarithm of the raw data)
#' @param cols Colors to use for plotting.When cols=NULL,use the default color matching
#' @return vizPAplot() reutrns box/violin/point plot
#' @examples
#' #way1.plot the chromosomal region
#' allPASC<-subsetPACds(scPACds,group = "celltype",chrs="chr12",noIntergenic=TRUE,pool=TRUE)
#' vizPAplot(scPACds =allPASC,datatype="count",group="celltype",object="chr",type="point",log=TRUE)
#' #way2.plot the gene
#' PAgene<-subsetPACds(scPACds,group = "celltype",noIntergenic=TRUE,genes="ENSMUSG00000021180")
#' vizPAplot(scPACds =PAgene,datatype="count",group="celltype",object="gene",type="box")
#' @name vizPAplot
#' @export
vizPAplot<-function(scPACds,datatype,group,object="chr",type="point",log=TRUE,cols=NULL){
  cbPalette<-c("#FF6666","#CC79A7","#00CC99","#006633","#99CCCC","#FFCC00","#993333","#FF6633","#CC3333",
               "#999999","#666666","#99CC99","#000033","#FF0033","#3366CC","#FF3300","#6be64d","#e66066", "#CC0099","#0000FF")
  if(object=="chr"){
    if(ncol(scPACds@counts)>=16){
      message("Warning:There are more than 16 cell types,too many cell types may cause the drawing graphics to be not beautiful enough,please choose the cell type through subsetPACds()")
    }

    allPA<-reshape2::melt(scPACds@counts,measure.vars=colnames(scPACds@counts),variable.name="celltype",value.name = "expr")
    if(datatype=="count"){
      ylabs<-"count"
      if(log==TRUE){
        allPA$expr <- log(allPA$expr)
        ylabs<-"log(count)"
      }

    }
    else{ylabs<-"ratio"}

    if(type=="box")
    {p<-ggplot(allPA, aes(x=celltype,y=expr,fill=celltype)) +
      geom_boxplot(alpha=0.8, outlier.colour="#FF3333",
                   outlier.shape=8,
                   outlier.size=2
      ) +  #使用outlier标注异常值
      stat_boxplot(geom = "errorbar",
                   lwd=0.5,
                   width=0.2)+ #添加误差线
      labs(y=ylabs ,fill = "Cell Type")+# fill为修改图例标题
      scale_fill_manual(values=cbPalette)+
      theme_classic()  #将背景修改为白色且没有灰色网格线

    if(!is.null(cols)){
      p<-ggplot(allPA, aes(x=celltype,y=expr,fill=celltype)) +
        geom_boxplot(alpha=0.8, outlier.colour="#FF3333",
                     outlier.shape=8,
                     outlier.size=2
        ) +  #使用outlier标注异常值
        stat_boxplot(geom = "errorbar",
                     lwd=0.5,
                     width=0.2)+ #添加误差线
        labs(y=ylabs ,fill = "Cell Type")+# fill为修改图例标题
        scale_fill_manual(values=cols)+
        theme_classic()
    }
    }
    if(type=="violin")
    {p<-ggplot(allPA, aes(x=celltype,y=expr)) +
      geom_violin(aes(fill=celltype),alpha=0.8 ) +
      stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                   geom="pointrange", color = "red")+  #加均值及方差
      labs(y=ylabs ,fill = "Cell Type")+
      scale_fill_manual(values=cbPalette)+
      theme_classic()
    if(!is.null(cols)){
      p<-ggplot(allPA, aes(x=celltype,y=expr)) +
        geom_violin(aes(fill=celltype),alpha=0.8 ) +
        stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
                     geom="pointrange", color = "red")+  #加均值及方差
        labs(y=ylabs ,fill = "Cell Type")+
        scale_fill_manual(values=cols)+
        theme_classic()
    }
    }
    if(type=="point")
    {p<-ggplot(allPA, aes(x=celltype,y=expr,color=celltype)) +#color是边框的颜色，fill是图形填充的颜色
      geom_jitter(alpha=0.3,
                  position=position_jitterdodge(jitter.width=0.35,
                                                jitter.height=0,
                                                dodge.width=0.8))+
      geom_boxplot(alpha=0.3,width=0.45,
                   position=position_dodge(width=0.8),
                   linewidth=0.6,outlier.colour=NA)+
      geom_violin(alpha=0.3,width=0.9,
                  position=position_dodge(width=0.8),
                  linewidth=0.6)+
      scale_color_manual(values=cbPalette)+
      labs(y=ylabs ,color = "Cell Type")+
      theme_classic()
    if(!is.null(cols)){
      p<-ggplot(allPA, aes(x=celltype,y=expr,color=celltype)) +#color是边框的颜色，fill是图形填充的颜色
        geom_jitter(alpha=0.3,
                    position=position_jitterdodge(jitter.width=0.35,
                                                  jitter.height=0,
                                                  dodge.width=0.8))+
        geom_boxplot(alpha=0.3,width=0.45,
                     position=position_dodge(width=0.8),
                     linewidth=0.6,outlier.colour=NA)+
        geom_violin(alpha=0.3,width=0.9,
                    position=position_dodge(width=0.8),
                    linewidth=0.6)+
        scale_color_manual(values=cols)+
        labs(y=ylabs ,color = "Cell Type")+
        theme_classic()
    }
    }
    return(p)

  }
  if(object=="gene"){

    if(nrow(scPACds@counts)>=2){
      PA1<-scPACds@counts[1,]
      PA1<-as.data.frame(t(PA1))
      PA1["celltype"]<-scPACds@colData[group]#添加细胞类型这一列
      PA1<-PA1%>%mutate(PA=rownames(scPACds@counts[1,]))#添加PA类型这一列
      colnames(PA1)<-c("expr","celltype","PA")

      for(i in 2:nrow(scPACds@counts)){
        PAi<-scPACds@counts[i,]
        PAi<-as.data.frame(t(PAi))
        celltype<-group
        PAi["celltype"]<-scPACds@colData[group]#添加细胞类型这一列
        PAi<-PAi%>%mutate(PA=rownames(scPACds@counts[i,]))#添加PA类型这一列
        colnames(PAi)<-c("expr","celltype","PA")
        PA1<-rbind(PA1,PAi)
      }
    }else{PA1<-scPACds@counts[1,]
    PA1<-as.data.frame(t(PA1))
    PA1["celltype"]<-scPACds@colData[group]#添加细胞类型这一列
    PA1<-PA1%>%mutate(PA=rownames(scPACds@counts[1,]))#添加PA类型这一列
    colnames(PA1)<-c("expr","celltype","PA")
    }
    if(datatype=="count"){
      ylabs<-"count"
    }else{
      ylabs<-"ratio"
    }
    if(length(table(scPACds@colData$celltype))>=7){
      message("Warning:There are more than 7 cell types,too many cell types may cause the drawing graphics to be not beautiful enough, please choose the cell type through subsetPACds()")
    }

    if(type=="box"){
      p<-ggplot(PA1, aes(x=celltype,y=expr,color=PA)) +
        geom_boxplot()+
        scale_color_manual(values=cbPalette)+
        labs(y=ylabs ,color = "PA")+
        theme_classic()
      if(!is.null(cols)){
        p<-ggplot(PA1, aes(x=celltype,y=expr,color=PA)) +
          geom_boxplot()+
          scale_color_manual(values=cols)+
          labs(y=ylabs ,color = "PA")+
          theme_classic()
      }
    }
    if(type=="violin"){
      p<- ggplot(PA1, aes(x=celltype,y=expr)) +
        geom_violin(aes(color=PA) ) +
        labs(y=ylabs ,color = "PA")+
        scale_color_manual(values=cbPalette)+
        theme_classic()
      if(!is.null(cols)){
        p<- ggplot(PA1, aes(x=celltype,y=expr)) +
          geom_violin(aes(color=PA) ) +
          labs(y=ylabs ,color = "PA")+
          scale_color_manual(values=cols)+
          theme_classic()
      }
    }
    if(type=="point"){
      p<-ggplot(PA1, aes(x=celltype,y=expr,color=PA)) +
        geom_jitter(alpha=0.3,
                    position=position_jitterdodge(jitter.width=0.2,
                                                  jitter.height=0,
                                                  dodge.width=0.8))+
        scale_color_manual(values=cbPalette)+
        labs(y=ylabs ,color = "PA")+
        theme_classic()
      if(!is.null(cols)){
        p<-ggplot(PA1, aes(x=celltype,y=expr,color=PA)) +
          geom_jitter(alpha=0.3,
                      position=position_jitterdodge(jitter.width=0.2,
                                                    jitter.height=0,
                                                    dodge.width=0.8))+
          scale_color_manual(values=cols)+
          labs(y=ylabs ,color = "PA")+
          theme_classic()
      }
    }
    return(p)

  }
}

#########################################vizBubble
library(Seurat)
# -------------- *** vizBubble() **** ---------------
#' Bubble plot
#'
#'Draws a bubble plot of single cell data (PA count,APA index)
#'
#' @usage vizBubble(PAdata,datatype,genes,celltype.anno=NULL,group.by,conds=NULL,col=c("lightgrey", "#9b393b"))
#' @param PAdata A dataset of the PACdataset class(eg:scPACds) or APA index(uses all samples to filter APA sites, and then gets index for each sample):should be a data frame
#' @param datatype PA/APA.When datatype=PA,the PA data provided is PA count;When datatype=APA,the PA data provided is APA index
#' @param genes Select the genes
#' @param celltype.anno A data frame denoting the cell classification information. Each row is one cell, and the column is the sample group
#' Row names in the celltype.anno need to match the column names of the PACdataset@counts or APA index.
#' when PAdata=PACdataset and colDataFile is provided(the information is stored in the colData at this time), celltype.anno=NULL
#' @param group.by Factor to group the cells by(the column name containing cell classification information, stored in PACdataset@colData or celltype.anno)
#' @param conds Select the type of cells needs to be plotted
#' @param col Colors to plot:the name of a palette from RColorBrewer::brewer.pal.info or a pair of colors defining a gradient
#' @return vizBubble() reutrns bubble plot
#' @examples
#' #show  PA expression data(PA count) in the gene across different cell types
#' vizBubble(PAdata=scPACds,datatype="PA",genes="ENSMUSG00000019969",group.by="celltype")
#' #change the color:the color name from RColorBrewer::brewer.pal.info
#' vizBubble(PAdata=scPACds,datatype="PA",genes="ENSMUSG00000019969",group.by="celltype",col="Set2")
#' #select the cell type
#' vizBubble(PAdata=scPACds,datatype="PA",genes="ENSMUSG00000019969",group.by="celltype",conds=c("ES","RS"),col=c("lightgrey", "green"))
#' @name vizPAplot
#' @export
vizBubble<-function(PAdata,datatype,genes,group.by,celltype.anno=NULL,conds=NULL,col=c("lightgrey", "#9b393b")){
  if(datatype=="PA"){
    PAS<-subsetPACds(PAdata,noIntergenic=TRUE,genes=genes)
    PAlist<-rownames(PAS@counts)
    Seurat<-CreateSeuratObject(PAS@counts,meta.data = PAS@colData,assay = "PA")
    if(!is.null(celltype.anno)){
      Seurat<-CreateSeuratObject(PAS@counts,meta.data = celltype.anno,assay = "PA")
    }

    Idents(Seurat)<-group.by
    if(!is.null(conds)){
      Idents(Seurat)<-conds
    }


    p<- DotPlot(Seurat,features = PAlist,cols=col,col.min=0,col.max=1)+
      theme(panel.grid.major = element_blank(), #主网格线
            panel.grid.minor = element_blank(), #次网格线
            panel.border = element_blank(), #边框
            axis.title = element_blank(),  #轴标题
            panel.background = element_rect(fill = 'white'), #背景色
            plot.background=element_rect(fill="white"))+
      coord_flip()
  }


  if(datatype=="APA"){
    PAdata<-as.data.frame(PAdata)
    genelist<-genes
    Seurat<-CreateSeuratObject(PAdata,meta.data = celltype.anno,assay = "APA") #仅考虑APA是数据框的形式
    Idents(Seurat)<-group.by
    if(!is.null(conds)){
      Idents(Seurat)<-conds
    }

    p<- DotPlot(Seurat,features = genelist,cols=col,col.min=0,col.max=1)+
      theme(panel.grid.major = element_blank(), #主网格线
            panel.grid.minor = element_blank(), #次网格线
            panel.border = element_blank(), #边框
            axis.title = element_blank(),  #轴标题
            panel.background = element_rect(fill = 'white'), #背景色
            plot.background=element_rect(fill="white"))+
      coord_flip()
  }
  return(p)
}

##########################################vizUMAP
# -------------- *** vizUMAP() **** ---------------
#' UMAP plot
#'
#' Graphs the output of a dimensional reduction technique on a 2D scatter plot where each point is a cell and it's positioned based on the cell embeddings determined by the reduction technique.
#' By default, cells are colored by their celltype (can be changed with the group.by parameter)
#'
#' @usage vizUMAP(PAdata,celltype.anno=NULL,group.by,PA = NULL,gene=NULL,dims=1:10,order=NULL,cols=NULL,order.cols=NULL,pt.cols= c("lightgrey", "#185827"),size=NULL)
#' @param PAdata A dataset of the PACdataset class(eg:scPACds) or APA index(uses all samples to filter APA sites, and then gets index for each sample):should be a data frame
#' @param celltype.anno A data frame denoting the cell classification information. Each row is one cell, and the column is the sample group
#' Row names in the celltype.anno need to match the column names of the PACdataset@counts or APA index.
#' when PAdata=PACdataset and colDataFile is provided(the information is stored in the colData at this time), celltype.anno=NULL
#' @param group.by Factor to group the cells by(the column name containing cell classification information, stored in PACdataset@colData or celltype.anno)
#' @param PA Select the poly(A) site to show it's expression in UMAP Plot
#' @param gene Select the APA gene to show the usage of PACs in UMAP Plot
#' @param dims Dimensions to plot, must be a two-length numeric vector specifying x- and y-dimensions
#' @param order Select the cell types that highlight in UMAP Plot
#' @param cols Vector of colors, each color corresponds to an identity class
#' @param order.cols The colors of the cell types that highlight(if order is not NULL)
#' @param pt.cols When PA or gene is not NULL ,the color of the point in UMAP Plot:the two colors to form the gradient over
#' @param size Adjust point size for plotting
#' @return vizUMAP() reutrns UMAP plot
#' @examples
#' vizUMAP(scPACds,group.by = "celltype",dims=1:5)
#' #Select the cell types that highlight in UMAP Plot
#' vizUMAP(scPACds,group.by = "celltype",dims=1:5,order=c("ES","RS"))
#' #Select the poly(A) site to show it's expression in UMAP Plot
#' vizUMAP(scPACds,group.by = "celltype",dims=1:5,PA="PA3346")
#' @name vizUMAP
#' @export
vizUMAP<-function(PAdata,celltype.anno=NULL,group.by,PA = NULL,gene=NULL,dims=1:10,order=NULL,cols=NULL,order.cols=NULL,pt.cols= c("lightgrey", "#185827"),size=NULL){

  palette.2 <- c("#7fb80e", "#1d953f",
                 "#2a5caa", "#a565ef",
                 "#FFA500", "#F08080",
                 "#5F9EA0", "#f9e264", "#009db2", "#375830",
                 "#FFE4B5","#d71345","#FF6347",
                 "#EE82EE","#a9ddd4","#9ec3db","#ad91cb","#7b52ae","#6b86ff")


  if(!is.null(celltype.anno)){
    Seurat<-CreateSeuratObject(PAdata,meta.data = celltype.anno,assay = "PA")
  }else{
    Seurat<-CreateSeuratObject(PAdata@counts,meta.data =PAdata@colData,assay = "PA")
  }
  ####标准化
  Seurat <- NormalizeData(Seurat, normalization.method = "LogNormalize", scale.factor = 10000,verbose = FALSE)
  Seurat<-FindVariableFeatures(Seurat,selection.method = "vst",nfeatures=100,verbose = FALSE)
  ####选择所有PA/gene进行归一化
  allPAs <- rownames(Seurat)
  umap<- ScaleData(Seurat, features = allPAs,verbose = FALSE)
  #######先PCA在UMAP,不能直接UMAP
  umap <- RunPCA(umap, features = allPAs)
  umap<-RunUMAP(object = umap,dims = dims)
  p<-DimPlot(umap,reduction ="umap" ,group.by=group.by,label=TRUE,cols= palette.2,pt.size = size)+
    theme_bw()+
    theme( panel.grid.major = element_blank(),#去掉主网格线
           panel.grid.minor = element_blank(),#去掉次网格线
           axis.text = element_blank(),#去掉坐标轴上数字
           axis.ticks = element_blank(),#去掉坐标轴刻度
           plot.title = element_blank())#去掉标题
  if(!is.null(cols)){
    p<-DimPlot(umap,reduction ="umap" ,group.by=group.by,label=TRUE,cols= cols,pt.size = size)+
      theme_bw()+
      theme( panel.grid.major = element_blank(),#去掉主网格线
             panel.grid.minor = element_blank(),#去掉次网格线
             axis.text = element_blank(),#去掉坐标轴上数字
             axis.ticks = element_blank(),#去掉坐标轴刻度
             plot.title = element_blank())#去掉标题
  }
  if(!is.null(order)){

    p<-DimPlotOrder(umap,reduction="umap",group.by = group.by,order=order,palette=palette.2,size = size)+
      theme_bw()+
      theme( panel.grid.major = element_blank(),#去掉主网格线
             panel.grid.minor = element_blank(),#去掉次网格线
             axis.text = element_blank(),#去掉坐标轴上数字
             axis.ticks = element_blank(),#去掉坐标轴刻度
             plot.title = element_blank())#去掉标题

    if(!is.null(order.cols)){
      p<- DimPlotOrder(umap,reduction="umap",group.by = group.by,order=order,palette=order.cols,size = size)+
        theme_bw()+
        theme( panel.grid.major = element_blank(),#去掉主网格线
               panel.grid.minor = element_blank(),#去掉次网格线
               axis.text = element_blank(),#去掉坐标轴上数字
               axis.ticks = element_blank(),#去掉坐标轴刻度
               plot.title = element_blank())#去掉标题
    }
  }
  if(!is.null(PA)){
    p<- FeaturePlot(umap,features=PA,reduction="umap",cols = pt.cols,pt.size = size)+
      theme_bw()+
      theme( panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             axis.text = element_blank(),
             axis.ticks = element_blank(),
             plot.title = element_text(hjust=0.5)
      )
  }
  if(!is.null(gene)){
    p<- FeaturePlot(umap,features=gene,reduction="umap",cols = pt.cols,pt.size = size)+
      theme_bw()+
      theme( panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             axis.text = element_blank(),
             axis.ticks = element_blank(),
             plot.title = element_text(hjust=0.5)
      )
  }
  return(p)
}

##################
DimPlotOrder <- function(object,reduction,group.by,order=NULL,palette=NULL,size=NULL){
  groups <- object@meta.data[,group.by]
  if(!is.null(order)){
    object@meta.data[,group.by] <- factor(groups, levels = order[sort(match(unique(groups),order))])
  }
  if(!is.null(palette)){
    color <- palette[sort(match(unique(groups),order))]
    figure <- DimPlot(object, reduction = reduction, group.by = group.by, cols = color,pt.size=size)
  }else{
    figure <- DimPlot(object, reduction = reduction, group.by = group.by,pt.size=size)
  }
  return(figure)
}



#####################subsetAPAMarkers
# -------------- *** subsetAPAMarkers() **** ---------------
#' APA markers of individual cell types
#'
#' Finds APA markers among different cell type
#'
#' @usage subsetAPAMarkers(APAindex,celltype.anno=NULL,group.by,cluster.1=NULL,cluster.2=NULL,every.cluster=FALSE,genes=NULL,only.pos=FALSE,method= "wilcox",order="down")
#' @param APAindex A data frame with rownames denoting genes, columns are samples, and values are the APA index values. Also includes other method-specific columns
#' @param celltype.anno A data frame denoting the cell classification information. Each row is one cell, and the column is the sample group
#' Row names in the celltype.anno need to match the column names of  APA index.
#' @param group.by Factor to group the cells by(the column name containing cell classification information, stored in celltype.anno)
#' @param cluster.1 Vector of cell names belonging to group 1
#' @param cluster.2 Vector of cell names belonging to group 2
#' @param every.cluster Finds APA markers for all cell types.When every.cluster=TRUE,cluster.1 and cluster.2 must be NULL
#' @param genes APA genes to test. Default is to use all APA genes
#' @param only.pos Only return positive markers (FALSE by default)
#' @param method Denotes which test to use. Available options are:
#' "wilcox" : Identifies differentially expressed APA genes between two groups of cells using a Wilcoxon Rank Sum test (default)
#' "t" : Identify differentially expressed APA genes between two groups of cells using the Student's t-test.
#' "LR" : Uses a logistic regression framework to determine differentially expressed APA genes.
#'  Constructs a logistic regression model predicting group membership based on each feature individually and compares this to a null model with a likelihood ratio test.
#' @param order When order=down,the returned values are sorted in descending order by avg_log2FC(The most significant marker in each group comes first);
#' When order=up,the returned values are sorted in ascending order by avg_log2FC(The most significant marker in each group ranks last)
#' @return subsetAPAMarkers() reutrns a data frame containing APA markers information
#' @examples
#' #find all markers of ES
#' PAmakers.1<- subsetAPAMarkers(APAindex,celltype.anno=scPACds@colData ,group.by="celltype",cluster.1 = "ES", only.pos=FALSE,method= "wilcox")
#' #find all markers distinguishing ES from RS
#' PAmakers.2<- subsetAPAMarkers(APAindex,celltype.anno=scPACds@colData ,group.by="celltype",cluster.1 = "ES", cluster.2="RS",only.pos=FALSE,method= "wilcox")
#' #find markers for every cell cluster compared to all remaining cell clusters, report only the positive ones
#' APAmarkers<- subsetAPAMarkers(APAindex,celltype.anno=scPACds@colData ,group.by="celltype",every.cluster = TRUE, only.pos=TRUE,method= "wilcox")
#' @name subsetAPAMarkers
#' @export
subsetAPAMarkers<-function(APAindex,celltype.anno=NULL,group.by,cluster.1=NULL,cluster.2=NULL,every.cluster=FALSE,
                           genes=NULL,only.pos=FALSE,method= "wilcox",order="down"){
  Seurat<-CreateSeuratObject(APAindex,meta.data =celltype.anno,assay = "APA")

  ####标准化
  Seurat <- NormalizeData(Seurat, normalization.method = "LogNormalize", scale.factor = 10000,verbose = FALSE)
  Seurat<-FindVariableFeatures(Seurat,selection.method = "vst",nfeatures=100,verbose = FALSE)
  ####选择所有PA/gene进行归一化
  allgenes <- rownames(Seurat)
  Seurat<- ScaleData(Seurat, features = allgenes,verbose = FALSE)
  #####获取操作对象
  Idents(Seurat)<-group.by

  if(!is.null(cluster.1)){
    APAmarkers <- FindMarkers( Seurat, only.pos = only.pos, min.pct = 0.25,logfc.threshold = 0.25,ident.1=cluster.1,features = genes, test.use = method)
    if(!is.null(cluster.2)){
      APAmarkers <- FindMarkers( Seurat, only.pos = only.pos, min.pct = 0.25,logfc.threshold = 0.25,

                                 ident.1=cluster.1, ident.2=cluster.2,features = genes, test.use = method)
    }
    if(order=="down"){
      APAmarkers<- APAmarkers[order(APAmarkers$avg_log2FC,decreasing = T),]
    }else{
      APAmarkers<- APAmarkers[order(APAmarkers$avg_log2FC,decreasing = F),]
    }
  }

  if(every.cluster==TRUE){
    APAmarkers <- FindAllMarkers(Seurat, only.pos = only.pos, min.pct = 0.25, logfc.threshold = 0.25,features = genes,test.use =method)
    if(order=="down"){
      APAmarkers<- APAmarkers[order(APAmarkers$cluster,-APAmarkers$avg_log2FC),]
    }else{
      APAmarkers<- APAmarkers[order(APAmarkers$cluster,APAmarkers$avg_log2FC),]
    }

  }
  return(APAmarkers)
}


##############################vizAPAMarker(test)
# -------------- *** vizAPAMarker() **** ---------------
#' Visualize APA markers in individual cell types.
#'
#' Draws a violin point/heatmap/bubble  plot to show APA markers
#'
#' @usage vizAPAMarkers(APAindex,celltype.anno,group.by,markers=NULL,type="heatmap",vln.flip=TRUE,vln.fill.by="celltype",vln.cols=NULL,dot.cols=c("lightgrey", "#9b393b"),hm.group.cols=NULL)
#' @param APAindex A data frame with rownames denoting genes, columns are samples, and values are the APA index values. Also includes other method-specific columns
#' @param celltype.anno A data frame denoting the cell classification information. Each row is one cell, and the column is the sample group
#' Row names in the celltype.anno need to match the column names of  APA index.
#' @param group.by Factor to group the cells by(the column name containing cell classification information, stored in celltype.anno)
#' @param markers APA markers to plot
#' @param type violin/heatmap/dot:specified drawing type.When type=violin/heatmap/dot,show APA markers in the form of violin plot/heatmap/bubble plot
#' @param vln.flip flip plot orientation.If type=violin and vln.flip=TRUE,cell classes on the x-axis
#' @param vln.fill.by celltype/markers:If type=violin,color violins based on either 'celltype' or 'markers'
#' @param vln.cols If type=violin,colors to use for plotting
#' @param dot.cols If type=dot,colors to use for plotting:the two colors to form the gradient over
#' @param hm.group.cols If type=heatmap,colors to use for the color bar:add a color bar showing group status for cells
#' @return vizAPAMarkers() reutrns violin plot/heatmap/bubble plot
#' @examples
#' #specified drawing type
#' topAPAM<-APAmarkers %>% group_by(cluster) %>% top_n(n =2, wt = avg_log2FC)
#' vizAPAMarkers(APAindex,celltype.anno=scPACds@colData ,group.by="celltype",type="heatmap",markers=topAPAM$gene)
#' #If type=violin,color violins based on 'celltype'
#' vizAPAMarkers(APAindex,celltype.anno=scPACds@colData ,group.by="celltype",type="violin",markers=topAPAM$gene)
#' #If type=violin,color violins based on 'markers'
#' vizAPAMarkers(APAindex,celltype.anno=scPACds@colData ,group.by="celltype",type="violin",markers=topAPAM$gene,vln.fill.by="markers")
#' @name vizAPAMarkers
#' @export
vizAPAMarkers<-function(APAindex,celltype.anno,group.by,markers=NULL,type="heatmap",
                        vln.flip=TRUE,vln.fill.by="celltype",vln.cols=NULL,dot.cols=c("lightgrey", "#9b393b"),hm.group.cols=NULL){

  Seurat<-CreateSeuratObject(APAindex,meta.data =scPACds@colData,assay = "APA")

  ####标准化
  Seurat <- NormalizeData(Seurat, normalization.method = "LogNormalize", scale.factor = 10000,verbose = FALSE)
  Seurat<-FindVariableFeatures(Seurat,selection.method = "vst",nfeatures=100,verbose = FALSE)
  ####选择所有PA/gene进行归一化
  allgenes <- rownames(Seurat)
  Seurat<- ScaleData(Seurat, features = allgenes,verbose = FALSE)
  #####获取操作对象
  Idents(Seurat)<-group.by

  palette <-c("#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887",
              "#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98",
              "#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072",
              "#7B68EE","#9400D3","#800080","#A0522D",
              "#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0")

  if(type=="violin"){

    if(vln.fill.by=="celltype"){
      vln.fill.by<-"ident"
    }else{
      vln.fill.by<-"feature"
    }
    if(is.null(vln.cols)){

      if(vln.flip==TRUE){
        p<- VlnPlot(Seurat,features=markers,pt.size=0,stack=TRUE,fill.by = vln.fill.by,
                    flip=TRUE,cols=palette
        )+
          theme(axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                legend.position = "none"
          )
      }else{
        p<- VlnPlot(Seurat,features=markers,pt.size=0,stack=TRUE,log=TRUE,fill.by = vln.fill.by,
                    cols=palette
        )+
          theme(axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                legend.position = "none"
          )
      }
    }
    if(!is.null(vln.cols)){
      if(vln.flip==TRUE){
        p<- VlnPlot(Seurat,features=markers,pt.size=0,stack=TRUE,fill.by = vln.fill.by,
                    flip=TRUE,cols=vln.cols
        )+
          theme(axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                legend.position = "none"
          )
      }else{
        p<-VlnPlot(Seurat,features=markers,pt.size=0,stack=TRUE,log=TRUE,fill.by = vln.fill.by,
                   cols=vln.cols
        )+
          theme(axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                legend.position = "none"
          )
      }

    }
  }
  if(type=="heatmap"){
    p<-DoHeatmap(Seurat,features =markers ,label = TRUE, disp.min=0,disp.max=1,group.colors=palette)
    if(!is.null(hm.group.cols)){
      p<-DoHeatmap(Seurat,features =markers ,label = TRUE, disp.min=0,disp.max=1,group.colors=hm.group.cols)
    }
  }
  if(type=="dot"){
    p<-DotPlot(Seurat,features = markers,col.min=0,col.max=1,col=dot.cols)+
      theme(panel.grid.major = element_blank(), #主网格线
            panel.grid.minor = element_blank(), #次网格线
            panel.border = element_blank(), #边框
            axis.title = element_blank(),  #轴标题
            panel.background = element_rect(fill = 'white'), #背景色
            plot.background=element_rect(fill="white"))+
      coord_flip()

  }
  return(p)
}

# -------------- *** subsetPACds() **** ---------------
#' Subset a PACdataset
#'
#' subsetPACds returns a subset of PACdataset.
#'
#' @param group a sample group, must be present in PACds@colData.
#' If NULL, then parameters of cond1, cond2, and avg are not in use.
#' @param cond1 a condition, must be present in PACds@colData.
#' @param cond2 a condition, must be present in PACds@colData.
#' If cond1 or cond2 is NULL, then both cond1 and cond2 are set NULL.
#' @param conds to subset multiple conditions.
#' @param avgPACtag if >0, then filter PACs with average PAT number across all samples >= avgPACtag, after subseting by group and conds.
#' @param avgGeneTag similar to avgPACtag, but for gene.
#' @param totPACtag filter PACs with total PAT number of all samples >=totPACtag.
#' @param noIntergenic whether to remove intergenic PACs.
#' @param avg TRUR/FALSE. If TRUE then take average of replicates of each of subset conditions.
#' @param pool TRUR/FALSE. If TRUE then take pool value of replicates of each of subset conditions.
#' @param choosePA value can e NULL/distal/proximal/APA. Only if PACds is all 3UTRAPA, then choosePA can be distal/proximal.
#'  APA means choose PACs from APA genes.
#' @param PAs to filter PACs by a PAC name list according to rownames in PAs.
#' @param genes to filter by a gene name list according to PACds@anno$gene.
#' @param chrs to filter by a chr name list according to PACds@anno$chr.
#' @param clearPAT If >0 then set @counts[<clearPAT]=0 and remove blank lines and reset PACds@anno.
#' @return a subset of PACdataset with 0-count lines removed.
#' @examples
#' data(scPACds)
#' PACds1=subsetPACds(scPACds, group='group', pool=TRUE)
#' ## Subset two conditions.
#' pacds=subsetPACds(scPACds, group='group',cond1='anther', cond2='embryo')
#' ## Subset PACs from APA genes in two conditions.
#' PACds1=subsetPACds(scPACds, group='group', cond1='anther', cond2='embryo', choosePA='apa')
#' ## Subset PACs of given genes.
#' subsetPACds(scPACds, genes=scPACds@anno$gene[1:5], verbose=TRUE)
#' ## Subset PACs of given PA names.
#' swPAC=subsetPACds(scPACds, PAs=rownames(scPACds@counts)[1:5], verbose=TRUE)

#' @name subsetPACds
#' @family PACdataset functions
#' @seealso [samplePACds()] to get sampled PACds.
#' @export
subsetPACds<-function(pacds,
                      group=NULL, cond1=NULL, cond2=NULL, conds=NULL,
                      avgPACtag=0, avgGeneTag=0, totPACtag=0,
                      choosePA=NULL, PAs=NULL, genes=NULL, chrs=NULL,
                      noIntergenic=FALSE, avg=FALSE, pool=FALSE,
                      clearPAT=0,
                      verbose=FALSE) {

  if (avg & pool) stop('avg or pool')
  n1=nrow(pacds@counts)
  txt=c('before subsetPACds'=n1)

  if (is.null(group)) {
    cond1=NULL; cond2=NULL; avg=FALSE; conds=NULL

    if(avg | pool) stop("group is NULL but avg or pool is true, don't know how to summarize")
  }

  if (!is.null(choosePA)) {
    choosePA=tolower(choosePA)
    if (!(choosePA %in% c('distal','proximal','apa'))) {
      stop("choosePA must be distal/proximal/APA")
    }
  }

  if (noIntergenic) {
    ii=getNonItgFtrId(pacds@anno$ftr)
    pacds=pacds[ii]
    txt=c(txt,'noItg'=nrow(pacds@counts))
  }

  if (is.null(conds)) {
    if (!is.null(group) & (is.null(cond1) | is.null(cond2))) {
      conds=levels(pacds@colData[,group])
    } else if (!is.null(group)) {
      conds=c(cond1, cond2)
    }
  }


  if (length(conds)>0) {
    if (!AinB(conds, pacds@colData[,group]) ) {
      stop(cat(paste(conds[1], collapse = ','),'not all in pacds@colData\n'))
    }

    smps=rownames(pacds@colData)[which(pacds@colData[,group] %in% conds)]
    pacds@counts=pacds@counts[,smps,drop=F]
    pacds@colData=pacds@colData[smps,,drop=F]
    for (i in 1:ncol(pacds@colData)) {
      if (is.factor(pacds@colData[,i])) pacds@colData[,i]=droplevels(pacds@colData[,i])
    }
    row0=which(rowSums(pacds@counts)==0)
    if (length(row0)>0) {
      pacds=pacds[-row0]
    }
    txt=c(txt,'After filter conds'=nrow(pacds@counts))
  }

  if (avgPACtag>0) {
    pacds=pacds[rowMeans(pacds@counts)>=avgPACtag]
    n2=nrow(pacds@counts)
    txt=c(txt,nrow(pacds@counts))
    names(txt)[length(txt)]=paste0('avgPACtag>=',avgPACtag)
  }

  if (totPACtag>0) {
    tot=rowSums(pacds@counts)
    pacds=pacds[tot>=totPACtag]
    txt=c(txt,nrow(pacds@counts))
    names(txt)[length(txt)]=paste0('totPACtag>=',totPACtag)
  }

  if (avgGeneTag>0) {
    pacds@anno=pacds@anno[rownames(pacds@counts),,drop=F]
    rs=rowSums(pacds@counts)
    b=aggregate(rs, list(gene=pacds@anno$gene), sum)
    gs=b$gene[b$x/ncol(pacds@counts)>=avgGeneTag]
    pacds=pacds[pacds@anno$gene %in% gs]
    txt=c(txt,nrow(pacds@counts))
    names(txt)[length(txt)]=paste0('avgGeneTag>=',avgGeneTag)
  }

  if (!is.null(PAs)) {
    rn=rownames(pacds@counts)
    pacds=pacds[rn[rn %in% PAs]]
    txt=c(txt,nrow(pacds@counts))
    names(txt)[length(txt)]='PAs'
  }

  if (!is.null(genes)) {
    gcol=NA
    if ('gene_id' %in% colnames(pacds@anno)) gcol='gene_id'
    if ('gene' %in% colnames(pacds@anno)) gcol='gene'
    if (is.na(gcol)) stop("subsetPACds by gene, but gene_id or gene not in pacds@anno\n")
    pacds=pacds[pacds@anno[, gcol] %in% genes]
    txt=c(txt,nrow(pacds@counts))
    names(txt)[length(txt)]='genes'
  }


  if (!is.null(chrs)) {
    pacds=pacds[pacds@anno$chr %in% chrs]
    txt=c(txt,nrow(pacds@counts))
    names(txt)[length(txt)]='chrs'
  }

  if (avg | pool) {
    if (avg) txt=c(txt,"averaging")
    if (pool) txt=c(txt,"pool")
    d=matrix()
    conds=unique(pacds@colData[,group])
    for (cond in conds) {
      if (avg) {
        di=rowMeans(pacds@counts[,rownames(pacds@colData)[which(pacds@colData[,group]==cond)], drop=F])
      } else {
        di=rowSums(pacds@counts[,rownames(pacds@colData)[which(pacds@colData[,group]==cond)], drop=F])
      }
      if (conds[1]==cond) {
        d=di
      } else {
        d=cbind(d,di)
      }
    }
    colnames(d)=conds
    rownames(d)=rownames(pacds@counts)
    pacds@counts=as.data.frame(d)
    pacds@colData=data.frame(group=conds)
    colnames(pacds@colData)=group
    rownames(pacds@colData)=conds
    pacds@counts=round(pacds@counts)
    for (i in 1:ncol(pacds@colData)) {
      pacds@colData[,i]=factor(pacds@colData[,i])
    }
  }

  #Note the filtering order!
  #If avg/pool=T, then do that first and then clearPAT
  #if a line is all 0, then remove the PAC.
  if (!is.null(clearPAT)) {
    if (clearPAT>0) {
      pacds@counts[pacds@counts<clearPAT]=0
      idx=rowSums(pacds@counts)!=0
      pacds=pacds[idx]
      txt=c(txt,'After clearPAT'=nrow(pacds@counts))
    }
  }

  idx=rowSums(pacds@counts)!=0
  if (sum(idx)!=length(pacds)) {
    pacds=pacds[idx]
    txt=c(txt,'After removing 0 lines'=nrow(pacds@counts))
  }

  #After clearPAT, then choosePA, because clearPAT will remove all-0 lines.
  if (!is.null(choosePA)) {
    if (choosePA %in% c('distal','proximal')) {
      pacds=get3UTRAPAds(pacds, sortPA=TRUE, choose2PA='PD')
      if (choosePA=='distal') {
        pacds=pacds[seq(2,nrow(pacds@counts),2)]
      } else {
        pacds=pacds[seq(1,nrow(pacds@counts),2)]
      }
    } else if (choosePA=='apa') {
      genes1=unique(pacds@anno$gene[duplicated(pacds@anno$gene)])
      pacds=pacds[pacds@anno$gene %in% genes1]
    }

    txt=c(txt,nrow(pacds@counts))
    names(txt)[length(txt)]=choosePA
  }


  if (verbose) {
    txt=cbind(txt)
    print(txt)
  }

  return(pacds)
}

# Given the ftr vectorget id that are not ^inter
# - ftr: character vector
# - return: idx of non-itg ftrs
getNonItgFtrId<-function(ftr) {
  id=grep('^inter',ftr)
  if (length(id)>0) {
    return((1:length(ftr))[-id])
  } else {
    return(1:length(ftr))
  }
}
