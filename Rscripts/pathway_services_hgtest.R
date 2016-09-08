#!/usr/bin/env Rscript
# Sept 1, 2016
# msun
# pathway/interaction analysis modules

# Load libraries

suppressPackageStartupMessages(library("igraph"))
suppressPackageStartupMessages(library("dnet"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("jsonlite"))
options(mc.cores=detectCores())

# Functions used
prepareInput <- function(in_dir,in_file, g2id_file, ob2id_file){
  # input data file
  data <- read.table(paste(in_dir, in_file, sep=""), sep=',', header=TRUE)[,1:3]  
  g2id <- read.table(paste(in_dir, g2id_file, sep=""), header=TRUE)  
  ob2id <- read.table(paste(in_dir, ob2id_file, sep=""), header=TRUE)
  
  # joins gene_name, entrez_id, score, pval, object_id
  data1 <- data[order(data$max.statistic., decreasing=TRUE),]
  g2id$score <- data1$max.statistic.[match(g2id$gene_name, data1$gene_name)]
  g2id$score[is.na(g2id$score)] <- 0 
  data2 <- data[order(data$min.pValue.),]
  g2id$pval <- data2$min.pValue.[match(g2id$gene_name, data2$gene_name)]
  g2id$pval[is.na(g2id$pval)] <- 1 
  gene_ob_df <- merge(g2id, ob2id, by='entrez_id')
  gene_ob_df1 <- gene_ob_df[order(gene_ob_df$score, decreasing=TRUE),]
  return(gene_ob_df1)
}


prepareInput2 <- function(in_dir,in_file, g2ob){
  # input data file
  data <- read.table(paste(in_dir, in_file, sep=""), sep=',', header=TRUE)[,1:3]  
  data1 <- data[order(data$max.statistic., decreasing=TRUE),]
  g2ob$score <- data1$max.statistic.[match(g2ob$gene_name, data1$gene_name)]
  g2ob$score[is.na(g2ob$score)] <- 0 
  data2 <- data[order(data$min.pValue.),]
  g2ob$pval <- data2$min.pValue.[match(g2ob$gene_name, data2$gene_name)]
  g2ob$pval[is.na(g2ob$pval)] <- 1 
  gene_ob_df1 <- g2ob[order(g2ob$score, decreasing=TRUE),]
  return(gene_ob_df1)
}


prepareGS <- function(in_dir, gs_file){
  gs <- fread(paste(in_dir, gs_file, sep=""), sep='|')
  gs_grouped <- gs[,.(ob_set=list(unique(object_id))), by=imid]
  gs_list <- as.list(gs_grouped[['ob_set']])
  names(gs_list) <- as.character(as.factor(gs_grouped[,.(imid)][[1]]))
  id2names <- unique(gs[,.(imid, mapname)])
  gs_maps <- list()
  gs_maps$list <- gs_list
  gs_maps$id2names <- id2names
  gs_maps$objects <- gs$object_id
  return(gs_maps)
}

run_phyper2 <- function(objects_ordered_df, genesets, ltype=c('entrez_id','gene_name','object_id'), 
                        ob_num=0, genelist=NULL, p_adjust_method='fdr', min_overlap=3){
  gs_list <- genesets$list
  gs_objects <- genesets$objects
  
  if (ob_num>0) {
    objects_ordered_df1 <- objects_ordered_df[!duplicated(objects_ordered_df$object_id),]
    ob_drawn <- objects_ordered_df1[1:ob_num, 'object_id']
  } else if (!is.null(genelist)){
    if (ltype=='entrez_id'){
      objects_ordered_df1 <- objects_ordered_df[!duplicated(objects_ordered_df$entrez_id),]
      ob_drawn <- objects_ordered_df1[objects_ordered_df1$entrez_id%in%genelist, 'object_id']
    } else if (ltype=='gene_name'){
      objects_ordered_df1 <- objects_ordered_df[!duplicated(objects_ordered_df$gene_name),]
      ob_drawn <- objects_ordered_df1[objects_ordered_df1$gene_name%in%genelist, 'object_id']
    } else if (ltype=='object_id'){
      objects_ordered_df1 <- objects_ordered_df[!duplicated(objects_ordered_df$object_id),]
      ob_drawn <- genelist[genelist%in%objects_ordered_df1$object_id]
    } else print('please specify ltype: entrez_id|gene_name|object_id')
  } else return()
  
  gs_hyper <- mclapply(gs_list, function(x) {
    x_in <- x%in%ob_drawn
    hg_out <- list()
    hg_out$hyper.p <- 1
    hg_out$num <- 0
    hg_out$overlap <- NA
    len_x <- sum(x_in)
    len_m <- length(x)
    len_n <- sum(!gs_objects%in%x)
    len_k <- length(ob_drawn)
    hyper_p <- phyper(len_x-1, len_m, len_n, len_k, lower.tail = FALSE)
    overlapped_objects <- x[x_in]
    if (len_x > 0){
      hg_out$hyper.p <- hyper_p
      hg_out$num <- len_x
      hg_out$overlap <- paste(overlapped_objects, collapse=',')
    }
    return(hg_out)
  })  
  
  gs_hyper_df <- as.data.frame(do.call(rbind, gs_hyper))
  gs_hyper_df$imid <- as.numeric(names(gs_list))
  gs_hyper_df$mapname <- unlist(genesets$id2names[match(gs_hyper_df$imid, genesets$id2names$imid), .(mapname)])
  gs_hyper_df$hyper.p <- unlist(gs_hyper_df$hyper.p)
  gs_hyper_df$num <- unlist(gs_hyper_df$num)
  gs_hyper_df$qval <- unlist(p.adjust(gs_hyper_df$hyper.p, method=p_adjust_method))
  gs_hyper_df$overlap <- as.character(gs_hyper_df$overlap)
  
  gs_hyper_passed <- gs_hyper_df[gs_hyper_df$num>=min_overlap,c('imid','mapname','hyper.p', 'qval', 'num', 'overlap')]
  return(gs_hyper_passed)
}

# Load arguments

parser <- ArgumentParser(description='Pathway-Interactome analysis modules')

# this should be pulled from the gene/protein tab
# parser$add_argument("-i", "--infile", type="character", nargs=1, help="gene_name_score_pval.csv", required=T)

# this could came from annotation table - refseq_genes
# parser$add_argument("-g", "--g2idfile", type="character", nargs=1, help="gene_name_to_entrez_id.tsv", required=T)

# this came from TR metabase tables, currently on themata - object_to_entrez_human
# parser$add_argument("-o", "--ob2idfile", type="character", nargs=1, lesshelp="object_id_to_entrez_id.tsv", required=T)

# user-input
parser$add_argument("-n", "--obnum", type="integer", nargs=1, help="num_of_objects_to_run_test", required=T)
parser$add_argument("-m", "--padj", type="character", nargs=1, help="p_adjust_method", default='fdr')
args <- parser$parse_args()

in_dir <- ''
in_file <- 'DM_comparative_FeEu_gsaTest1.csv'
g2id_file <- 'prejoint_gene_names_to_entrez_id_06172016.tsv'
ob2id_file <- 'object_to_entrez_human.tsv'

# in_file <- args$infile
# g2id_file <- args$g2idfile
# ob2id_file <- args$ob2idfile
# in_dir <- '/Users/msun/Documents/Working/Pathways'

gs_file <- 'pathway_maps_objects_08252016.csv'

# outdir <- paste(in_dir, '/pw_service_test/', strftime(Sys.time(),"%Y-%m-%d_%H-%M-%S"),sep='')
# dir.create(outdir, showWarnings = F, recursive = T, mode = "0755")

# sink(paste(outdir,'/parameters.txt',sep=''), append=FALSE,split=FALSE)
# params = sort(names(args));
# for(i in 1:length(params)){
#   cat(paste(params[i],": ",args[[params[i]]],"\n",sep=""));
# }
# sink();

# Prepare input files

cat('Loading input files...\n')
genesets <- prepareGS(in_dir, gs_file)

g2id <- read.table(paste(in_dir, g2id_file, sep=""), header=TRUE)  
ob2id <- read.table(paste(in_dir, ob2id_file, sep=""), header=TRUE)
g2ob <- merge(g2id, ob2id, by='entrez_id')


time1 <- proc.time()
# objects_ordered_df <- prepareInput(in_dir, in_file, g2id_file, ob2id_file)
objects_ordered_df <- prepareInput2(in_dir, in_file, g2ob)
object_list <- objects_ordered_df$object_id


# Run test

cat('Start test...\n')
gs_hyper_passed <- run_phyper2(objects_ordered_df, genesets, ltype='object_id', ob_num=args$obnum, p_adjust_method=args$padj, min_overlap=3)
# gs_hyper_passed <- run_phyper2(objects_ordered_df, genesets, ltype='object_id', ob_num=200, p_adjust_method='fdr', min_overlap=3)

top_hg <- gs_hyper_passed[order(gs_hyper_passed$qval),]
toJSON(top_hg)
# write.table(top_hg, file=paste(outdir, '/hgout.txt', sep=''), sep='\t', quote=FALSE, row.names=FALSE)


# toJSON(top_hg)
proc.time()-time1
cat('Done\n')
quit('no')