#!/usr/bin/env Rscript
# Nov 2, 2016
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


# modified from gsea, short version with wilcox.mann-whitney u test, object-based
runGSEA_short <- function(objects_ordered_df, genesets, p_adjust_method='fdr', pval_cutoff=0.05){
  gs_list <- genesets$list
  # ob_score <- objects_ordered_df[, c('object_id', 'score')]
  
  score <- -log(objects_ordered_df$pval)
  # score  <- dchisq(objects_ordered_df$pval, df=1)
  score[is.infinite(score)] <- max(score[is.finite(score)])+1
  ob_score <- data.frame(object_id=objects_ordered_df$object_id, score)
  
  gsea_short_outlist <- mclapply(gs_list, function(x){
    inMap <- ob_score$object_id%in%x
    inScore <- ob_score[inMap, 'score']
    outScore <- ob_score[!inMap, 'score']
    gsea_df <- list()
    # gsea_df$gs_size <- length(x)
    # gsea_df$num <- length(inScore)
    gsea_df$statistic <- 0
    gsea_df$pval <- 1
    if (length(inScore) > 0){
      test_out <- wilcox.test(inScore, outScore, alternative='greater')
      gsea_df$statistic <- test_out$statistic[[1]]
      gsea_df$pval <- test_out$p.value
    }
    return(gsea_df)
  })
  
  gsea_short_df <- as.data.frame(do.call(rbind, gsea_short_outlist))
  gsea_short_df$imid <- names(gs_list)
  gsea_short_df$mapname <- genesets$id2names[match(gsea_short_df$imid, genesets$id2names$imid), .(mapname)]
  gsea_short_df$qval <- unlist(p.adjust(gsea_short_df$pval, method=p_adjust_method))
  gsea_short_passed <-gsea_short_df[gsea_short_df$pval<=pval_cutoff,c('imid','mapname','pval', 'qval')]
  return(gsea_short_passed)
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
# parser$add_argument("-n", "--obnum", type="integer", nargs=1, help="num_of_objects_to_run_test", required=T)
parser$add_argument("-m", "--padj", type="character", nargs=1, help="p_adjust_method", default='fdr')
parser$add_argument("-p", "--pval", type="double", nargs=1, help="pval_cutoff", default=0.05)
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

# cat('Loading input files...\n')
genesets <- prepareGS(in_dir, gs_file)

g2id <- read.table(paste(in_dir, g2id_file, sep=""), header=TRUE)  
ob2id <- read.table(paste(in_dir, ob2id_file, sep=""), header=TRUE)
g2ob <- merge(g2id, ob2id, by='entrez_id')


time1 <- proc.time()
# objects_ordered_df <- prepareInput(in_dir, in_file, g2id_file, ob2id_file)
objects_ordered_df <- prepareInput2(in_dir, in_file, g2ob)
object_list <- objects_ordered_df$object_id


# Run test

# cat('Start test...\n')
gsea_short_passed <- runGSEA_short(objects_ordered_df, genesets, p_adjust_method=args$padj, pval_cutoff=args$pval)
top_ea <- gsea_short_passed[order(gsea_short_passed$qval),]

toJSON(top_ea)
# write.table(top_hg, file=paste(outdir, '/hgout.txt', sep=''), sep='\t', quote=FALSE, row.names=FALSE)


# proc.time()-time1
# cat('Done\n')
quit('no')