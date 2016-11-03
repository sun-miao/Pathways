#!/usr/bin/env Rscript
# Oct 31, 2016
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
  data <- read.table(paste(in_dir, in_file, sep=""), sep=',', header=TRUE)
  # data1 <- data[order(data$max.statistic., decreasing=TRUE),]
  # g2ob$score <- data1$max.statistic.[match(g2ob$gene_name, data1$gene_name)]
  # g2ob$score[is.na(g2ob$score)] <- 0
  
  # just use the p-values
  data2 <- data[order(data$min.pValue.),]
  g2ob$pval <- data2$min.pValue.[match(g2ob$gene_name, data2$gene_name)]
  g2ob$pval[is.na(g2ob$pval)] <- 1
  gene_ob_df1 <- g2ob[order(g2ob$pval),]
  return(gene_ob_df1)
}


prepareIG <- function(in_dir, inter_file, objects_ordered_df, include_low=TRUE){
  all <- fread(paste(in_dir, inter_file, sep=""), sep=",", header=TRUE)
  if (include_low==TRUE){
    el_dt <- all[(all$trust!='no link')&(all$object_id1%in%objects_ordered_df$object_id)&
                 (all$object_id2%in%objects_ordered_df$object_id)&
                 (all$object_id1!=all$object_id2),
                 .(object_id1, object_id2, link_id, mechanism, effect, trust)]
  } else {
    el_dt <- all[(all$trust%in%c('high','medium'))&(all$object_id1%in%objects_ordered_df$object_id)&
                 (all$object_id2%in%objects_ordered_df$object_id)&
                 (all$object_id1!=all$object_id2),
                 .(object_id1, object_id2, link_id, mechanism, effect, trust)]
  }

  for (col in c("link_id", "object_id1", "object_id2")) el_dt[, (col) := as.character(el_dt[[col]])]
  G <- graph.data.frame(el_dt, directed=TRUE) # this will assign all the edge attributes
  # G <- graph.edgelist(apply(as.matrix(el_dt[,.(object_id1, object_id2)]),1:2,as.character), directed=TRUE)

  # V(G)$score <- objects_ordered_df[match(names(V(G)),as.character(objects_ordered_df$object_id)), 'score']
  V(G)$pval <- objects_ordered_df[match(names(V(G)),as.character(objects_ordered_df$object_id)), 'pval']
  V(G)$gene_name <- as.character(as.factor(objects_ordered_df[match(names(V(G)),as.character(objects_ordered_df$object_id)), 'gene_name']))
  
  
  
  return(G)
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
parser$add_argument("-n", "--obnum", type="integer", nargs=1, help="num_of_objects_to_run_test")
# parser$add_argument("-l", "--genefile", type="character", nargs=1, help="user_defined_gene_list.csv")
args <- parser$parse_args()

in_dir <- ''
in_file <- 'DM_comparative_FeEu_gsaTest1.csv'
g2id_file <- 'prejoint_gene_names_to_entrez_id_06172016.tsv'
ob2id_file <- 'object_to_entrez_human.tsv'

# in_file <- args$infile
# g2id_file <- args$g2idfile
# ob2id_file <- args$ob2idfile
# in_dir <- '/Users/msun/Documents/Working/Pathways'

# gs_file <- 'pathway_maps_objects_08252016.csv'
inter_file <- 'interactions_all_06302016.csv'

# Prepare input files
# cat('Loading input files...\n')
# genesets <- prepareGS(in_dir, gs_file)

g2id <- read.table(paste(in_dir, g2id_file, sep=""), header=TRUE)  
ob2id <- read.table(paste(in_dir, ob2id_file, sep=""), header=TRUE)
g2ob <- merge(g2id, ob2id, by='entrez_id')

time1 <- proc.time()
# objects_ordered_df <- prepareInput(in_dir, in_file, g2id_file, ob2id_file)
objects_ordered_df <- prepareInput2(in_dir, in_file, g2ob)
object_list <- objects_ordered_df$object_id

# visualize selected objects in interactome
G <- prepareIG(in_dir, inter_file, objects_ordered_df, include_low=FALSE)

# if given obnum
subG <- dNetInduce(G, objects_ordered_df[1:args$obnum,'object_id'], knn = 0, remove.loops = F, largest.comp = F)
subG_df <- as_long_data_frame(subG)[,3:12]
toJSON(subG_df, digit=NA)

# if given custom list of entrez_id
# entrez_gene_list <- read.table(paste(in_dir, args$genefile, sep=""))[[1]]
# subG <- dNetInduce(G, objects_ordered_df[objects_ordered_df$entrez_id%in%entrez_gene_list,'object_id'],
#                    knn = 0, remove.loops = F, largest.comp = F)
# subG_df <- as_long_data_frame(subG)[,3:12]

# proc.time()-time1
# cat('Done\n')
quit('no')