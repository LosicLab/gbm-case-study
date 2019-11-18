# ---------------------------------------------------
# Title: Clonal Evolution Analysis
#        for GBM Case Study (Part I)
# Author: Paula Restrepo, paularstrpo@gmail.com
# Last updated: 02/01/2019
# ---------------------------------------------------

libs <- c("tidyverse", "clonevol","sciClone", 'fishplot')
lapply(libs, require, character.only=TRUE)
rm(libs)

load('_rdata/cosmic.cgc.rdata')
load('_rdata/parsed.variants.dnaseq.rdata')
load('_rdata/cnv.segments.rdata')
#load('_rdata/master.tables.wxs.rdata')

# ---------------------------------------------------
# 1) prep inputs for bamreadcount which
#    will then get the proper sciclone input lists.
# ---------------------------------------------------
names <- levels(factor(cns_table$sample))


sample.names <- c("10647A_vs_10647B", "10647A_vs_9260B_Blood", "10647A_vs_9260B", "10647B_vs_10647A", "10647B_vs_10647C", "10647B_vs_9260B_Blood", "10647B_vs_9260B", "10647C_vs_10647A", "10647C_vs_10647B", "10647C_vs_9260B_Blood", "10647C_vs_9260B", "9260B_vs_9260B_Blood" )


somatic.samples <- c("10647A_vs_9260B_Blood","10647B_vs_9260B_Blood","10647C_vs_9260B_Blood" ,"9260B_vs_9260B_Blood" )
vcf_filtered_df$sample <- as.factor(vcf_filtered_df$sample)
levels(vcf_filtered_df$sample) <- sample.names

master_som_filtered <- vcf_filtered_df %>% as_tibble() %>% filter(sample %in% somatic.samples)

query.bed <- master_som_filtered %>% as_tibble() %>%
             mutate(chr=`CHROM`) %>%
             select(`chr`, `POS`) %>%
             mutate(key=paste(`chr`, `POS`, sep='_')) %>%
             filter(! duplicated(`key`)) %>%
             select(- `key`) %>%
             mutate(POS2= `POS`)

write_tsv(query.bed, path='_reports/tables/sciclone/readCountQuerySites.txt', col_names = FALSE)


query.bed2 <- master_som_filtered %>% as_tibble() %>%
    mutate(chr=`CHROM`) %>%
    select(`chr`, `POS`, `REF`, `ALT`) %>%
    mutate(variant=paste(`chr`, `POS`, `REF`, `ALT`,sep='_')) %>%
    filter(! duplicated(`variant`)) %>%
    mutate(POS2= `POS`)

# ---------------------------------------------------
# 2) Use the resulting query list of mutation sites
#    from all samples to pull all the read counts at
#    those sites
# ---------------------------------------------------

# ***** use bam-readcount - results in minerva:
# * results: /sc/orga/projects/losicb01a/common_folder/restrp01/GBM/raw_data/wxs/hg38_gbm-case/bamreadcount/
# * script: /sc/orga/projects/losicb01a/common_folder/restrp01/GBM/raw_data/wxs/hg38_gbm-case/_scripts/run_bamreadcount.sh

# ---------------------------------------------------
# 3) read in the results of bam-readcount and parse
#    them into something sciclone can use.
# ---------------------------------------------------
path='~/Minerva/losicb01a/common_folder/restrp01/GBM/raw_data/wxs/hg38_gbm-case/readcounts'
flist <- list.files(path=path, pattern='*.txt$', full.names = TRUE)

# from the docs for bam-readcount at https://github.com/genome/bam-readcount
info_col <- 'base:count:avg_mapping_quality:avg_basequality:avg_se_mapping_quality:num_plus_strand:num_minus_strand:avg_pos_as_fraction:avg_num_mismatches_as_fraction:avg_sum_mismatch_qualities:num_q2_containing_reads:avg_distance_to_q2_start_in_q2_reads:avg_clipped_length:avg_distance_to_effective_3p_end'
info_col <- strsplit(x=info_col, split=':', fixed=TRUE)
info_col <- as.data.frame(t(do.call(rbind, info_col)))
info_col <- as.character(info_col$V1)

readCounts <- lapply(flist, function(x){
    x %>%
    read_tsv(col_names=c('chr','position', 'reference_base', 'depth', paste('info', 1:5, sep=''))) %>%
    as_tibble() %>%
    gather(key, info, paste('info', 1:5, sep='')) %>%
    filter(key != 'info1') %>%
    select(- key) %>%
    separate(info, info_col, sep = ':', convert = TRUE) %>%
    mutate(VAF = (`count` / `depth`) * 100,
           variant=paste(chr, position, reference_base, base, sep='_'))
})

vaf_list <- lapply(readCounts, function(x){
    x %>%
    mutate(chromNum = as.numeric(as.character(gsub(pattern='chr', replacement='', x=`chr`))),
           ref_reads = as.numeric(`depth`) - as.numeric(`count`),
           var_reads = as.numeric(`count`)) %>%
    filter(! is.na(`chromNum`)) %>% inner_join(query.bed2, by=c('variant')) %>%
        select(`chromNum`, `position`, `ref_reads`, `var_reads`, `VAF`) %>% as.data.frame()
})

# purity correction
vaf_list2 <- lapply(vaf_list, function(x){
    if(mean(x$VAF) <= 30){
    purity <-
    x <- x %>% mutate(VAF = VAF * 2 )
    }else{
        return(x)
    }
})


cnv_list <- lapply(cns_list, function(x){
    as_tibble(x) %>%
    mutate(chromosome=as.numeric(gsub(pattern='chr', replacement='', x=`chromosome`))) %>%
    select(`chromosome`, `start`, `end`, `log2`) %>%
    filter(! is.na(`chromosome`)) %>% as.data.frame()
    } )

# save everything you need to run sciclone.
# and maybe a little more just in case...
# It will take a *while* to run.
save(vaf_list, cnv_list, readCounts, master_som_filtered, names, file='_rdata/clonality.inputs.rdata')