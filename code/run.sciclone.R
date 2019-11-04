# ---------------------------------------------------
# Title: Clonal Evolution Analysis
#        for GBM Case Study (Part II)
# Author: Paula Restrepo, paularstrpo@gmail.com
# Last updated: 11/04/2019
# ---------------------------------------------------

# ---------------------------------------------------
# 1) Run sciclone without CNV information
# ---------------------------------------------------

library(sciClone)
load('data/clonality.inputs.rdata')

vaf_table <- do.call(rbind, vaf_list)


ggplot(data=vaf_table) + 
  aes(x=as.factor(chromNum), y=VAF, color=as.factor(chromNum)) + 
  geom_violin() + geom_jitter(alpha=0.5) + 
  theme_classic()

vaf_list_clean <- lapply(vaf_list, function(x){
    x[! (as.factor(x$chromNum) %in% c(1,7,9,12,14)),]})


# all samples
sc2 <- sciClone(vafs=vaf_list_clean,
               copyNumberCalls=cnv_list,
               cnCallsAreLog2=TRUE,
               sampleNames=names,
               minimumDepth = 300,
               useSexChrs=FALSE)


sc.variants <- data.frame(sc@vafs.merged)
#save(sc, sc.variants, file='_rdata/sciclone.results.rdata')


## output results table and cluster plots
writeClusterTable(sc, "_reports/tables/sciclone/sciclone.clusters.tsv")

sc.plot1d(sco=sc,
          outputFile="_reports/figures/sciclone.clusters.1d.pdf", cnToPlot=c(2))

sc.plot2d(sco=sc,
         outputFile="_reports/figures/sciclone.clusters.2d.pdf", singlePage = TRUE, ylim = 40, xlim=40)

